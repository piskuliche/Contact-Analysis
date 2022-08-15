! Compile with  gfortran contact.f90 -I $lb_gmx_inc -L $lb_gmx_lib -lgmxfort -Ofast -fopenmp -o bin/contact 
! Contact: Zeke Piskulich piskulichz@gmail.com
! Copyright Aug 2022, Zeke Piskulich Boston University
! Part of git@github.com:piskuliche/Contact-Analysis.git

module constants
  implicit none
  real,parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: ifile
  integer :: nframes,n_calc_frames
  integer :: fstart, fstop
  integer :: frame
  integer :: natoms, n
  integer :: ncorr,nsep
  integer :: num_g1,num_g2
  integer :: count_g1, count_g2
  real    :: r_cut, r_cut_sq, rz_thick
  
  real, dimension(3) :: L
  real, dimension(3,3) :: dims

  character(len=40) :: fname,iname
  
  common ifile
  common nframes,n_calc_frames
  common fstart,fstop
  common frame
  common natoms,n
  common ncorr,nsep
  common num_g1, num_g2
  common count_g1, count_g2
  common r_cut, r_cut_sq, rz_thick


  common L
  common dims
  common fname, iname

end module parameters

program contacts
  use gmxfort_trajectory
  use constants
  use parameters
  implicit none
  type (Trajectory) :: trj
  integer :: test,i,j
  integer :: cP,cO,cH
  integer :: nstop
  real    :: edote
  real    :: start,finish

  integer, allocatable :: is_g1(:), is_g2(:)
  


  ! Read Input File
  call Read_Input()
  write(*,*) trim(fname),trim(iname)

  

  call trj%open(trim(fname),trim(iname))
  nframes = trj%nframes
  natoms = trj%natoms()
  write(*,*) "There are ", nframes, " frames"
  write(*,*) "There are ", natoms, " atoms"

  ! Read is_atoms.txt and assign the is_g1 and is_g2 vectors
  ! These vectors store 1 if it is a match
  ! and stores a 0 if it is not. 
  allocate(is_g1(natoms)); allocate(is_g2(natoms))
  call Find_Atoms(is_g1,is_g2)
  

  ! Read Frame and call contacts
  call cpu_time(start)
  open(21,file='contacts.out')
  do i=1,n_calc_frames
    if (mod(i,500) == 0)  write(*,*) "Reached frame ",i
    call Read_Frame(is_g1,is_g2,trj)
  enddo
  call cpu_time(finish)
  print '("Time = ", f6.3, " seconds.")', finish-start
  call trj%close()
end program contacts

subroutine Find_Atoms(is_g1,is_g2)
  use parameters
  use constants
  implicit none

  integer :: i
  integer,dimension(natoms) :: is_g1, is_g2
  
  is_g1 = 0; is_g2=0
  num_g1=0; num_g2=0
  open(11,file="is_atom.txt",status='old')
  do i=1, natoms
    read(11,*) is_g1(i), is_g2(i)
  enddo
  close(11)
  num_g1 = sum(is_g1); num_g2 = sum(is_g2)
end subroutine

subroutine Read_Frame(is_g1,is_g2,trj)
  use gmxfort_trajectory
  use parameters
  use constants
  implicit none

  type (Trajectory) :: trj

  integer :: i,j

  real ::  maxz, minz

  integer, dimension(natoms) :: is_g1, is_g2
  integer, dimension(num_g1) :: g1_id
  real(8), dimension(3) :: coord
  real(8), dimension(num_g1,3) :: r_g1
  real(8), dimension(num_g2,3) :: r_g2
  real(8), dimension(num_g2) :: zs
  real(8), dimension(3,3) :: mybox
  
  ! Zero vectors
  r_g1 = 0.0; r_g2 = 0.0
  mybox=0.0; g1_id = 0

  ! Read frame
  n = trj%read_next(1)

  ! Loop over atoms and assign coordinates
  ! for the group that is the reference group
  count_g1 = 0; count_g2 = 0
  coord=0; zs=0.0
  do j=1, natoms
    coord = trj%x(1,j)
    if (is_g2(j) == 1) then
      r_g2(count_g2+1,:) = coord
      count_g2 = count_g2 + 1
      zs(count_g2) = coord(3)
    endif
  enddo

  maxz = maxval(zs)
  minz = minval(zs)
  write(*,*) "Maximum Value of Z for G2: ", maxz
  write(*,*) "Minimum Value of Z for G2: ", minz

  ! Loop over the 
  do j=1, natoms
    if (is_g1(j) == 1) then
      ! Read coordinates
      coord = trj%x(1,j)
      ! Identifies whether z coordinate is within the shell of interest (5 angstroms above/below)
      ! This cuts down on the distance calculation required!
      if (coord(3) < maxz+.5 .and. coord(3) > minz-.5) then
        r_g1(count_g1+1,:) = coord
        count_g1 = count_g1 + 1
        g1_id(count_g1) = j 
      endif
    endif
  enddo

  write(*,*) "G1 N: ", num_g1
  write(*,*) "G2 N: ", num_g2
  write(*,*) "G1 N in region: ", count_g1

  mybox = trj%box(1)

  ! Clean up with final distance calculation
  call dist_calc(r_g1,r_g2,g1_id,mybox)
end subroutine

subroutine dist_calc(r_g1,r_g2,g1_id,box)
  use gmxfort_utils
  use parameters
  use constants
  implicit none

  integer :: i,j
  integer :: tot_sum

  real(8) :: dval
  integer, dimension(num_g1) :: g1_id
  real(8), dimension(3,3) :: box
  real(8), dimension(num_g1, 3) :: r_g1
  real(8), dimension(num_g2, 3) :: r_g2
  !real(8), dimension(num_g1,num_g2) :: dist_array
  real(8), dimension(count_g1) :: contact
 
  contact = 0.0
  write(*,*) r_g1(1,3)
  !$OMP PARALLEL DO PRIVATE(j,dval)
  do i=1,count_g1
    do j=1,count_g2
      dval = distance2(r_g1(i,:),r_g2(j,:),box)
      if (dval < r_cut_sq) then
        contact(i) = 1
        exit
      end if
    enddo
  enddo
  !$OMP END PARALLEL DO
  tot_sum =  sum(contact)
  do i=1,num_g1
    if (contact(i) == 1) then
      write(21,*) g1_id(i)
    else
      write(21,*) 0
    endif
  enddo

end subroutine



subroutine Read_Input()
  use parameters
  use constants
  open(10,file='contacts.in',status='old')
  read(10,*)
  ! File Name, 
  read(10,*) fname, iname
  read(10,*)
  ! Cutoff distance and z layer thickness above/below min value
  read(10,*) r_cut, rz_thick
  read(10,*)
  ! correlation length, separation of tos, num frames
  read(10,*) ncorr, nsep,n_calc_frames
  close(10)
  r_cut_sq = r_cut*r_cut
end subroutine

