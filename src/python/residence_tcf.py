#!/usr/bin/env python
"""
This code is part of the repository git@github.com:piskuliche/Contact-Analysis.git

It calculates the residence time TCF by using the file output previously by contact.f90

It reads the contacts.in file used by other parts of the program. 
"""


import numpy as np
import pandas


def calc_tcf(data,ncorr, nsep, ntos, nframes):
    data = data.reshape(nframes,int(len(data)/nframes))
    occ = (data>0)*1
    np.savetxt("sums.vals",np.c_[np.sum(occ,axis=0)])
    tcf = []
    for outer in range(ntos):
        t_out = outer*nsep
        tcf_temp = np.zeros(ncorr)
        counts = np.zeros(ncorr)
        #prev = occ[t_out]
        prev = data[t_out]
        prev = prev[prev!=0]
        for inner in range(ncorr):
            t_in = t_out + inner
            # Select and remove zeros
            curr = data[t_in]
            curr = curr[curr!=0]
            # Get Intersect
            the_intersect = np.intersect1d(prev,curr) # Don't tell the ring
            counts[inner] = len(the_intersect) # or approximately five seasons, dammit bartowski!
            # Calculate contribution to TCF
            tcf_temp[inner] = counts[inner]/counts[0]
            prev = the_intersect.copy()
        tcf.append(tcf_temp)
    final = np.average(tcf,axis=0)
    np.savetxt("final_tcf.dat",np.c_[final])
    return

def read_input():
    ncorr, nsep, n_cal_frames
    with open("contacts.in",'r') as f:
        f.readline()
        fname, iname = f.readline().strip().split()
        f.readline()
        r_cut, rz_thick = f.readline().strip().split()
        f.readline()
        ncorr, nsep, n_calc_frames = f.readline().strip().split()

        ncorr = int(ncorr)
        nsep = int(nsep)
        n_calc_frames = int(n_calc_frames)
    return ncorr, nsep, n_calc_frames


if __name__ == "__main__":
    ncorr, nsep, nframes = read_input()
    ntos = nframes/nsep - ncorr/nsep

    dvals = pandas.read_csv("contacts.out",header=None)
    dvals=dvals.to_numpy()

    calc_tcf(dvals,ncorr,nsep,ntos,nframes)
