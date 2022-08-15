# Contact-Analysis
This is a code for calculating contacts within a certain range of distances from simulations in GROMACS.

## Installation

Note - this code requires the package

gfortran contact.f90 -I $lb\_gmx\_inc -L $lb\_gmx\_lib -lgmxfort -o contact -Ofast -fopenmp


