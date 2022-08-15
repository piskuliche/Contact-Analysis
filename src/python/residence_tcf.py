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
            # absorbing boundary conditions
            #curr = occ[t_in]
            curr = data[t_in]
            curr = curr[curr!=0]
            #counts[inner] = np.sum(prev*curr)
            the_intersect = np.intersect1d(prev,curr) # Don't tell the ring
            counts[inner] = len(the_intersect) # or approximately five seasons, dammit bartowski!
            tcf_temp[inner] = counts[inner]/counts[0]
            prev = the_intersect.copy()
        tcf.append(tcf_temp)
    final = np.average(tcf,axis=0)
    np.savetxt("final_tcf.dat",np.c_[final])
    return

#dvals = np.loadtxt("contacts.out",usecols=0,dtype=int)
dvals = pandas.read_csv("contacts.out",header=None)
dvals=dvals.to_numpy()
calc_tcf(dvals,4000,25,40,5000)
