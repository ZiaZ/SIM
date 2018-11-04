
N = 20 #initial length of the system
M = 20 # number of extra cols to add when we extend, #the extra bit
rc = 1.6 #cut off for breaking a bond
k = 4 #spring constant 
a = 1.3  #atom spacing

vacuum = 30.0 #how much space should be around the box
delta = 3 #remember the greek delta in the paper, 
#ratio between the stress intensity factor at the crack tip to the critical stress intensity factor. p.9.

dt = 0.025  #time step #in seconds
beta = 0.5  #damping factor
strain_rate = -1e-6

lm = 3 # length multiplier
T = 2000 #Temperature in Kelvin


p = {
    "N" : N,
    "M" : M,
    "rc" : rc,
    "k" : k,
    "a" : a,
    "vacuum" : vacuum,
    "delta" : delta,
    "dt" : dt,
    "beta" : beta,
    "strain_rate" : strain_rate,
    "lm" : lm,
    "T" : T,
}