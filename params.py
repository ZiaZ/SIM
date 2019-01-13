
N = 20 #initial length of the system
M = 20 # number of extra cols to add when we extend, #the extra bit
k = 0.5 #spring constant 

rc = 1.05 #cut off for breaking a bond
a = 1 #lattice spacing

vacuum = 30.0 #how much space should be around the box
delta = 1.6 #remember the greek delta in the paper, 
#ratio between the stress intensity factor at the crack tip to the critical stress intensity factor. p.9.

dt = 0.25451264177891814  #time step in fm (femto seconds)
beta = 0.05  #damping factor
strain_rate = 1e-6

lm = 8 # length multiplier
T = 0 #Temperature in Kelvin
seed_lm = 0.2 #seed length as a fraction of the whole length
v_verlet = True

keep_test = False #save output traj?
desc = "test"

# this makes it easier to output parameters to a log file.
def compose_params():
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
    "Seed lm" : seed_lm,
    "descrition": desc
    }
    
    return p