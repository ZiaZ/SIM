#!/usr/bin/python
# -*- coding: utf-8 -*-

import run_ideal_brittle_solid as ribs
import params
import numpy as np
import time

start_time = time.time()
#-----------------Delta Section Sim-------------------

# deltas = np.arange(2.1, 2.3, 0.1)
# for D in deltas:
#     params.delta = D
#     ribs.ribs(params)

#-------------------Delta Section Sim -----------------

#-----------------Temperature Section Sim---------------
# deltas = [1.4,1.5,1.6,1.8,2.0]
# params.keep_test = True
# params.k = 5
# for delta in deltas:
#     temps = [-50,0,10,25,35,45,65,85,100,130]
#     temps_kelvin = [ el+273 for el in temps]
#     params.delta = delta
#     for temp in temps_kelvin:
#         if delta == 1.55 and temp < (65+273):
#             continue
#         params.desc = '(T='+str(temp-273)+')_(D='+str(delta)+')'
#         params.T = temp
#         ribs.ribs(params, frame_count = 1200)
#-----------------Temperature Section Sim---------------

#-----------------Void Section Sim-------------------



#-------------------Void Section Sim -----------------
    
end_time = time.time()
time_diff = (end_time - start_time)
time_unit = "seconds"
if time_diff > 60:
    time_diff /= 60
    time_unit = 'minutes'
print("the simulation took |%s %s| in total." %(time_diff, time_unit))       


import subprocess
subprocess.call(["shutdown", "/s"])