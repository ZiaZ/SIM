import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import csv
import ase.units as units

def get_xy(sim_index = 125):
    dir = ".simout/sim_"+str(sim_index)+"_true_deltas"
    values = open(dir+'/tip_x.txt').readlines()
    values2 = [val.split(',') for val in values]
    #x_pos = [float(val[1].strip('\n')) for val in values2]
    y_pos = [ float(val[0]) for val in values2]
    x_pos = [0]
    for i in range(len(y_pos) - 1):
        x_pos.append(round((x_pos[-1]+0.025), 3))
    if (sim_index == 133):
        plt.plot(x_pos,y_pos, 'or')
        plt.show()
    return [x_pos, y_pos]




x_pos_D = []
y_pos_D = []
errors =  []
sim_range = []



# y_pos_D /= 0.398101

with open(".simout/DandSimNo.csv",'r') as file:
    rows = csv.reader(file)
    for row in rows:
        sim_range.append([float(row[0]),int(row[1])])

for i in sim_range:
    x_pos_D.append(i[0])
    xy = get_xy(i[1])
    m, c, r_value, p_value, std_err = stats.linregress(xy[0], xy[1])
    # m, c, err = np.polyfit(xy[0],xy[1], 1, full = True)
    y_pos_D.append(m)
    errors.append(std_err)

errors = np.array(errors)
x_pos_D = np.array(x_pos_D)
y_pos_D = np.array(y_pos_D)

# del(x_pos_D[:4]); del(x_pos_D[-13:])
# del(y_pos_D[:4]); del(y_pos_D[-13:])

plt.xlabel('Strain (Δ)')
plt.ylabel('Crack Velocity')
plt.title('Relationship of crack speed against strain.')
plt.plot(x_pos_D, y_pos_D,'bo',x_pos_D, y_pos_D,'r')

plt.plot(x_pos_D, 5*(1 - 1/x_pos_D**2),'g-')


plt.xticks(np.arange(min(x_pos_D) - 0.1, max(x_pos_D) + 0.1, 0.05))
# plt.yticks(np.arange(min(y_pos_D), max(y_pos_D), 0.2))
plt.yticks(np.arange(-1, max(y_pos_D) + 0.5, 0.2))
plt.errorbar(x_pos_D,y_pos_D,yerr=errors*100, fmt="g*")
plt.xlim(1.3,2.1)
plt.ylim(-1,max(y_pos_D) + 0.5)
plt.grid()
plt.show()
