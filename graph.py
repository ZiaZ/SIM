import matplotlib.pyplot as plt

values = open('tip_x.txt').readlines()
values2 = [val.split(',') for val in values]
x_pos = [float(val[1].strip('\n')) for val in values2]
y_pos = [ float(val[0]) for val in values2]


plt.plot(x_pos, y_pos)
plt.show()