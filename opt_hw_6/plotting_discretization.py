import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
links = np.array([1,2,4,8,10,100,1000,10000])
qdd = np.array([[681.375, 672.016],
    [1127.22, 1109.09],
    [1306.67, 1286.19],
    [1375.61, 1355.14],
    [1388.65, 1368.17],
    [1435.63, 1414.93],
    [1440.43, 1419.69], 
    [1440.92, 1420.17]])
solve_time = np.array([325., 508., 1025., 1967., 2198., 22681., 421416., 27775798.])*10**-6
plt.semilogx(links,qdd)

plt.xlabel('# of Discrete Segments')
plt.ylabel('Joint Acceleration ($rad/s^{2}$)')
plt.savefig('discrete_accuracy.png')
plt.show()
