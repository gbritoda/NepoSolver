import lib.nepo_tools as nepo
import numpy as np
import cmath

np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)

#Test Ybus
Ybus = np.array([
    [-5.5, 1/0.4, 1/0.5, 0, 0],
    [1/0.4, -11.5, 1/0.25, 0, 1/0.2],
    [1/0.5, 1/0.25, -14, 1/0.125, 0],
    [0, 0, 1/0.125, -10, 1/0.5],
    [0, 1/0.2, 0, 1/0.5, -7.8]
])
Ybus = 1j*Ybus

np.set_printoptions(precision=4)

zbus = np.array([
    [0.2, 0.2, 0.2, 0.2],
    [0.2, 0.6, 0.2, 0.6],
    [0.2, 0.2, 0.8, 0.2],
    [0.2, 0.6, 0.2, 1.1]
])
zbus = 1j*zbus

zbus2 = 1j*np.array([
    [0.717, 0.61, 0.533, 0.58],
    [0.61, 0.732, 0.64, 0.697],
    [0.533, 0.64, 0.717, 0.67],
    [0.58, 0.697, 0.67, 0.763]
])

#ex 14
zbus3 = 1j*np.array([
    [0.15, 0.08, 0.04, 0.07],
    [0.08, 0.15, 0.06, 0.09],
    [0.04, 0.06, 0.13, 0.05],
    [0.07, 0.09, 0.05, 0.12]
])
# zbus3[-1][-1] = 0.12

#ex 15
ybus15 = 1j*np.array([
    [-12, 5, 2],
    [5, -7.5, 2.5],
    [2, 2.5, -8.5]
])

zbus15 = nepo.inv(ybus15)

#ex 17
zbus17 = 1j*np.array([
    [0.244, 0.194, 0.156, 0.146],
    [0.194, 0.23, 0.15, 0.15],
    [0.155, 0.15, 0.2, 0.105],
    [0.146, 0.15, 0.105, 0.196]
])