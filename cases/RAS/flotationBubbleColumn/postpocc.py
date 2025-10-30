import openfoamparser_mai as Ofpp
import numpy as np
import matplotlib.pyplot as plt
import os


def is_float(value):
    try:
        a = float(value)
        if a > 0:
            return True
        else:
            return False
    except ValueError:
        return False


def find_time(dir_name=''):
    dir_list = os.listdir(dir_name)
    T = []
    T_name = []

    for name in dir_list:
        if is_float(name):
            T.append(float(name))
            T_name.append(name)

    combined = list(zip(T, T_name))
    sorted_combined = sorted(combined, key=lambda x: x[0])

    return zip(*sorted_combined)

path = "/home/user/OpenFOAM/user-v2412/applications/solvers/flotationReactingTwoPhaseEulerFoam/cases/RAS/flotationBubbleColumn"
T, T_name = find_time(path)
X = np.zeros(len(T))
Y = np.zeros(len(T))
for i, t in enumerate(T):
    X[i] = np.sum(Ofpp.parse_internal_field(path + "/" + T_name[i] + "/alpha.gas") * 
                  1 * 
                  Ofpp.parse_internal_field(path + "/" + T_name[i] + "/particle.gas"))
    Y[i] = np.sum(Ofpp.parse_internal_field(path + "/" + T_name[i] + "/alpha.liquid") * 
                  1000 * 
                  Ofpp.parse_internal_field(path + "/" + T_name[i] + "/particle.liquid"))

plt.plot(T, X / np.max(X), "-o")
plt.plot(T, Y / np.max(Y), "-o")
plt.show()

