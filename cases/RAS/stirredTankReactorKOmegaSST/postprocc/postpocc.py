import openfoamparser_mai as Ofpp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


delay = 10.6

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

path = "/home/user/OpenFOAM/user-v2412/applications/solvers/flotationReactingTwoPhaseEulerFoam/cases/RAS/stirredTankReactorKOmegaSST"
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

# plt.plot(T, X / np.max(X), "-o")
T = np.array(T)
Y = Y[T > delay]
Y = Y / np.max(Y)
# Y = 1 - (1 - Y) * 7
T = T[T > delay] - delay

plt.plot(T, Y, "--", label="30 micron curent work")
# plt.xlim(delay)

df_240 = pd.read_csv(
    "240mic.csv",
    sep=' ',          
    decimal=',',      
    header=None,     
    engine='python'   
)

df_60 = pd.read_csv(
    "60mic.csv",
    sep=' ',          
    decimal=',',      
    header=None,     
    engine='python'   
)

df_30 = pd.read_csv(
    "30mic.csv",
    sep=' ',          
    decimal=',',      
    header=None,     
    engine='python'   
)

df_30_2 = pd.read_csv(
    "30mic2.csv",
    sep=' ',          
    decimal=',',      
    header=None,     
    engine='python'   
)

df_15 = pd.read_csv(
    "15mic.csv",
    sep=' ',          
    decimal=',',      
    header=None,     
    engine='python'   
)

plt.plot(df_240[0].array * 7, df_240[1].array, label="240 micron KohSchwarz")
plt.plot(df_60[0].array * 7, df_60[1].array, label="60 KohSchwarz")
# plt.plot(df_30[0].array, df_30[1].array, label="30 KohSchwarz")
plt.plot(df_30_2[0].array * 7, df_30_2[1].array, label="30 KohSchwarz")
plt.plot(df_15[0].array * 7, df_15[1].array, label="15 KohSchwarz")

# plt.xscale('log')  # Логарифмический масштаб по оси X
plt.yscale('log')  # Логарифмический масштаб по оси Y

# Настройки графика
plt.xlabel('Время работы установки')
plt.ylabel('относительная массовая доля частиц')
plt.grid(True, which="both", ls="--")  # Сетка для логарифмического масштаба
plt.legend()

plt.xlim(0, 200)
plt.ylim(0.4, 1)

plt.show()

