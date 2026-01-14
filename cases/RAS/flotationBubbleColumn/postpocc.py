import openfoamparser_mai as Ofpp
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


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

def find_slope(data, dx) -> float:

    # Вычисление производной
    derivative = np.diff(data)

    # Нахождение индекса максимальной производной
    max_derivative_index = np.argmax(np.abs(derivative)[:-10])

    # Координата скачка
    jump_coordinate = dx * max_derivative_index

    return jump_coordinate

path = "/home/user/OpenFOAM/user-v2412/applications/solvers/flotationReactingTwoPhaseEulerFoam/cases/RAS/flotationBubbleColumn"
T, T_name = find_time(path)
shape = (150, 50)
alpha_0 = 5.2
# X = np.zeros(len(T))
# Y = np.zeros(len(T))
# for i, t in enumerate(T):
#     X[i] = np.sum(Ofpp.parse_internal_field(path + "/" + T_name[i] + "/alpha.gas") * 
#                   1 * 
#                   Ofpp.parse_internal_field(path + "/" + T_name[i] + "/particle.gas"))
#     Y[i] = np.sum(Ofpp.parse_internal_field(path + "/" + T_name[i] + "/alpha.liquid") * 
#                   1000 * 
#                   Ofpp.parse_internal_field(path + "/" + T_name[i] + "/particle.liquid"))

# plt.plot(T, X / np.max(X), "-o")
# plt.plot(T, Y / np.max(Y), "-o")
# plt.show()

x = np.zeros(len(T))
a = np.zeros(shape[1])

for i, t in enumerate(T):

    alpha = Ofpp.parse_internal_field(path + "/" + T_name[i] + "/alpha.gas").reshape(shape)

    for j in range(shape[1]):
        a[j] = find_slope(alpha[:, j], 1)
    x[i] = np.mean(a)

plt.plot(T, x)
plt.show()

# df = pd.DataFrame({f"T_{alpha_0}": T, f"x_{alpha_0}": x/shape[0]})
df = pd.read_csv(path + f"/df_10.csv")

df[f"T_{alpha_0}"] = T
df[f"x_{alpha_0}"] = x/shape[0]


df.to_csv(path + f"/df_10.csv", index=False)
n = len([col for col in df.columns if col.startswith('X_')])
for column in df.columns:
    if column[0] == "T":
        plt.plot(df[column].array, df["x" + column[1:]].array, label=column[2:])

plt.legend()
plt.show()