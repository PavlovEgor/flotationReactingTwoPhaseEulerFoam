import numpy as np

N = 1000
alpha = np.linspace(-np.pi, np.pi, N)
CL = np.zeros_like(alpha)
CD = np.zeros_like(alpha)

for a in alpha:
    print(f"({(a * 180 / np.pi)} {np.pi * (np.sin(2 * (a - np.pi/2))) / (4 + np.pi * np.sin((a - np.pi/2)))} {2 * np.pi * (np.sin((a - np.pi/2)) ** 2) / (4 + np.pi * np.sin((a - np.pi/2)))})")