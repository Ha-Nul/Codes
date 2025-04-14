import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import scipy.stats as stats
import seaborn as sns

import Equilb as eq

coord = np.linspace(0,1,100)
momen = np.linspace(0,1,100)

Energy = 1
mass = 1
spconst = 1

sample_count = 2000

test = []
for i in range(sample_count):
  result = eq.Hamconst(coord, momen, mass, spconst, Energy)
  #print(result)
  if result is not None:  # Only append if result is not None
    test.append(result)

tup_test = tuple(map(tuple,test))
unique_tup = set(tup_test)

# 새로운 그림 생성
fig = plt.figure(figsize=(16, 7))

# 두 개의 서브플롯 생성 (1x2 그리드)
ax1 = fig.add_subplot(1, 2, 1, projection='3d')  # 왼쪽 3D 플롯
ax2 = fig.add_subplot(1, 2, 2)  # 오른쪽 colormap 플롯

# 3D 플롯 데이터 (ax1)
x = [t[0] for t in unique_tup if t is not None]
y = [t[1] for t in unique_tup if t is not None]
z = [t[2] for t in unique_tup if t is not None]

ax1.scatter(x, y, z, c=z, cmap='viridis')  # Use z values for color mapping

x_2 = [-t[0] for t in unique_tup if t is not None]
y_2 = [t[1] for t in unique_tup if t is not None]
z_2 = [t[2] for t in unique_tup if t is not None]

ax1.scatter(x_2, y_2, z_2, c=z_2, cmap='viridis')

x_3 = [-t[0] for t in unique_tup if t is not None]
y_3 = [-t[1] for t in unique_tup if t is not None]
z_3 = [t[2] for t in unique_tup if t is not None]

ax1.scatter(x_3, y_3, z_3, c=z_3, cmap='viridis')

x_4 = [t[0] for t in unique_tup if t is not None]
y_4 = [-t[1] for t in unique_tup if t is not None]
z_4 = [t[2] for t in unique_tup if t is not None]

ax1.scatter(x_4, y_4, z_4, c=z_4, cmap='viridis')

ax1.set_xlabel('qnorm')
ax1.set_ylabel('mnorm')
ax1.set_zlabel('Hamiltonian')
ax1.view_init(elev=90, azim=-180)
ax1.set_title("3D Hamiltonian Plot (Sample count): "+str(sample_count))

# Colormap 플롯 데이터 (ax2)
x_coords = [t[0] for t in unique_tup]
y_coords = [t[1] for t in unique_tup]

kernel = stats.gaussian_kde(np.vstack([x_coords, y_coords]))

x_grid = np.linspace(min(x_coords), max(x_coords), len(x_coords))
y_grid = np.linspace(min(y_coords), max(y_coords), len(y_coords))
X, Y = np.meshgrid(x_grid, y_grid)
Z = np.reshape(kernel(np.vstack([X.ravel(), Y.ravel()])).T, X.shape)

contour = ax2.contourf(X, Y, Z, cmap='viridis')
fig.colorbar(contour, ax=ax2)
ax2.set_xlabel("qnorm")
ax2.set_ylabel("mnorm")
ax2.set_title("Probability Density Function (Sample count): " + str(sample_count))

plt.tight_layout()
plt.show()

