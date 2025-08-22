import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

xdata = np.linspace(0, 10, 100)

fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(-1.5, 1.5)

line, = ax.plot([],[], lw=2)

def update(frame):
    ydata = np.sin(xdata + frame * 0.1)
    line.set_data(xdata,ydata)

    return line,

ani = FuncAnimation(fig, update, frames=100,interval=50, blit=True)

plt.show()