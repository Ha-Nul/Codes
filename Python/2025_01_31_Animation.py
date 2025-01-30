import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Figure와 Axes 객체 생성
fig, ax = plt.subplots()
x = np.linspace(0, 2 * np.pi, 200)
line, = ax.plot(x, np.sin(x))

# 업데이트 함수
def animate(i):
    line.set_ydata(np.sin(x + i / 10.0))  # y 데이터 업데이트
    return line,

# FuncAnimation 객체 생성
ani = FuncAnimation(fig, animate, frames=100, interval=20, blit=True)

# 애니메이션 표시
plt.show()

# 애니메이션 저장 (필요한 경우)
# ani.save('sine_wave.gif', writer='imagemagick')