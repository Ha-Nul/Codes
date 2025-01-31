import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# 초기값
a_init = 1.0

# Figure와 Axes 객체 생성
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)  # 슬라이더를 위한 공간 확보

# x축 데이터
x = np.linspace(0, 2 * np.pi, 200)

# 초기 그래프
line, = ax.plot(x, np.sin(a_init * x))
ax.set_ylim(-1.1, 1.1)

# 슬라이더 위치 및 크기
axcolor = 'lightgoldenrodyellow'
ax_a = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

# 슬라이더 생성
s_a = Slider(ax_a, 'a', 0.1, 5.0, valinit=a_init, valstep=0.1)

# 업데이트 함수
def update(val):
    a = s_a.val
    line.set_ydata(np.sin(a * x))
    ax.set_title(f'sin({a:.1f}x)')
    fig.canvas.draw_idle()

# 슬라이더 값 변경 이벤트에 업데이트 함수 연결
s_a.on_changed(update)

plt.show()