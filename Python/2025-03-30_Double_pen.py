import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

# 물리적 상수
g = 9.8  # 중력 가속도 (m/s^2)
L1, L2 = 1.0, 1.0  # 진자 길이 (m)
m1, m2 = 1.0, 1.0  # 질량 (kg)

# 이중 진자의 운동 방정식
def double_pendulum(t, y):
    theta1, omega1, theta2, omega2 = y
    
    # 운동 방정식에 사용할 변수들
    delta = theta2 - theta1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta) * np.cos(delta)
    den2 = (L2 / L1) * den1
    
    # 미분방정식
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    
    domega1_dt = ((m2 * L1 * omega1 * omega1 * np.sin(delta) * np.cos(delta)
                  + m2 * g * np.sin(theta2) * np.cos(delta)
                  + m2 * L2 * omega2 * omega2 * np.sin(delta)
                  - (m1 + m2) * g * np.sin(theta1))
                  / den1)
    
    domega2_dt = ((-m2 * L2 * omega2 * omega2 * np.sin(delta) * np.cos(delta)
                  + (m1 + m2) * g * np.sin(theta1) * np.cos(delta)
                  - (m1 + m2) * L1 * omega1 * omega1 * np.sin(delta)
                  - (m1 + m2) * g * np.sin(theta2))
                  / den2)
    
    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt]

# 초기 조건 설정
# [theta1, omega1, theta2, omega2]
# 약간의 차이를 둔 두 가지 초기 조건
y0_1 = [np.pi/2, 0, np.pi/2, 0]  # 첫 번째 초기 조건
y0_2 = [np.pi/2 + 0.01, 0, np.pi/2, 0]  # 두 번째 초기 조건 (약간 다름)

# 시간 범위
t_span = (0, 30)
t_eval = np.linspace(t_span[0], t_span[1], 2000)

# 미분방정식 풀기
solution1 = solve_ivp(double_pendulum, t_span, y0_1, method='RK45', t_eval=t_eval)
solution2 = solve_ivp(double_pendulum, t_span, y0_2, method='RK45', t_eval=t_eval)

# 결과 추출
t = solution1.t
theta1_1, omega1_1, theta2_1, omega2_1 = solution1.y
theta1_2, omega1_2, theta2_2, omega2_2 = solution2.y

# 위상 공간 시각화
plt.figure(figsize=(15, 10))

# 첫 번째 진자의 위상 공간 (θ1, ω1)
plt.subplot(2, 2, 1)
plt.plot(theta1_1, omega1_1, 'b-', alpha=0.7, label='Initial Condition 1')
plt.plot(theta1_2, omega1_2, 'r-', alpha=0.7, label='Initial Condition 2')
plt.grid(True)
plt.xlabel(r'$\theta_1$ [rad]')
plt.ylabel(r'$\omega_1$ [rad/s]')
plt.title('Phase Space of the First Pendulum')
plt.legend()

# 두 번째 진자의 위상 공간 (θ2, ω2)
plt.subplot(2, 2, 2)
plt.plot(theta2_1, omega2_1, 'b-', alpha=0.7, label='Initial Condition 1')
plt.plot(theta2_2, omega2_2, 'r-', alpha=0.7, label='Initial Condition 2')
plt.grid(True)
plt.xlabel(r'$\theta_2$ [rad]')
plt.ylabel(r'$\omega_2$ [rad/s]')
plt.title('Phase Space of the Second Pendulum')
plt.legend()

# 첫 번째 진자 vs 두 번째 진자 각도 (θ1, θ2)
plt.subplot(2, 2, 3)
plt.plot(theta1_1, theta2_1, 'b-', alpha=0.7, label='Initial Condition 1')
plt.plot(theta1_2, theta2_2, 'r-', alpha=0.7, label='Initial Condition 2')
plt.grid(True)
plt.xlabel(r'$\theta_1$ [rad]')
plt.ylabel(r'$\theta_2$ [rad]')
plt.title('θ1 vs θ2')
plt.legend()

# 시간에 따른 진자 각도 변화
plt.subplot(2, 2, 4)
plt.plot(t, theta1_1, 'b-', alpha=0.7, label='θ1 (IC1)')
plt.plot(t, theta2_1, 'g-', alpha=0.7, label='θ2 (IC1)')
plt.plot(t, theta1_2, 'r--', alpha=0.7, label='θ1 (IC2)')
plt.plot(t, theta2_2, 'm--', alpha=0.7, label='θ2 (IC2)')
plt.grid(True)
plt.xlabel('Time [s]')
plt.ylabel('Angle [rad]')
plt.title('Angles vs Time')
plt.legend()

plt.tight_layout()
plt.show()

# 애니메이션을 위한 이중 진자 좌표 계산 함수
def get_coords(theta1, theta2):
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)
    
    return x1, y1, x2, y2

# 애니메이션 생성 (첫 번째 초기 조건만 사용)
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-2.2, 2.2)
ax.grid(True)

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], 'r-', alpha=0.3)  # 두 번째 진자의 궤적

x1_data, y1_data, x2_data, y2_data = get_coords(theta1_1, theta2_1)
trace_x, trace_y = [], []

def init():
    line.set_data([], [])
    trace.set_data([], [])
    return line, trace

def animate(i):
    x1, y1, x2, y2 = x1_data[i], y1_data[i], x2_data[i], y2_data[i]
    
    trace_x.append(x2)
    trace_y.append(y2)
    
    # 마지막 50개 위치만 그림
    line.set_data([0, x1, x2], [0, y1, y2])
    trace.set_data(trace_x[-100:], trace_y[-100:])
    
    return line, trace

ani = animation.FuncAnimation(fig, animate, frames=len(t), 
                             interval=20, blit=True, init_func=init)

plt.title('Double Pendulum Animation')
plt.tight_layout()
plt.show()

# 특정 시점을 저장하고 싶다면 아래와 같이 할 수 있습니다:
# ani.save('double_pendulum.gif', writer='pillow', fps=30)