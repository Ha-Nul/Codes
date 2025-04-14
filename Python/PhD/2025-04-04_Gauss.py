import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import integrate
import matplotlib.cm as cm

class QuantumHarmonicOscillator:
    def __init__(self, m=1.0, omega=1.0, hbar=1.0):
        """
        조화진동자 파라미터 초기화
        m: 질량
        omega: 각진동수
        hbar: 플랑크 상수/2π
        """
        self.m = m
        self.omega = omega
        self.hbar = hbar
        
    def propagator(self, x_final, x_initial, t):
        """
        조화진동자의 전파자(propagator) 계산
        수치적 안정성을 위해 계산 방식 조정
        """
        # 수치적 안정성을 위해 작은 값 추가 (t가 0에 가까울 때 문제 방지)
        epsilon = 1e-10
        sin_term = np.sin(self.omega * t) + epsilon
        
        prefactor = np.sqrt(self.m * self.omega / (2 * np.pi * self.hbar * sin_term))
        exponent = (self.m * self.omega / (2 * self.hbar * sin_term)) * \
                   ((x_final**2 + x_initial**2) * np.cos(self.omega * t) - 2 * x_final * x_initial)
        
        # 지수부가 너무 크면 오버플로우 방지
        if np.abs(exponent) > 700:  # exp(709)이 약 최대 부동소수점 한계
            exponent = np.sign(exponent) * 700
            
        return prefactor * np.exp(1j * exponent)
    
    def evolve_wavefunction(self, psi_initial, x_values, t):
        """
        초기 파동함수를 시간 t만큼 발전시킴
        수치 적분 방식 개선
        """
        psi_evolved = np.zeros(len(x_values), dtype=complex)
        
        # 적분 그리드 정의 (더 조밀하게)
        x_initial_values = np.linspace(-10, 10, 1000)  # 적분 영역을 더 세밀하게 나눔
        dx = x_initial_values[1] - x_initial_values[0]
        
        # 초기 파동함수 값들
        psi_initial_values = np.array([psi_initial(x) for x in x_initial_values])
        
        for i, x_final in enumerate(x_values):
            # 각 초기 위치에서 전파자 계산
            propagators = np.array([self.propagator(x_final, x_i, t) for x_i in x_initial_values])
            
            # 수치 적분 (심프슨 법칙)
            integrand = psi_initial_values * propagators
            result = integrate.simpson(integrand, x_initial_values)
            psi_evolved[i] = result
            
        return psi_evolved
    
    def gaussian_wavepacket(self, x, x0, sigma, p0=0):
        """
        중심이 x0, 운동량이 p0, 너비가 sigma인 가우시안 파동 패킷
        """
        normalization = 1.0 / (sigma * np.sqrt(2 * np.pi))
        return normalization * np.exp(-(x - x0)**2 / (2 * sigma**2)) * np.exp(1j * p0 * x / self.hbar)
    
    def create_superposition(self, gaussian_params):
        """
        여러 가우시안의 중첩 상태 생성
        gaussian_params: (x0, sigma, weight) 튜플의 리스트
        """
        def psi(x):
            result = 0
            for x0, sigma, weight in gaussian_params:
                result += weight * self.gaussian_wavepacket(x, x0, sigma)
            return result
        return psi

# 시뮬레이션 파라미터 설정
x_min, x_max = -10, 10
nx = 100  # x 공간의 격자점 수
x_values = np.linspace(x_min, x_max, nx)
dt = 0.01  # 시간 단계 (약간 키움)
tmax = 3.0  # 최대 시간
times = np.arange(0, tmax, dt)

# 조화진동자 객체 생성
ho = QuantumHarmonicOscillator(m=1.0, omega=1.0, hbar=1.0)

# 가우시안 파라미터 설정: (중심위치, 편차, 가중치)
gaussian_params = [
    (0.0, 0.4, 1/3),  # 중심 0, 편차 0.4
    (0.0, 0.8, 1/3),  # 중심 0, 편차 0.8
    (0.0, 0.6, 1/3)   # 중심 0, 편차 0.6
]

# 초기 파동함수 (가우시안들의 중첩)
psi_initial = ho.create_superposition(gaussian_params)

# 초기 파동함수 값 계산
psi_values_initial = np.array([psi_initial(x) for x in x_values])
prob_density_initial = np.abs(psi_values_initial)**2

# 그래프 설정
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(x_min, x_max)
ax.set_ylim(0, 1.5 * np.max(prob_density_initial))
ax.set_xlabel('Position (x)')
ax.set_ylabel('Probability Density')
ax.set_title('Quantum Harmonic Oscillator: Time Evolution')
ax.grid(True, alpha=0.3)

# 각 가우시안 컴포넌트 개별 플롯 (초기상태)
colors = cm.viridis(np.linspace(0, 1, len(gaussian_params) + 1))
for i, (x0, sigma, weight) in enumerate(gaussian_params):
    gaussian = lambda x: weight * ho.gaussian_wavepacket(x, x0, sigma)
    g_values = np.array([gaussian(x) for x in x_values])
    g_prob = np.abs(g_values)**2
    ax.plot(x_values, g_prob, '--', color=colors[i], alpha=0.5, 
            label=f'Gaussian {i+1}: σ={sigma}')

# 초기 확률 밀도 플롯
line, = ax.plot(x_values, prob_density_initial, 'k-', lw=2, label='Total')
time_text = ax.text(0.02, 0.95, 'Time: 0.0', transform=ax.transAxes)
norm_text = ax.text(0.02, 0.90, 'Norm: 1.0', transform=ax.transAxes)

ax.legend()

# 시간 진화 애니메이션
def update(frame):
    t = times[frame]
    print(f"Processing frame {frame}, time = {t}")  # 진행 상황 출력
    
    # 파동함수 시간 진화
    psi_evolved_values = ho.evolve_wavefunction(psi_initial, x_values, t)
    
    # 확률 밀도 계산
    prob_density = np.abs(psi_evolved_values)**2
    
    # 정규화 전에 norm 계산 및 검사
    norm = integrate.simpson(prob_density, x_values)
    print(f"Calculated norm: {norm}")  # 디버깅용
    
    # nan 또는 매우 작은 값이면 경고 출력 후 기본값 사용
    if np.isnan(norm) or np.abs(norm) < 1e-10:
        print("Warning: Invalid norm detected, using default value")
        prob_density = prob_density_initial.copy()  # 초기 상태로 돌아감
        norm = 1.0
    else:
        prob_density /= norm

    # 정규화 후 확인
    pb_norm = integrate.simpson(prob_density, x_values)
    
    # 그래프 업데이트
    line.set_ydata(prob_density)
    time_text.set_text(f'Time: {t:.3f}')
    norm_text.set_text(f'Norm: {pb_norm:.2f}')
    
    return line, time_text

# 애니메이션 생성
ani = FuncAnimation(fig, update, frames=len(times), interval=50, blit=True)

plt.tight_layout()
plt.show()

# 특정 시간에서의 파동함수 진화를 별도로 확인하고 싶을 때
def plot_evolution_at_times(times_to_plot):
    plt.figure(figsize=(12, 8))
    
    for i, t in enumerate(times_to_plot):
        print(f"Calculating for time t = {t}")
        psi_evolved_values = ho.evolve_wavefunction(psi_initial, x_values, t)
        prob_density = np.abs(psi_evolved_values)**2
        
        # 정규화
        norm = integrate.simpson(prob_density, x_values)
        if np.isnan(norm) or np.abs(norm) < 1e-10:
            print(f"Warning: Invalid norm at t={t}, using default")
            prob_density = prob_density_initial.copy()
        else:
            prob_density /= norm
        
        plt.plot(x_values, prob_density, label=f't = {t:.3f}')
    
    plt.xlim(x_min, x_max)
    plt.xlabel('Position (x)')
    plt.ylabel('Probability Density')
    plt.title('Quantum Harmonic Oscillator: Wavefunction at Different Times')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

# 특정 시간에서의 진화 확인 (시간 간격을 더 크게)
# plot_evolution_at_times([0.0, 1.0, 2.0, 3.0, 4.0])