import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import integrate
import matplotlib.cm as cm

class QuantumBrownianMotion:
    def __init__(self, m=1.0, omega=1.0, hbar=1.0, gamma=0.1, kT=0.5):
        """
        양자 브라운 운동 모델 파라미터 초기화
        m: 질량
        omega: 각진동수
        hbar: 플랑크 상수/2π
        gamma: 감쇠 상수 (환경과의 결합 강도)
        kT: 환경의 온도 (볼츠만 상수 × 온도)
        """
        self.m = m
        self.omega = omega
        self.hbar = hbar
        self.gamma = gamma
        self.kT = kT
        
    def propagator(self, x_final, x_initial, t):
        """
        양자 브라운 운동 모델의 전파자(propagator) 계산
        Caldeira-Leggett 모델 기반 근사 전파자
        """
        # 수치적 안정성을 위한 작은 값
        epsilon = 1e-10
        
        # 감쇠된 진동수 계산
        omega_d = np.sqrt(max(self.omega**2 - (self.gamma/2)**2, 0))
        
        if omega_d > epsilon:  # 약한 감쇠 경우 (underdamped)
            # 지수 감쇠 인자
            exp_factor = np.exp(-self.gamma * t / 2)
            
            # 사인 항에 0 나눔 방지
            sin_term = np.sin(omega_d * t) + epsilon
            
            # 유효 질량 (시간 의존)
            m_eff = self.m * exp_factor / sin_term
            
            # 위상 인자
            phase = (self.m / (2 * self.hbar)) * (
                (x_final**2 + x_initial**2) * (self.omega / np.tan(omega_d * t + epsilon)) - 
                2 * x_final * x_initial * (self.omega / sin_term)
            ) * exp_factor
            
            # 확산 인자 (온도에 의한 효과)
            diffusion = self.kT * (1 - np.exp(-self.gamma * t)) / self.hbar
            diff_term = (-(x_final - x_initial * np.exp(-self.gamma * t / 2))**2 * diffusion)
            
            # 전파자 계산 (복소수 형태)
            prefactor = np.sqrt(m_eff * omega_d / (2 * np.pi * self.hbar))
            
            # 수치적 안정성을 위한 지수 제한
            if np.abs(phase) > 700:
                phase = np.sign(phase) * 700
            if np.abs(diff_term) > 700:
                diff_term = -700  # 항상 음수여야 함
            
            return prefactor * np.exp(1j * phase + diff_term)
        
        else:  # 강한 감쇠 경우 (overdamped)
            # 이 경우 다른 수식 사용 필요
            # 간단한 근사로 대체
            lambda1 = self.gamma/2 + np.sqrt((self.gamma/2)**2 - self.omega**2)
            lambda2 = self.gamma/2 - np.sqrt((self.gamma/2)**2 - self.omega**2)
            
            prefactor = np.sqrt(self.m / (2 * np.pi * self.hbar * t))
            
            # 고전적 작용 계산
            action = (self.m / (2 * self.hbar)) * (
                (x_final**2 + x_initial**2) * (self.gamma / (1 - np.exp(-self.gamma * t))) - 
                2 * x_final * x_initial / (np.sinh(self.gamma * t / 2) + epsilon)
            )
            
            # 확산 인자
            diffusion = self.kT * (1 - np.exp(-self.gamma * t)) / self.hbar
            diff_term = (-(x_final - x_initial * np.exp(-self.gamma * t / 2))**2 * diffusion)
            
            # 수치적 안정성을 위한 지수 제한
            if np.abs(action) > 700:
                action = np.sign(action) * 700
            if np.abs(diff_term) > 700:
                diff_term = -700
            
            return prefactor * np.exp(1j * action + diff_term)
    
    def evolve_density_matrix(self, rho_initial, x_values, t):
        """
        초기 밀도 행렬을 시간 t만큼 발전시킴
        양자 브라운 운동에서는 밀도 행렬 접근이 더 적합
        """
        nx = len(x_values)
        rho_evolved = np.zeros((nx, nx), dtype=complex)
        
        # 적분 그리드 정의
        x_grid = np.linspace(-10, 10, 100)
        dx = x_grid[1] - x_grid[0]
        
        for i, x_final in enumerate(x_values):
            for j, y_final in enumerate(x_values):
                # 이중 적분 수행 (x_initial, y_initial에 대해)
                integral = 0.0
                
                for x_initial in x_grid:
                    for y_initial in x_grid:
                        # 위치 인덱스 계산
                        idx_x = np.argmin(np.abs(x_values - x_initial))
                        idx_y = np.argmin(np.abs(x_values - y_initial))
                        
                        # 초기 밀도 행렬 값
                        rho_init_val = rho_initial[idx_x, idx_y]
                        
                        # 전파자 곱 계산
                        propagator_term = (
                            self.propagator(x_final, x_initial, t) * 
                            np.conj(self.propagator(y_final, y_initial, t))
                        )
                        
                        # 적분에 기여
                        integral += rho_init_val * propagator_term * dx**2
                
                rho_evolved[i, j] = integral
        
        return rho_evolved
    
    def reduce_grid_density_matrix(self, rho, x_values):
        """
        계산 효율성을 위해 밀도 행렬의 그리드 크기 줄이기
        """
        nx = len(x_values)
        if nx <= 20:  # 이미 작은 그리드면 그대로 반환
            return rho, x_values
        
        # 그리드 크기 줄이기
        step = nx // 20
        reduced_x = x_values[::step]
        reduced_rho = rho[::step, ::step]
        
        return reduced_rho, reduced_x
    
    def evolve_wavefunction(self, psi_initial, x_values, t):
        """
        파동함수 시간 발전 (브라운 운동에 맞게 조정)
        simpson 적분 사용하여 계산 속도 향상
        """
        psi_evolved = np.zeros(len(x_values), dtype=complex)
        
        # 적분 범위 설정
        x_min, x_max = -10, 10
        
        for i, x_final in enumerate(x_values):
            # 적분할 함수 정의 (클로저 사용)
            def integrand(x_initial):
                # 초기 파동함수 값 계산
                psi_init_val = psi_initial(x_initial)
                # 전파자 계산
                prop = self.propagator(x_final, x_initial, t)
                return psi_init_val * prop
            
            # scipy.integrate.simpson로 적분 수행
            # simpson는 복소수 함수를 직접 적분할 수 없으므로 실수부와 허수부를 분리
            real_part, _ = integrate.quad(lambda x: np.real(integrand(x)), x_min, x_max)
            imag_part, _ = integrate.quad(lambda x: np.imag(integrand(x)), x_min, x_max)
            
            psi_evolved[i] = real_part + 1j * imag_part
        
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
        gaussian_params: (x0, sigma, p0, weight) 튜플의 리스트
        """
        def psi(x):
            result = 0
            for x0, sigma, p0, weight in gaussian_params:
                result += weight * self.gaussian_wavepacket(x, x0, sigma, p0)
            return result
        return psi
    
    def calculate_probability_density(self, psi_values):
        """
        파동함수로부터 확률 밀도 계산
        """
        return np.abs(psi_values)**2
    
    def calculate_density_matrix_diag(self, rho):
        """
        밀도 행렬의 대각선 요소 (확률 밀도) 추출
        """
        return np.real(np.diag(rho))

# 시뮬레이션 파라미터
x_min, x_max = -3, 3
nx = 3  # x 공간의 격자점 수
x_values = np.linspace(x_min, x_max, nx)
dt = 0.001  # 시간 단계
tmax = 3.0  # 최대 시간
times = np.arange(0, tmax, dt)

# 양자 브라운 운동 객체 생성 (감쇠 및 온도 추가)
qbm = QuantumBrownianMotion(m=1.0, omega=1.0, hbar=1.0, gamma=0.2, kT=0.5)

# 가우시안 파라미터 설정: (중심위치, 편차, 운동량, 가중치)
gaussian_params = [
    (0.0, 0.5, 1.0, 0.6),  # 중심 0, 편차 0.5, 운동량 1.0
    (-2.0, 0.4, -0.5, 0.4)  # 중심 -2, 편차 0.4, 운동량 -0.5
]

# 초기 파동함수 (가우시안들의 중첩)
psi_initial = qbm.create_superposition(gaussian_params)

# 초기 파동함수 값 계산
psi_values_initial = np.array([psi_initial(x) for x in x_values])
prob_density_initial = qbm.calculate_probability_density(psi_values_initial)

# 초기 밀도 행렬 생성 (순수 상태)
rho_initial = np.zeros((nx, nx), dtype=complex)
for i in range(nx):
    for j in range(nx):
        rho_initial[i, j] = psi_values_initial[i] * np.conj(psi_values_initial[j])

# 그래프 설정
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# 상단 그래프: 확률 밀도
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(0, 1.5 * np.max(prob_density_initial))
ax1.set_xlabel('Position (x)')
ax1.set_ylabel('Probability Density')
ax1.set_title('Quantum Brownian Motion: Probability Distribution')
ax1.grid(True, alpha=0.3)

# 하단 그래프: 위상 공간 표현 (불확정성 표현)
extent = [x_min, x_max, x_min, x_max]
im = ax2.imshow(np.abs(rho_initial), extent=extent, origin='lower', cmap='viridis', 
                vmin=0, vmax=np.max(np.abs(rho_initial)))
ax2.set_xlabel('Position (x)')
ax2.set_ylabel('Position (y)')
ax2.set_title('Density Matrix Components')
plt.colorbar(im, ax=ax2, label='|ρ(x,y)|')

# 눈금 추가
x_ticks = np.linspace(x_min, x_max, 5)
ax2.set_xticks(x_ticks)
ax2.set_yticks(x_ticks)
ax2.grid(False)  # 그리드는 행렬 시각화에 방해가 될 수 있어 비활성화

# 각 가우시안 컴포넌트 개별 플롯 (초기상태)
colors = cm.tab10(np.linspace(0, 1, len(gaussian_params)))
for i, (x0, sigma, p0, weight) in enumerate(gaussian_params):
    gaussian = lambda x: weight * qbm.gaussian_wavepacket(x, x0, sigma, p0)
    g_values = np.array([gaussian(x) for x in x_values])
    g_prob = qbm.calculate_probability_density(g_values)
    ax1.plot(x_values, g_prob, '--', color=colors[i], alpha=0.5, 
            label=f'Gaussian {i+1}: x₀={x0}, p₀={p0}')

# 초기 확률 밀도 플롯
line, = ax1.plot(x_values, prob_density_initial, 'k-', lw=2, label='Total')
time_text = ax1.text(0.02, 0.95, 'Time: 0.0', transform=ax1.transAxes)
norm_text = ax1.text(0.02, 0.90, 'Norm: 1.0', transform=ax1.transAxes)
purity_text = ax1.text(0.02, 0.85, 'Purity: 1.0', transform=ax1.transAxes)

ax1.legend()

# 시간 진화 애니메이션
def update(frame):
    t = times[frame]
    print(f"Processing frame {frame}, time = {t}")
    
    # 파동함수 시간 진화
    try:
        psi_evolved_values = qbm.evolve_wavefunction(psi_initial, x_values, t)
        prob_density = qbm.calculate_probability_density(psi_evolved_values)
        
        # 밀도 행렬 진화 (효율성을 위해 그리드 크기 줄임)
        reduced_rho, reduced_x = qbm.reduce_grid_density_matrix(rho_initial, x_values)
        rho_evolved = qbm.evolve_density_matrix(reduced_rho, reduced_x, t)
        
        # 정규화 전에 norm 계산
        norm = integrate.simpson(prob_density, x_values)
        print(f"Calculated norm: {norm}")  # 디버깅용
        
        # 순도(purity) 계산 - Tr(ρ²)
        purity = np.real(np.trace(np.matmul(rho_evolved, rho_evolved)))
        print(f"Calculated purity: {purity}")  # 디버깅용
        
        # nan이면 경고 출력 후 기본값 사용
        if np.isnan(norm) or np.abs(norm) < 1e-10:
            print("Warning: Invalid norm detected")
            prob_density = prob_density_initial.copy()
            norm = 1.0
            purity = 1.0
            rho_display = np.abs(rho_initial)
        else:
            # norm이 변할 것으로 예상됨 (브라운 운동 특성상)
            # 그러나 그래프 비교를 위해 정규화
            prob_density = prob_density / norm
            
            # 밀도 행렬 표시용 데이터 준비
            rho_display = np.abs(rho_evolved)
            
            # 밀도 행렬 범위 설정 (필요시 수정)
            # reduced_x를 사용하는 경우 extent도 업데이트 필요
            if len(reduced_x) != len(x_values):
                im.set_extent([reduced_x.min(), reduced_x.max(), reduced_x.min(), reduced_x.max()])
        
        # 그래프 업데이트
        line.set_ydata(prob_density)
        im.set_data(rho_display)
        im.set_clim(0, np.max(rho_display))
        
        time_text.set_text(f'Time: {t:.3f}')
        norm_text.set_text(f'Norm: {norm:.4f}')
        purity_text.set_text(f'Purity: {purity:.4f}')
        
        return line, im, time_text, norm_text, purity_text
    
    except Exception as e:
        print(f"Error in frame {frame}: {e}")
        return line, im, time_text, norm_text, purity_text

# 애니메이션 생성
ani = FuncAnimation(fig, update, frames=len(times), interval=100, blit=True)

plt.tight_layout()
plt.show()

# 특정 시간에서의 파동함수 진화를 별도로 확인하고 싶을 때
def plot_evolution_at_times(times_to_plot):
    fig, axes = plt.subplots(len(times_to_plot), 2, figsize=(14, 4*len(times_to_plot)))
    
    for i, t in enumerate(times_to_plot):
        print(f"Calculating for time t = {t}")
        
        try:
            # 파동함수 진화
            psi_evolved_values = qbm.evolve_wavefunction(psi_initial, x_values, t)
            prob_density = qbm.calculate_probability_density(psi_evolved_values)
            
            # 밀도 행렬 진화
            reduced_rho, reduced_x = qbm.reduce_grid_density_matrix(rho_initial, x_values)
            rho_evolved = qbm.evolve_density_matrix(reduced_rho, reduced_x, t)
            
            # Norm 계산
            norm = integrate.simpson(prob_density, x_values)
            
            # 순도(purity) 계산
            purity = np.real(np.trace(np.matmul(rho_evolved, rho_evolved)))
            
            # 확률 밀도 플롯
            axes[i, 0].plot(x_values, prob_density / np.max(prob_density), 'b-', lw=2)
            axes[i, 0].set_xlim(x_min, x_max)
            axes[i, 0].set_ylabel('Prob. Density')
            axes[i, 0].set_title(f't = {t:.1f}, Norm = {norm:.4f}, Purity = {purity:.4f}')
            axes[i, 0].grid(True, alpha=0.3)
            
            # 밀도 행렬 플롯
            im = axes[i, 1].imshow(np.abs(rho_evolved), origin='lower', cmap='viridis',
                              extent=[reduced_x.min(), reduced_x.max(), reduced_x.min(), reduced_x.max()])
            axes[i, 1].set_title(f'Density Matrix at t = {t:.1f}')
            plt.colorbar(im, ax=axes[i, 1])
            
        except Exception as e:
            print(f"Error at time {t}: {e}")
            axes[i, 0].text(0.5, 0.5, f"Error: {e}", ha='center', transform=axes[i, 0].transAxes)
            axes[i, 1].text(0.5, 0.5, f"Error: {e}", ha='center', transform=axes[i, 1].transAxes)
    
    axes[-1, 0].set_xlabel('Position (x)')
    axes[-1, 1].set_xlabel('Position (x)')
    axes[-1, 1].set_ylabel('Position (y)')
    
    plt.tight_layout()
    plt.show()

# 특정 시간에서의 진화 확인
# plot_evolution_at_times([0.0, 1.0, 2.0, 5.0])