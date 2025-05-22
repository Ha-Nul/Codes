import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

class QuantumSimulator:
    def __init__(self, dimension=2, hamiltonian=None, gamma=0.1):
        """
        양자 시스템 시뮬레이터 초기화
        
        Args:
            dimension: 힐버트 공간의 차원
            hamiltonian: 시스템의 해밀토니안 (None인 경우 기본값 생성)
            gamma: 완화 계수 (dissipation rate)
        """
        self.dimension = dimension
        
        # 기본 해밀토니안 생성 (만약 제공되지 않은 경우)
        if hamiltonian is None:
            # 랜덤 에르미트 행렬 생성
            h = np.random.rand(dimension, dimension) + 1j * np.random.rand(dimension, dimension)
            self.hamiltonian = (h + h.conj().T) / 2
        else:
            self.hamiltonian = hamiltonian
            
        self.gamma = gamma
        
        # 초기 밀도행렬 생성 (기본값: |0⟩⟨0|)
        self.rho = np.zeros((dimension, dimension), dtype=complex)
        self.rho[0, 0] = 1.0
        
        # 깁스 상태 계산 (thermal equilibrium state)
        self._compute_gibbs_state()
        
        # Lindblad 연산자들 생성
        self._create_lindblad_operators()
        
    def _compute_gibbs_state(self, beta=1.0):
        """깁스 상태 계산 - 양자 상세 균형 조건의 기준점"""
        # e^(-βH) 계산
        eigenvalues, eigenvectors = np.linalg.eigh(self.hamiltonian)
        exp_term = np.zeros_like(self.hamiltonian)
        
        for i in range(self.dimension):
            exp_term += np.exp(-beta * eigenvalues[i]) * np.outer(
                eigenvectors[:, i], eigenvectors[:, i].conj()
            )
            
        # 정규화
        self.gibbs_state = exp_term / np.trace(exp_term)
        
    def _create_lindblad_operators(self):
        """상세 균형 조건을 만족하는 Lindblad 연산자 생성"""
        # 고유값 분해
        eigenvalues, eigenvectors = np.linalg.eigh(self.hamiltonian)
        
        # 전이 연산자들 생성 (quantum jumps)
        self.lindblad_ops = []
        
        for i in range(self.dimension):
            for j in range(i+1, self.dimension):
                # 에너지 차이
                omega = eigenvalues[j] - eigenvalues[i]
                
                # 전이 확률 (detailed balance)
                rate = np.exp(-omega/2) if omega > 0 else 1
                
                # Lindblad 연산자: |i⟩⟨j|
                L_ij = np.sqrt(rate) * np.outer(
                    eigenvectors[:, i], eigenvectors[:, j].conj()
                )
                
                # Hermitian conjugate: |j⟩⟨i|
                L_ji = np.sqrt(rate) * np.outer(
                    eigenvectors[:, j], eigenvectors[:, i].conj()
                )
                
                self.lindblad_ops.append(L_ij)
                self.lindblad_ops.append(L_ji)
    
    def set_initial_state(self, rho):
        """초기 밀도행렬 설정"""
        if rho.shape != (self.dimension, self.dimension):
            raise ValueError("밀도행렬의 차원이 일치하지 않습니다")
        self.rho = rho
        
    def lindblad_dissipator(self, rho):
        """Lindblad 형태의 소산자 계산"""
        dissipator = np.zeros((self.dimension, self.dimension), dtype=complex)
        
        for L in self.lindblad_ops:
            # D[L]ρ = LρL† - 1/2{L†L, ρ}
            LrhoL_dag = L @ rho @ L.conj().T
            L_dag_L = L.conj().T @ L
            anti_commutator = L_dag_L @ rho + rho @ L_dag_L
            
            dissipator += LrhoL_dag - 0.5 * anti_commutator
            
        return self.gamma * dissipator
    
    def quantum_master_equation(self, rho, t):
        """양자 마스터 방정식 계산"""
        # Liouville-von Neumann 항: -i[H, ρ]
        commutator = -1j * (self.hamiltonian @ rho - rho @ self.hamiltonian)
        
        # Lindblad 소산자 항
        dissipator = self.lindblad_dissipator(rho)
        
        # 전체 시간 발전
        return commutator + dissipator
    
    def evolve(self, t_max, dt, return_times=False):
        """밀도행렬 시간 발전"""
        steps = int(t_max / dt)
        rho_history = np.zeros((steps, self.dimension, self.dimension), dtype=complex)
        times = np.linspace(0, t_max, steps)
        
        rho_history[0] = self.rho
        
        # 룽게-쿠타 4차 방법으로 시간 발전
        for i in range(1, steps):
            k1 = self.quantum_master_equation(self.rho, times[i-1])
            k2 = self.quantum_master_equation(self.rho + 0.5*dt*k1, times[i-1] + 0.5*dt)
            k3 = self.quantum_master_equation(self.rho + 0.5*dt*k2, times[i-1] + 0.5*dt)
            k4 = self.quantum_master_equation(self.rho + dt*k3, times[i-1] + dt)
            
            self.rho = self.rho + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
            
            # 수치적 오류로 인한 미세한 에르미트성 깨짐 보정
            self.rho = (self.rho + self.rho.conj().T) / 2
            
            # 트레이스가 1이 되도록 보정
            self.rho = self.rho / np.trace(self.rho)
            
            rho_history[i] = self.rho
            
        if return_times:
            return rho_history, times
        return rho_history
    
    def check_detailed_balance(self):
        """양자 상세 균형 조건 검증"""
        # 깁스 상태에서 시간 미분이 0에 가까운지 확인
        derivative = self.quantum_master_equation(self.gibbs_state, 0)
        error = np.linalg.norm(derivative)
        
        print(f"깁스 상태에서의 시간 미분 노름: {error:.6e}")
        if error < 1e-10:
            print("양자 상세 균형 조건이 만족됩니다.")
        else:
            print("양자 상세 균형 조건이 완벽히 만족되지 않습니다.")
        
        return error
    
    def visualize_matrix(self, rho, ax=None, title="밀도행렬"):
        """밀도행렬 시각화"""
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
            
        # 실수부와 허수부
        real_part = np.real(rho)
        imag_part = np.imag(rho)
        
        # 색상맵 범위 설정
        max_val = max(np.max(np.abs(real_part)), np.max(np.abs(imag_part)))
        
        # 실수부 시각화
        im1 = ax.imshow(real_part, cmap='RdBu_r', vmin=-max_val, vmax=max_val)
        
        # 행렬 값 표시
        for i in range(self.dimension):
            for j in range(self.dimension):
                text = f"{real_part[i, j]:.2f}"
                if imag_part[i, j] >= 0:
                    text += f"+{imag_part[i, j]:.2f}j"
                else:
                    text += f"{imag_part[i, j]:.2f}j"
                ax.text(j, i, text, ha='center', va='center', 
                        color='black' if abs(real_part[i, j]) < 0.5*max_val else 'white',
                        fontsize=9)
        
        ax.set_title(title)
        return im1
    
    def animate_evolution(self, rho_history, times, interval=100):
        """밀도행렬 진화 애니메이션"""
        fig, axs = plt.subplots(1, 2, figsize=(16, 6))
        
        # 밀도행렬 실수부 시각화
        max_val = np.max(np.abs(np.real(rho_history)))
        im1 = axs[0].imshow(np.real(rho_history[0]), cmap='RdBu_r', 
                           vmin=-max_val, vmax=max_val)
        axs[0].set_title("밀도행렬 실수부")
        plt.colorbar(im1, ax=axs[0])
        
        # 대각성분 (점유수) 시각화
        populations = np.real([np.diag(rho) for rho in rho_history])
        lines = []
        for i in range(self.dimension):
            line, = axs[1].plot(times, populations[:, i], 
                               label=f"상태 {i}", lw=2)
            lines.append(line)
        
        axs[1].set_title("상태 점유수 (대각성분)")
        axs[1].set_xlabel("시간")
        axs[1].set_ylabel("확률")
        axs[1].legend()
        axs[1].grid(True)
        
        # 시간 표시
        time_text = axs[0].text(0.02, 0.95, '', transform=axs[0].transAxes)
        
        def update(frame):
            # 밀도행렬 업데이트
            im1.set_array(np.real(rho_history[frame]))
            time_text.set_text(f'시간: {times[frame]:.2f}')
            
            # 대각성분 하이라이트
            for i, line in enumerate(lines):
                line.set_data(times[:frame+1], populations[:frame+1, i])
            
            axs[1].relim()
            axs[1].autoscale_view()
            
            return [im1, time_text] + lines
        
        ani = FuncAnimation(fig, update, frames=len(times),
                           interval=interval, blit=True)
        
        plt.tight_layout()
        return ani
    
    def bloch_vector(self, rho):
        """2차원 시스템에서 블로흐 벡터 계산"""
        if self.dimension != 2:
            raise ValueError("블로흐 벡터는 2차원 시스템에서만 정의됩니다")
        
        # 파울리 행렬
        sigma_x = np.array([[0, 1], [1, 0]])
        sigma_y = np.array([[0, -1j], [1j, 0]])
        sigma_z = np.array([[1, 0], [0, -1]])
        
        # 블로흐 벡터 성분 계산
        x = np.real(np.trace(rho @ sigma_x))
        y = np.real(np.trace(rho @ sigma_y))
        z = np.real(np.trace(rho @ sigma_z))
        
        return np.array([x, y, z])
    
    def animate_bloch_sphere(self, rho_history, times, interval=100):
        """블로흐 구 애니메이션 (2차원 시스템용)"""
        if self.dimension != 2:
            print("블로흐 구 애니메이션은 2차원 시스템에서만 가능합니다")
            return None
        
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # 블로흐 구 그리기
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = 0.99 * np.outer(np.cos(u), np.sin(v))
        y = 0.99 * np.outer(np.sin(u), np.sin(v))
        z = 0.99 * np.outer(np.ones(np.size(u)), np.cos(v))
        
        # 약간 투명한 블로흐 구
        ax.plot_surface(x, y, z, color='b', alpha=0.1)
        
        # 축 그리기
        ax.plot([-1, 1], [0, 0], [0, 0], 'k-', lw=1, alpha=0.5)  # x축
        ax.plot([0, 0], [-1, 1], [0, 0], 'k-', lw=1, alpha=0.5)  # y축
        ax.plot([0, 0], [0, 0], [-1, 1], 'k-', lw=1, alpha=0.5)  # z축
        
        # 레이블
        ax.text(1.1, 0, 0, r'$x$', fontsize=15)
        ax.text(0, 1.1, 0, r'$y$', fontsize=15)
        ax.text(0, 0, 1.1, r'$z$', fontsize=15)
        
        # 블로흐 벡터들 계산
        bloch_vectors = np.array([self.bloch_vector(rho) for rho in rho_history])
        
        # 초기 블로흐 벡터 플롯
        line, = ax.plot([], [], [], 'r-', lw=2, alpha=0.7)
        point, = ax.plot([], [], [], 'ro', ms=10)
        
        # 깁스 상태의 블로흐 벡터
        gibbs_bloch = self.bloch_vector(self.gibbs_state)
        ax.plot([0, gibbs_bloch[0]], [0, gibbs_bloch[1]], 
                [0, gibbs_bloch[2]], 'g--', lw=2, alpha=0.7)
        ax.plot([gibbs_bloch[0]], [gibbs_bloch[1]], [gibbs_bloch[2]], 
                'go', ms=10, label='깁스 상태')
        
        # 시간 표시
        time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes)
        
        # 축 설정
        ax.set_box_aspect([1, 1, 1]
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])
        ax.set_zlim([-1.1, 1.1])
        ax.set_title("블로흐 구에서의 양자 상태 진화", fontsize=16)
        
        def init():
            line.set_data([], [])
            line.set_3d_properties([])
            point.set_data([], [])
            point.set_3d_properties([])
            time_text.set_text('')
            return line, point, time_text
        
        def update(frame):
            # 경로 업데이트
            x_data = bloch_vectors[:frame+1, 0]
            y_data = bloch_vectors[:frame+1, 1]
            z_data = bloch_vectors[:frame+1, 2]
            
            line.set_data(x_data, y_data)
            line.set_3d_properties(z_data)
            
            # 현재 점 업데이트
            point.set_data([bloch_vectors[frame, 0]], [bloch_vectors[frame, 1]])
            point.set_3d_properties([bloch_vectors[frame, 2]])
            
            # 시간 표시 업데이트
            time_text.set_text(f'시간: {times[frame]:.2f}')
            
            return line, point, time_text
        
        ani = FuncAnimation(fig, update, frames=len(times), 
                           init_func=init, interval=interval, blit=True)
        
        ax.legend()
        plt.tight_layout()
        return ani

# 메인 실행 코드
if __name__ == "__main__":
    # 2차원 시스템 시뮬레이션 (쿼비트)
    dim = 2
    
    # 해밀토니안 설정 (예: 에너지 차이가 있는 2준위 시스템)
    H = np.array([
        [0.848, np.sqrt(0.152*0.848)],
        [np.sqrt(0.152*0.848), 0.152]
    ])
    
    # 시뮬레이터 초기화
    sim = QuantumSimulator(dimension=dim, hamiltonian=H, gamma=0.1)
    
    # 양자 상세 균형 조건 확인
    sim.check_detailed_balance()
    
    # 초기 상태 설정 (예: |0⟩ 상태)
    rho_init = np.zeros((dim, dim), dtype=complex)
    rho_init[0, 0] = 0.848
    rho_init[0, 1] = 0.359j
    rho_init[1, 0] = -0.359j
    rho_init[1, 1] = 0.152
    sim.set_initial_state(rho_init)
    
    # 시간 발전 시뮬레이션
    t_max = 50.0
    dt = 0.1
    rho_history, times = sim.evolve(t_max, dt, return_times=True)
    
    # 시각화
    plt.figure(figsize=(12, 8))
    
    # 초기 상태와 최종 상태 비교
    plt.subplot(2, 3, 1)
    sim.visualize_matrix(rho_history[0], plt.gca(), "초기 밀도행렬")
    
    plt.subplot(2, 3, 2)
    sim.visualize_matrix(rho_history[-1], plt.gca(), "최종 밀도행렬")
    
    plt.subplot(2, 3, 3)
    sim.visualize_matrix(sim.gibbs_state, plt.gca(), "깁스 상태 (평형)")
    
    # 대각성분(점유수) 진화 시각화
    plt.subplot(2, 1, 2)
    populations = np.real([np.diag(rho) for rho in rho_history])
    for i in range(dim):
        plt.plot(times, populations[:, i], label=f"상태 {i}", lw=2)
        
    plt.xlabel("시간")
    plt.ylabel("점유 확률")
    plt.title("상태 점유수 진화")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    
    # 애니메이션 생성
    ani = sim.animate_evolution(rho_history, times, interval=50)
    
    # 블로흐 구 애니메이션 (2차원 시스템의 경우)
    if dim == 2:
        bloch_ani = sim.animate_bloch_sphere(rho_history, times, interval=50)
    
    # 애니메이션 저장 (선택사항)
    # ani.save('density_matrix_evolution.mp4', writer='ffmpeg', fps=20)
    # if dim == 2:
    #     bloch_ani.save('bloch_sphere_evolution.mp4', writer='ffmpeg', fps=20)
    
    plt.show()