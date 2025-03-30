import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# 상수 정의
m = 1.0  # 질량
k = 1.0  # 스프링 상수
E_total = 10.0  # 총 에너지 (임의 값)

def hamiltonian(x, y, z, px, py, pz):
    """
    해밀토니안 함수: 총 에너지 = 운동 에너지 + 퍼텐셜 에너지
    여기서는 3D 조화 진동자(harmonic oscillator) 모델을 사용
    """
    kinetic = (px**2 + py**2 + pz**2) / (2 * m)  # 운동 에너지 T = p^2/(2m)
    potential = 0.5 * k * (x**2 + y**2 + z**2)  # 퍼텐셜 에너지 V = (1/2)k*r^2
    return kinetic + potential

# 1. 위치 공간에서의 등위면(constant energy surface) - 구 형태
def plot_position_sphere():
    # 위치 공간에서, 퍼텐셜 에너지가 일정하다면 r^2 = 상수
    # 따라서 위치 좌표들은 구 위에 있어야 함
    
    # 운동 에너지가 0일 때 (정지해 있을 때) 가능한 최대 퍼텐셜 에너지 = 총 에너지
    # 따라서 x^2 + y^2 + z^2 = 2*E_total/k
    
    r_max = np.sqrt(2 * E_total / k)
    
    # 구 표면 생성
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r_max * np.outer(np.cos(u), np.sin(v))
    y = r_max * np.outer(np.sin(u), np.sin(v))
    z = r_max * np.outer(np.ones(np.size(u)), np.cos(v))
    
    # 플롯 생성
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot_surface(x, y, z, color='b', alpha=0.3)
    
    # 좌표축 설정
    max_range = r_max * 1.2
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Constraint Energy surface, in coordinate space (V = E_total)')
    
    plt.show()

# 2. 운동량 공간에서의 등위면
def plot_momentum_sphere():
    # 운동량 공간에서, 운동 에너지가 일정하다면 p^2 = 상수
    # 따라서 운동량 좌표들은 구 위에 있어야 함
    
    # 퍼텐셜 에너지가 0일 때 가능한 최대 운동 에너지 = 총 에너지
    # 따라서 px^2 + py^2 + pz^2 = 2*m*E_total
    
    p_max = np.sqrt(2 * m * E_total)
    
    # 구 표면 생성
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    px = p_max * np.outer(np.cos(u), np.sin(v))
    py = p_max * np.outer(np.sin(u), np.sin(v))
    pz = p_max * np.outer(np.ones(np.size(u)), np.cos(v))
    
    # 플롯 생성
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot_surface(px, py, pz, color='r', alpha=0.3)
    
    # 좌표축 설정
    max_range = p_max * 1.2
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    
    ax.set_xlabel('p_x')
    ax.set_ylabel('p_y')
    ax.set_zlabel('p_z')
    ax.set_title('Constraint Energy surface, in momentum space (T = E_total)')
    
    plt.show()

# 3. 위상 공간에서의 등위면 샘플링
def sample_phase_space_points(n_points=1000):
    """
    에너지 조건을 만족하는 위상 공간의 점들을 샘플링
    """
    points = []
    
    # 무작위 샘플링 시도
    count = 0
    while len(points) < n_points and count < n_points * 100:
        count += 1
        
        # 위치 범위 (위치 에너지가 총 에너지보다 작아야 함)
        r_max = np.sqrt(2 * E_total / k)
        
        # 무작위 위치 생성
        x = (2 * np.random.random() - 1) * r_max
        y = (2 * np.random.random() - 1) * r_max
        z = (2 * np.random.random() - 1) * r_max
        
        # 남은 에너지로 가능한 운동 에너지 계산
        V = 0.5 * k * (x**2 + y**2 + z**2)
        T_remaining = E_total - V
        
        # 에너지가 음수면 이 위치는 불가능
        if T_remaining <= 0:
            continue
        
        # 가능한 최대 운동량
        p_max = np.sqrt(2 * m * T_remaining)
        
        # 운동량 공간에서 무작위 방향 선택 (구 표면 위의 점)
        phi = 2 * np.pi * np.random.random()
        cos_theta = 2 * np.random.random() - 1
        sin_theta = np.sqrt(1 - cos_theta**2)
        
        px = p_max * sin_theta * np.cos(phi)
        py = p_max * sin_theta * np.sin(phi)
        pz = p_max * cos_theta
        
        # 해밀토니안 검증 (부동소수점 오차를 고려하여 근사값 확인)
        H = hamiltonian(x, y, z, px, py, pz)
        if abs(H - E_total) < 1e-10:
            points.append((x, y, z, px, py, pz))
    
    return np.array(points)

# 4. 위상 공간에서의 투영 시각화
def plot_phase_space_projections():
    points = sample_phase_space_points(3000)
    
    if len(points) == 0:
        print("Cannot find availiable phase point.")
        return
    
    fig = plt.figure(figsize=(15, 10))
    
    # 위치 공간 투영 (x, y, z)
    ax1 = fig.add_subplot(231, projection='3d')
    ax1.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', s=1, alpha=0.5)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_title('projection of coordinate space on phase space')
    
    # 운동량 공간 투영 (px, py, pz)
    ax2 = fig.add_subplot(232, projection='3d')
    ax2.scatter(points[:, 3], points[:, 4], points[:, 5], c='r', s=1, alpha=0.5)
    ax2.set_xlabel('p_x')
    ax2.set_ylabel('p_y')
    ax2.set_zlabel('p_z')
    ax2.set_title('projection of momentum space on phase space')
    
    # x-px 위상 공간 투영
    ax3 = fig.add_subplot(233)
    ax3.scatter(points[:, 0], points[:, 3], c=points[:, 2], s=1, alpha=0.5)
    ax3.set_xlabel('x')
    ax3.set_ylabel('p_x')
    ax3.set_title('x-p_x phase space projection')
    
    # y-py 위상 공간 투영
    ax4 = fig.add_subplot(234)
    ax4.scatter(points[:, 1], points[:, 4], c=points[:, 2], s=1, alpha=0.5)
    ax4.set_xlabel('y')
    ax4.set_ylabel('p_y')
    ax4.set_title('y-p_y phase space projection')
    
    # z-pz 위상 공간 투영
    ax5 = fig.add_subplot(235)
    ax5.scatter(points[:, 2], points[:, 5], c=points[:, 0], s=1, alpha=0.5)
    ax5.set_xlabel('z')
    ax5.set_ylabel('p_z')
    ax5.set_title('z-p_z phase space projection')
    
    # 에너지 분포
    ax6 = fig.add_subplot(236)
    energies = [hamiltonian(*p) for p in points]
    ax6.hist(energies, bins=50)
    ax6.set_xlabel('Energy')
    ax6.set_ylabel('Frequent')
    ax6.set_title('Energy distribution')
    
    plt.tight_layout()
    plt.show()

# 5. 에너지 등고선에 따른 위상 공간 매니폴드 시각화
def plot_energy_ellipsoids():
    """
    다양한 에너지 값에 대해 위치 및 운동량 공간에서의 등위면 시각화
    """
    fig = plt.figure(figsize=(15, 8))
    
    # 위치 공간 등위면 (다양한 에너지)
    ax1 = fig.add_subplot(121, projection='3d')
    
    # 다양한 에너지 값에 대해 등위면 그리기
    energies = [E_total * 0.2, E_total * 0.4, E_total * 0.6, E_total * 0.8, E_total]
    colors = ['blue', 'cyan', 'green', 'orange', 'red']
    
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    
    for i, E in enumerate(energies):
        r = np.sqrt(2 * E / k)
        x = r * np.outer(np.cos(u), np.sin(v))
        y = r * np.outer(np.sin(u), np.sin(v))
        z = r * np.outer(np.ones(np.size(u)), np.cos(v))
        
        ax1.plot_surface(x, y, z, color=colors[i], alpha=0.2)
    
    r_max = np.sqrt(2 * E_total / k)
    max_range = r_max * 1.2
    ax1.set_xlim(-max_range, max_range)
    ax1.set_ylim(-max_range, max_range)
    ax1.set_zlim(-max_range, max_range)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_title('위치 공간 에너지 등위면')
    
    # 운동량 공간 등위면 (다양한 에너지)
    ax2 = fig.add_subplot(122, projection='3d')
    
    for i, E in enumerate(energies):
        p = np.sqrt(2 * m * E)
        px = p * np.outer(np.cos(u), np.sin(v))
        py = p * np.outer(np.sin(u), np.sin(v))
        pz = p * np.outer(np.ones(np.size(u)), np.cos(v))
        
        ax2.plot_surface(px, py, pz, color=colors[i], alpha=0.2)
    
    p_max = np.sqrt(2 * m * E_total)
    max_range = p_max * 1.2
    ax2.set_xlim(-max_range, max_range)
    ax2.set_ylim(-max_range, max_range)
    ax2.set_zlim(-max_range, max_range)
    ax2.set_xlabel('p_x')
    ax2.set_ylabel('p_y')
    ax2.set_zlabel('p_z')
    ax2.set_title('운동량 공간 에너지 등위면')
    
    plt.tight_layout()
    plt.show()

# 6. 조화 진동자 궤적 시뮬레이션 (시간에 따른 위상 공간 궤적)
def simulate_harmonic_oscillator_trajectory():
    """
    3D 조화 진동자의 위상 공간 궤적을 시뮬레이션
    """
    # 초기 조건 설정
    x0, y0, z0 = 1.0, 1.0, 1.0
    px0, py0, pz0 = 0.0, 0.0, 0.0
    
    # 초기 에너지 계산 및 조정
    E0 = hamiltonian(x0, y0, z0, px0, py0, pz0)
    scale_factor = np.sqrt(E_total / E0)
    
    # 초기 조건 조정 (에너지를 E_total로 맞춤)
    x0 *= scale_factor
    y0 *= scale_factor
    z0 *= scale_factor
    
    # 시간 단계
    dt = 0.05
    t_max = 50
    steps = int(t_max / dt)
    
    # 궤적 저장 배열
    trajectory = np.zeros((steps, 6))
    trajectory[0] = [x0, y0, z0, px0, py0, pz0]
    
    # 해밀턴 방정식에 따른 시간 발전
    for i in range(1, steps):
        x, y, z, px, py, pz = trajectory[i-1]
        
        # 해밀턴 방정식: 
        # dx/dt = ∂H/∂px, dy/dt = ∂H/∂py, dz/dt = ∂H/∂pz
        # dpx/dt = -∂H/∂x, dpy/dt = -∂H/∂y, dpz/dt = -∂H/∂z
        
        # 위치 업데이트
        new_x = x + px/m * dt
        new_y = y + py/m * dt
        new_z = z + pz/m * dt
        
        # 운동량 업데이트
        new_px = px - k * x * dt
        new_py = py - k * y * dt
        new_pz = pz - k * z * dt
        
        trajectory[i] = [new_x, new_y, new_z, new_px, new_py, new_pz]
    
    # 시간에 따른 에너지 변화 확인
    energies = np.array([hamiltonian(*trajectory[i]) for i in range(steps)])
    energy_error = (energies - E_total) / E_total
    
    # 궤적 시각화
    fig = plt.figure(figsize=(15, 10))
    
    # 위치 공간 궤적
    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2], 'b-', lw=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_title('Trajectory of coordinate space')
    
    # 운동량 공간 궤적
    ax2 = fig.add_subplot(222, projection='3d')
    ax2.plot(trajectory[:, 3], trajectory[:, 4], trajectory[:, 5], 'r-', lw=1)
    ax2.set_xlabel('p_x')
    ax2.set_ylabel('p_y')
    ax2.set_zlabel('p_z')
    ax2.set_title('Trajectory of momentum space')
    
    # x-px 위상 평면 궤적
    ax3 = fig.add_subplot(223)
    ax3.plot(trajectory[:, 0], trajectory[:, 3], 'g-', lw=1)
    ax3.set_xlabel('x')
    ax3.set_ylabel('p_x')
    ax3.set_title('x-p_x Trajectory of Phase space projection')
    
    # 에너지 오차
    ax4 = fig.add_subplot(224)
    ax4.plot(np.arange(steps) * dt, energy_error, 'k-', lw=1)
    ax4.set_xlabel('Time')
    ax4.set_ylabel('Relative variation of energy')
    ax4.set_title('Variation of energy conservation')
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()

# 모든 시각화 함수 실행
def run_all_visualizations():
    print("1. Constant Energy in coordinate space")
    plot_position_sphere()
    
    print("2. Constant Energy in momentum space")
    plot_momentum_sphere()
    
    print("3. Projection in phase space")
    plot_phase_space_projections()
    
    print("4. 다양한 에너지의 등위면")
    plot_energy_ellipsoids()
    
    print("5. Simulation of Harmonic oscillator trajectory")
    simulate_harmonic_oscillator_trajectory()

# 실행
if __name__ == "__main__":
    run_all_visualizations()