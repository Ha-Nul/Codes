import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_causality():
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 1. 설정값 (임의의 단위)
    z_obs = 4.0       # 관측점 P의 높이 z
    c = 1.0           # 빛의 속도
    t = 6.0           # 시간 t (ct > z 가 되도록 설정)
    ct = c * t        # 정보가 도달할 수 있는 거리 (구의 반지름)
    
    # 교차 원판의 반지름 (r_max) 계산
    if ct > z_obs:
        r_max = np.sqrt(ct**2 - z_obs**2)
    else:
        r_max = 0
        print("신호가 아직 도달하지 않았습니다.")
        return

    # 2. 전류판 (z=0 평면) 그리기 - 맥락을 위한 배경
    grid_range = ct * 1.2
    x = np.linspace(-grid_range, grid_range, 50)
    y = np.linspace(-grid_range, grid_range, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    ax.plot_surface(X, Y, Z, alpha=0.1, color='gray', label='Current Sheet (z=0)')

    # 3. 인과율의 구 (Causality Sphere) 그리기 - 중심 P(0,0,z), 반지름 ct
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    x_sph = ct * np.outer(np.cos(u), np.sin(v))
    y_sph = ct * np.outer(np.sin(u), np.sin(v))
    z_sph = ct * np.outer(np.ones(np.size(u)), np.cos(v)) + z_obs
    
    # 구는 내부가 보여야 하므로 wireframe으로 그립니다.
    ax.plot_wireframe(x_sph, y_sph, z_sph, color='blue', alpha=0.2, linewidth=0.5)

    # 4. 적분 영역 (Integration Disk) 그리기 - 가장 중요한 부분!
    # 구와 평면이 만나는 빨간색 원판
    theta_disk = np.linspace(0, 2*np.pi, 100)
    r_disk = np.linspace(0, r_max, 50)
    R, THETA = np.meshgrid(r_disk, theta_disk)
    X_disk = R * np.cos(THETA)
    Y_disk = R * np.sin(THETA)
    Z_disk = np.zeros_like(X_disk)
    
    ax.plot_surface(X_disk, Y_disk, Z_disk, color='red', alpha=0.6)

    # 5. 기하학적 요소 표시 (점 P, 높이 z, 반지름 ct, r_max)
    
    # 관측점 P
    ax.scatter([0], [0], [z_obs], color='black', s=100, label='Observer P(0,0,z)')
    
    # 높이 z 선 (점선)
    ax.plot([0, 0], [0, 0], [0, z_obs], color='black', linestyle='--')
    ax.text(0, 0, z_obs/2, '  z', color='black')

    # 빗변 ct 선 (점선) - 빨간 원판의 가장자리로 연결
    ax.plot([0, r_max], [0, 0], [z_obs, 0], color='blue', linestyle='--', linewidth=2)
    ax.text(r_max/2, 0, z_obs/2, ' ct', color='blue', fontweight='bold')

    # 밑변 r_max 선 (실선)
    ax.plot([0, r_max], [0, 0], [0, 0], color='red', linewidth=2)
    ax.text(r_max/2, 0, -0.5, r' $\r_{max}$', color='red', fontweight='bold', fontsize=12)

    # 축 설정
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Causality Geometry: Integration Area (Red Disk)\nct = {ct}, z = {z_obs}, r_max = {r_max:.2f}')
    
    # 비율 유지 (중요: 기하학적 왜곡 방지)
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), z_sph.max()]).max() / 2.0
    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = z_sph.max() * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(0, z_sph.max())

    plt.show()

if __name__ == "__main__":
    visualize_causality()