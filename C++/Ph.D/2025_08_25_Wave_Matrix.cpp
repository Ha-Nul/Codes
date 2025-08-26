#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

using matrix_t = Eigen::MatrixXcd;
using state_t = std::vector<std::complex<double>>;
using namespace boost::numeric::odeint;

// 크로네커 곱
matrix_t kron(const matrix_t& A, const matrix_t& B) {
    matrix_t C(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }
    return C;
}

// 행렬을 std::vector로 변환
state_t vectorize_matrix(const matrix_t& mat) {
    state_t vec(mat.rows() * mat.cols());
    int idx = 0;
    for (int j = 0; j < mat.cols(); ++j) {
        for (int i = 0; i < mat.rows(); ++i) {
            vec[idx++] = mat(i, j);
        }
    }
    return vec;
}

// std::vector를 행렬로 변환
matrix_t devectorize_to_matrix(const state_t& vec, int rows, int cols) {
    matrix_t mat(rows, cols);
    int idx = 0;
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            mat(i, j) = vec[idx++];
        }
    }
    return mat;
}

// H_hat 행렬 생성 함수 (6 큐비트)
matrix_t create_H_hat_matrix() {
    /*-------------Qubit count: 6 qubits (2^6 = 64)------------*/
    int dim = 64;
    /*-----------------------------------*/
    matrix_t swap_matrix = matrix_t::Zero(dim, dim);
    int s_pos = 5, h_pos = 2;  // 6큐비트: 비트 위치 0-5 범위
    for (int i = 0; i < dim; ++i) {
        int s_bit = (i >> s_pos) & 1;
        int h_bit = (i >> h_pos) & 1;
        int j = i;
        if (s_bit != h_bit) {
            int mask = (1 << s_pos) | (1 << h_pos);
            j = i ^ mask;
        }
        swap_matrix(j, i) = 1.0;
    }
    return swap_matrix;
}

// M 연산자 행렬 생성 함수 (6 큐비트)
matrix_t create_M_matrix() {
    /*-------------Qubit count: 6 qubits (2^6 = 64)------------*/
    int dim = 64;
    /*-----------------------------------*/
    int d_system = 2;
    
    // 1. SWAP(S, P) 연산자 행렬 생성
    matrix_t swap_sp_matrix = matrix_t::Zero(dim, dim);
    int s_pos = 5, p_pos = 1;  // 6큐비트: 비트 위치 0-5 범위
    for (int i = 0; i < dim; ++i) {
        int s_bit = (i >> s_pos) & 1;
        int p_bit = (i >> p_pos) & 1;
        int j = i;
        if (s_bit != p_bit) {
            int mask = (1 << s_pos) | (1 << p_pos);
            j = i ^ mask;
        }
        swap_sp_matrix(j, i) = 1.0;
    }

    // 2. 얽힘 상태 프로젝터 연산자 행렬 생성
    int pq_dim = d_system * d_system;  // 4
    Eigen::VectorXcd psi_gamma = Eigen::VectorXcd::Zero(pq_dim);
    for (int i = 0; i < d_system; ++i) {
        psi_gamma(i * d_system + i) = 1.0 / std::sqrt(d_system);
    }
    
    matrix_t projector_pq = psi_gamma * psi_gamma.adjoint();
    matrix_t I_srh = matrix_t::Identity(16, 16);  // 6큐비트에서 나머지 부분: 64/4 = 16
    matrix_t projector_full_matrix = kron(I_srh, projector_pq);

    // 3. 최종 M 행렬 계산
    matrix_t M_matrix = (1.0 / std::sqrt(d_system)) * projector_full_matrix * swap_sp_matrix;
    return M_matrix;
}

// Boost.odeint를 위한 동역학 정의
struct LindbladSystem {
    const matrix_t& m_L;
    LindbladSystem(const matrix_t& L) : m_L(L) {}
    
    void operator()(const state_t& rho_vec, state_t& drho_vec_dt, double /* t */) const {
        // std::vector → Eigen::VectorXcd 변환
        Eigen::VectorXcd rho_eigen = Eigen::Map<const Eigen::VectorXcd>(rho_vec.data(), rho_vec.size());
        
        // 행렬-벡터 곱셈
        /*--------------------------------------------------------------*/
        Eigen::VectorXcd result = m_L * rho_eigen;
        /*--------------------------------------------------------------*/
        // Eigen::VectorXcd → std::vector 변환
        drho_vec_dt.resize(result.size());
        for (int i = 0; i < result.size(); ++i) {
            drho_vec_dt[i] = result(i);
        }
    }
};

int main() {
    std::cout << "--- C++ Benchmark Start (6 Qubits) ---" << std::endl;

    try {
        // 1. 연산자 행렬 생성
        std::cout << "1. Creating operator matrices..." << std::endl;
        matrix_t H = create_H_hat_matrix();
        matrix_t M = create_M_matrix();
        std::cout << "H matrix size: " << H.rows() << "x" << H.cols() << std::endl;
        std::cout << "M matrix size: " << M.rows() << "x" << M.cols() << std::endl;

        // 2. 리우빌리안 초연산자 L 생성
        std::cout << "2. Creating Liouvillian super-operator..." << std::endl;
        /*----------------------Matrix Creation (6 qubits)----------------------------*/
        matrix_t I_64 = matrix_t::Identity(64, 64);  // 64x64 단위행렬
        /*----------------------Matrix Creation----------------------------*/
        matrix_t MtM = M.adjoint() * M;
        std::complex<double> I_im(0.0, -1.0);

        matrix_t L = I_im * (kron(I_64, H) - kron(H.transpose(), I_64)) +
                       (kron(M.conjugate(), M) - 
                        0.5 * kron(I_64, MtM) - 
                        0.5 * kron(MtM.transpose(), I_64));
        
        std::cout << "Liouvillian size: " << L.rows() << "x" << L.cols() << std::endl;

        // 3. 초기 상태 정의 및 벡터화
        std::cout << "3. Defining initial state..." << std::endl;
        /*----------------------Matrix Creation (6 qubits)----------------------------*/
        matrix_t rho_initial = matrix_t::Zero(64, 64);  // 64x64 밀도행렬
        /*----------------------Matrix Creation----------------------------*/
        rho_initial(0, 0) = 1.0;  // |000000><000000| 상태
        
        state_t rho_vec = vectorize_matrix(rho_initial);
        std::cout << "Initial state vector size: " << rho_vec.size() << std::endl;
        std::cout << "Initial trace: " << rho_initial.trace() << std::endl;

        // 4. Boost.odeint를 이용한 시뮬레이션
        std::cout << "4. Running ODE solver..." << std::endl;
        /*--------------------------time step adjust : calculation process -----------------------*/
        double t = 100;
        double delta = 0.01;
        int n_steps = 10000;
        /*--------------------------time step adjust : calculation process -----------------------*/
        
        LindbladSystem system(L);
        runge_kutta4<state_t> stepper;
        
        // 여러 스텝 실행 (원래 코드는 한 스텝만 했지만, 실제 시뮬레이션을 위해)
        for (int step = 0; step < 100; ++step) {  // 100 스텝만 실행
            stepper.do_step(system, rho_vec, t, delta);
            t += delta;
            if (step % 20 == 0) {
                std::cout << "Step " << step << " completed" << std::endl;
            }
        }
        
        std::cout << "After ODE steps, vector size: " << rho_vec.size() << std::endl;

        // 5. 결과 확인
        std::cout << "5. Simulation finished." << std::endl;
        matrix_t rho_final = devectorize_to_matrix(rho_vec, 64, 64);  // 64x64로 변경
        
        std::cout << "Final trace: " << rho_final.trace() << std::endl;
        std::cout << "Final state (0,0) element: " << rho_final(0, 0) << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "--- Benchmark Complete ---" << std::endl;
    return 0;
}
