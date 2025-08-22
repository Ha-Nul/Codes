import contextlib
import time
import qutip as qt
from qutip.animation import anim_wigner


import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from qutip import (
    about, basis, brmesolve, destroy, expect, liouvillian,
    qeye, sigmax, sigmaz, spost, spre, tensor
)
from qutip.core.environment import (
    DrudeLorentzEnvironment, ExponentialBosonicEnvironment, system_terminator
)
from qutip.solver.heom import HEOMSolver, HSolverDL

#%matplotlib inline

def cot(x):
    """Vectorized cotangent of x."""
    return 1.0 / np.tan(x)

def dl_matsubara_params(lam, gamma, T, nk):
    """Calculation of the real and imaginary expansions of the Drude-Lorenz
    correlation functions.
    """
    ckAR = [lam * gamma * cot(gamma / (2 * T))]
    ckAR.extend(
        8 * lam * gamma * T * np.pi * k * T / (
            (2 * np.pi * k * T) ** 2 - gamma**2
        )
        for k in range(1, nk + 1)
    )
    vkAR = [gamma]
    vkAR.extend(2 * np.pi * k * T for k in range(1, nk + 1))

    ckAI = [lam * gamma * (-1.0)]
    vkAI = [gamma]

    return ckAR, vkAR, ckAI, vkAI

def plot_result_expectations(plots, axes=None):
    """Plot the expectation values of operators as functions of time.

    Each plot in plots consists of (solver_result, measurement_operation,
    color, label).
    """
    if axes is None:
        fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8, 8))
        fig_created = True
    else:
        fig = None
        fig_created = False

    # add kw arguments to each plot if missing
    plots = [p if len(p) == 5 else p + ({},) for p in plots]
    for result, m_op, color, label, kw in plots:
        exp = np.real(expect(result.states, m_op))
        kw.setdefault("linewidth", 2)
        axes.plot(result.times, exp, color, label=label, **kw)

    if fig_created:
        axes.legend(loc=0, fontsize=12)
        axes.set_xlabel("t", fontsize=28)

    return fig

@contextlib.contextmanager
def timer(label):
    """Simple utility for timing functions:

    with timer("name"):
        ... code to time ...
    """
    start = time.time()
    yield
    end = time.time()
    print(f"{label}: {end - start}")

# Default solver options:

default_options = {
    "nsteps": 1500,
    "store_states": True,
    "rtol": 1e-12,
    "atol": 1e-12,
    "method": "vern9",
    "progress_bar": "enhanced",
}

# Defining the system Hamiltonian
eps = 0.5  # Energy of the 2-level system.
Del = 1.0  # Tunnelling term
Hsys = 0.5 * eps * sigmaz() + 0.5 * Del * sigmax()

# Initial state of the system.
rho0 = basis(2, 0) * basis(2, 0).dag()

# System-bath coupling (Drude-Lorentz spectral density)
Q = sigmaz()  # coupling operator

# Bath properties:
gamma = 0.5  # cut off frequency
lam = 0.1  # coupling strength
T = 0.5
beta = 1.0 / T

#################################### HEOM parameters for Calculation time ###############################
NC = 20  # cut off parameter for the bath
Nk = 8  # terms in the Matsubara expansion of the correlation function
#################################### HEOM parameters for Calculation time ###############################

# Times to solve for
tlist = np.linspace(0, 50, 1000)

# Define some operators with which we will measure the system
# 1,1 element of density matrix - corresonding to groundstate
P11p = basis(2, 0) * basis(2, 0).dag()
P22p = basis(2, 1) * basis(2, 1).dag()
# 1,2 element of density matrix  - corresonding to coherence
P12p = basis(2, 0) * basis(2, 1).dag()

def plot_spectral_density():
    """Plot the Drude-Lorentz spectral density"""
    w = np.linspace(0, 5, 1000)
    J = w * 2 * lam * gamma / (gamma**2 + w**2)

    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8, 8))
    axes.plot(w, J, "r", linewidth=2)
    axes.set_xlabel(r"$\omega$", fontsize=28)
    axes.set_ylabel(r"J", fontsize=28)


plot_spectral_density()

ckAR, vkAR, ckAI, vkAI = dl_matsubara_params(nk=Nk, lam=lam, gamma=gamma, T=T)

options = {**default_options}

with timer("RHS construction time"):
    env = ExponentialBosonicEnvironment(ckAR, vkAR, ckAI, vkAI)
    HEOMMats = HEOMSolver(Hsys, (env, Q), NC, options=options)

with timer("ODE solver time"):
    resultMats = HEOMMats.run(rho0, tlist)




'''
plot_result_expectations(
    [
        (resultMats, P11p, "b", "P11 Mats"),
        (resultMats, P12p, "r", "P12 Mats"),
    ]
)
'''

system_states = resultMats.states

export_arr = np.array([state.full().flatten() for state in system_states])
combined = np.column_stack((tlist, export_arr))
np.savetxt('./testdataNc20Nk8.txt', combined)


# HEOM 결과로부터 위그너 애니메이션 생성
'''
xvec = np.linspace(-5, 5, 100)  # 해상도를 낮춰서 빠르게
pvec = np.linspace(-5, 5, 100)

fig, ani = anim_wigner(system_states, xvec=xvec, yvec=pvec, 
                      cmap='RdBu', colorbar=True)

plt.show()
'''
