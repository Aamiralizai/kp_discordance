"""Core multi-scale ODE: enzyme -> periplasm -> bulk depletion -> population.

Importable version of the model validated in experiments/toy_model.py.
On real data, replace `genotype_to_params` with the output of
src/kinetics/param_map.py (constants fixed from data/raw/enzyme_kinetics/).
"""
import numpy as np
from scipy.integrate import solve_ivp

# --- default PD / population constants (override from config for real fits) ---
DEFAULTS = dict(
    PSI=1.4, NMAX=1e9, EMAX=4.0, EC50=1.0, HILL=2.0,
    N0_STD=5e5, KM_HYD=50.0, KHYD=0.9, T_END=24.0,
)

CONC_LADDER = np.array([0.25 * 2**k for k in range(14)])  # 0.25 .. 2048 mg/L


def periplasmic_fraction(E, P, c_enz=2.0, efflux_mult=1.0):
    """kappa: fraction of external drug reaching the PBP target.
    Low P AND high enzyme jointly collapse kappa -> super-additive synergy."""
    kappa = P / (P + c_enz * E)
    return kappa * efflux_mult


def _odes(t, y, E, kappa, p):
    A, N = max(y[0], 0.0), max(y[1], 1.0)
    A_peri = kappa * A
    kill = p["EMAX"] * A_peri**p["HILL"] / (p["EC50"]**p["HILL"] + A_peri**p["HILL"])
    dN = p["PSI"] * N * (1 - N / p["NMAX"]) - kill * N
    dA = -p["KHYD"] * E * (N / 1e8) * A / (p["KM_HYD"] + A)   # collective depletion
    return [dA, dN]


def simulate(A0, E, kappa, N0=None, params=None, n_eval=60):
    p = {**DEFAULTS, **(params or {})}
    N0 = p["N0_STD"] if N0 is None else N0
    t_eval = np.linspace(0, p["T_END"], n_eval)
    sol = solve_ivp(_odes, [0, p["T_END"]], [A0, N0], args=(E, kappa, p),
                    t_eval=t_eval, method="LSODA", rtol=1e-6, atol=1e-3)
    return sol.t, sol.y[0], sol.y[1]


def predict_MIC(E, kappa, N0=None, params=None):
    """Lowest ladder concentration preventing net growth over T_END."""
    p = {**DEFAULTS, **(params or {})}
    N0 = p["N0_STD"] if N0 is None else N0
    for A0 in CONC_LADDER:
        _, _, N = simulate(A0, E, kappa, N0=N0, params=p)
        if N[-1] <= N0:
            return A0
    return CONC_LADDER[-1] * 2  # censored high


if __name__ == "__main__":
    # sanity check: carbapenemase + porin loss should be super-additive
    for name, (E, P) in {
        "WT":                 (0.0, 1.0),
        "carbapenemase":      (1.0, 1.0),
        "porin_loss":         (0.0, 0.05),
        "carbapenemase+porin":(1.0, 0.05),
    }.items():
        kappa = periplasmic_fraction(E, P)
        print(f"{name:22s} MIC = {predict_MIC(E, kappa):7.1f} mg/L")
