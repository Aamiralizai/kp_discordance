"""
TOY MODEL — Multi-scale mechanistic genotype -> phenotype map for
carbapenem resistance in Klebsiella pneumoniae.

Purpose: proof-of-concept on SYNTHETIC data where we know the ground truth.
It demonstrates the pipeline's central claim before touching real data:

    mutations  ->  kinetic parameters  ->  multi-scale ODE  ->  phenotype (MIC / time-kill)

and the headline test:

    A MECHANISTIC model that COMPOSES effects generalizes to unseen
    mutation COMBINATIONS, where a CORRELATIONAL ML model (trained where those
    mutations never co-occur) fails to recover their super-additive interaction.

This is a deliberately simplified system. Numbers are illustrative, not fitted
to real isolates. It exists to check that the logic produces the expected
qualitative outcomes.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt
import itertools

rng = np.random.default_rng(7)

# ----------------------------------------------------------------------
# 1. FIXED "TRUE" BIOLOGY  (the synthetic ground truth generator)
# ----------------------------------------------------------------------
# Population / PD parameters (per hour, mg/L, CFU/mL)
PSI   = 1.4      # max growth rate
NMAX  = 1e9      # carrying capacity
EMAX  = 4.0      # max kill rate
EC50  = 1.0      # periplasmic conc for half-max kill (mg/L)
HILL  = 2.0
N0_STD = 5e5     # standard inoculum (CLSI-like)
KM_HYD = 50.0    # hydrolysis Km (mg/L)
KHYD   = 0.9     # collective hydrolysis scale
T_END  = 24.0    # hours

# genotype feature order
FEATURES = ["carbapenemase", "high_copy", "porin_loss", "efflux_up"]

def genotype_to_params(g, theta=None):
    """Map a binary genotype vector to mechanistic kinetic parameters.

    theta lets a *calibrated* mechanistic model tune a few scalars; when
    theta is None we use the TRUE generative values.
    """
    if theta is None:
        c_enz, P_loss, E_base, copy_mult, efflux_mult = 2.0, 0.05, 1.0, 3.0, 0.5
    else:
        c_enz, P_loss, E_base, copy_mult, efflux_mult = theta

    carb, copy_, porin, efflux = g
    E = E_base * carb                      # enzyme level (0 if no carbapenemase)
    if carb and copy_:
        E *= copy_mult                     # copy number amplifies enzyme
    P = P_loss if porin else 1.0           # porin permeability
    # periplasmic exposure fraction: low P AND high enzyme -> strong synergy
    kappa = P / (P + c_enz * E)
    if efflux:
        kappa *= efflux_mult               # efflux removes drug from target
    return E, kappa

def odes(t, y, A_ext_dep, E, kappa):
    A, N = y
    A = max(A, 0.0)
    N = max(N, 1.0)
    A_peri = kappa * A
    kill = EMAX * A_peri**HILL / (EC50**HILL + A_peri**HILL)
    dN = PSI * N * (1.0 - N / NMAX) - kill * N
    # collective enzymatic depletion scales with population size (inoculum effect)
    dA = -A_ext_dep * E * (N / 1e8) * A / (KM_HYD + A)
    return [dA, dN]

def simulate(A0, E, kappa, N0=N0_STD, t_end=T_END, n_eval=60):
    t_eval = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(odes, [0, t_end], [A0, N0], args=(KHYD, E, kappa),
                    t_eval=t_eval, method="LSODA", rtol=1e-6, atol=1e-3)
    return sol.t, sol.y[0], sol.y[1]

# doubling-dilution concentration ladder for MIC
CONC_LADDER = np.array([0.25 * 2**k for k in range(14)])   # 0.25 .. 2048 mg/L

def predict_MIC(E, kappa, N0=N0_STD):
    """Lowest concentration on the ladder that prevents net growth over T_END."""
    for A0 in CONC_LADDER:
        _, _, N = simulate(A0, E, kappa, N0=N0)
        if N[-1] <= N0:            # net static/cidal
            return A0
    return CONC_LADDER[-1] * 2     # off-scale (censored high)

def mic_of_genotype(g, theta=None, N0=N0_STD):
    E, kappa = genotype_to_params(g, theta)
    return predict_MIC(E, kappa, N0=N0)

# ----------------------------------------------------------------------
# 2. GROUND-TRUTH PHENOTYPES for all 2^4 genotypes
# ----------------------------------------------------------------------
all_genotypes = [np.array(g) for g in itertools.product([0, 1], repeat=4)]
true_mic = {tuple(g): mic_of_genotype(g) for g in all_genotypes}

print("Genotype (carb,copy,porin,efflux) -> true MIC (mg/L)")
for g in all_genotypes:
    print(f"  {tuple(int(x) for x in g)} -> {true_mic[tuple(g)]:7.1f}")

# ----------------------------------------------------------------------
# 3. TRAIN / OOD SPLIT  (mimic lineage co-occurrence structure)
#    In "observed" data, porin_loss and carbapenemase NEVER co-occur.
#    The clinically dangerous convergent genotypes are therefore OOD.
# ----------------------------------------------------------------------
def is_ood(g):
    carb, copy_, porin, efflux = g
    return bool(carb and porin)          # the unseen convergent combination

train_g = [g for g in all_genotypes if not is_ood(g)]
ood_g   = [g for g in all_genotypes if is_ood(g)]

log2 = lambda x: np.log2(np.array(x))
Xtr = np.array([g for g in train_g]); ytr = log2([true_mic[tuple(g)] for g in train_g])
Xood = np.array([g for g in ood_g]);  yood = log2([true_mic[tuple(g)] for g in ood_g])

# ----------------------------------------------------------------------
# 4a. CORRELATIONAL ML MODELS  (trained on log2 MIC of training genotypes)
# ----------------------------------------------------------------------
rf = RandomForestRegressor(n_estimators=400, random_state=0).fit(Xtr, ytr)
ridge = Ridge(alpha=0.1).fit(Xtr, ytr)
rf_ood_pred, ridge_ood_pred = rf.predict(Xood), ridge.predict(Xood)

# ----------------------------------------------------------------------
# 4b. MECHANISTIC MODEL  (calibrate a few scalars on the SAME training MICs,
#     then predict OOD by COMPOSITION — same information as ML, right structure)
# ----------------------------------------------------------------------
def mech_loss(theta):
    pred = log2([mic_of_genotype(g, theta=theta) for g in train_g])
    return np.mean((pred - ytr) ** 2)

theta0 = np.array([2.0, 0.05, 1.0, 3.0, 0.5])
res = minimize(mech_loss, theta0, method="Nelder-Mead",
               options={"maxiter": 300, "xatol": 1e-3, "fatol": 1e-4})
theta_hat = res.x
mech_ood_pred = log2([mic_of_genotype(g, theta=theta_hat) for g in ood_g])

def mae(a, b): return float(np.mean(np.abs(np.array(a) - np.array(b))))

print("\n=== OOD generalization (log2 MIC MAE, doubling dilutions) ===")
print(f"  Random Forest : {mae(rf_ood_pred, yood):.2f}")
print(f"  Ridge (linear): {mae(ridge_ood_pred, yood):.2f}")
print(f"  Mechanistic   : {mae(mech_ood_pred, yood):.2f}")

# ----------------------------------------------------------------------
# 5. FIGURES
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 11})
C = {"wt": "#2b8cbe", "carb": "#238b45", "porin": "#d95f0e", "both": "#cb181d",
     "mech": "#333333", "rf": "#e6550d", "ridge": "#9e9ac8"}

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

# --- Fig 1: time-kill curves at a fixed clinically relevant dose ---
ax = axes[0]
dose = 16.0
scenarios = [
    ("WT (no mechanisms)",        np.array([0,0,0,0]), C["wt"]),
    ("Carbapenemase only",        np.array([1,0,0,0]), C["carb"]),
    ("Porin loss only",           np.array([0,0,1,0]), C["porin"]),
    ("Carbapenemase + porin loss",np.array([1,0,1,0]), C["both"]),
]
for label, g, col in scenarios:
    E, kappa = genotype_to_params(g)
    t, A, N = simulate(dose, E, kappa)
    ax.plot(t, np.log10(np.clip(N, 1, None)), color=col, lw=2.2, label=label)
ax.axhline(np.log10(N0_STD), ls=":", color="grey", lw=1)
ax.set_title(f"Time-kill at {dose:.0f} mg/L", fontweight="bold")
ax.set_xlabel("Time (h)"); ax.set_ylabel("log10 CFU/mL")
ax.legend(fontsize=8.5, frameon=False); ax.set_ylim(2, 9.5)

# --- Fig 2: inoculum effect (apparent MIC vs starting inoculum) ---
ax = axes[1]
inocula = np.array([1e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8])
for label, g, col in [("Carbapenemase", np.array([1,0,0,0]), C["carb"]),
                      ("Carbapenemase + porin loss", np.array([1,0,1,0]), C["both"])]:
    E, kappa = genotype_to_params(g)
    mics = [predict_MIC(E, kappa, N0=n0) for n0 in inocula]
    ax.plot(np.log10(inocula), np.log2(mics), "o-", color=col, lw=2.2, label=label)
ax.set_title("Inoculum effect (emergent)", fontweight="bold")
ax.set_xlabel("log10 starting inoculum (CFU/mL)")
ax.set_ylabel("log2 apparent MIC (mg/L)")
ax.legend(fontsize=8.5, frameon=False)

# --- Fig 3: OOD generalization — predicted vs true (the money plot) ---
ax = axes[2]
lim = [min(yood)-1, max(yood)+1]
ax.plot(lim, lim, "k--", lw=1, alpha=0.6)
ax.scatter(yood, mech_ood_pred, s=90, color=C["mech"], label=f"Mechanistic (MAE {mae(mech_ood_pred,yood):.1f})", zorder=3)
ax.scatter(yood, rf_ood_pred, s=90, color=C["rf"], marker="s", label=f"Random Forest (MAE {mae(rf_ood_pred,yood):.1f})", zorder=3)
ax.scatter(yood, ridge_ood_pred, s=90, color=C["ridge"], marker="^", label=f"Ridge (MAE {mae(ridge_ood_pred,yood):.1f})", zorder=3)
ax.set_title("Unseen convergent genotypes (OOD)", fontweight="bold")
ax.set_xlabel("True log2 MIC"); ax.set_ylabel("Predicted log2 MIC")
ax.legend(fontsize=8.5, frameon=False, loc="upper left")

for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("results/figures/toy_results.png", dpi=150, bbox_inches="tight")
print("\nSaved figure -> kp_toy_results.png")

# summary numbers for the writeup
print("\n=== KEY TAKEAWAYS (synthetic) ===")
gc = tuple(np.array([1,0,0,0])); gp = tuple(np.array([0,0,1,0])); gb = tuple(np.array([1,0,1,0])); gw = tuple(np.array([0,0,0,0]))
print(f"  MIC WT                     : {true_mic[gw]:.1f} mg/L")
print(f"  MIC carbapenemase only     : {true_mic[gc]:.1f} mg/L")
print(f"  MIC porin loss only        : {true_mic[gp]:.1f} mg/L")
print(f"  MIC BOTH (convergent, OOD) : {true_mic[gb]:.1f} mg/L")
add_expect = true_mic[gc] * (true_mic[gp]/true_mic[gw])
print(f"  additive expectation       : {add_expect:.1f} mg/L  -> observed {true_mic[gb]:.1f} = super-additive")
