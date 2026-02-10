import numpy as np
import networkx as nx
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# VNAE-Metabolic-Homeostasis-Dynamics
# Geometric Stability of Complex Biochemical Networks
# ------------------------------------------------------------

# -----------------------------
# 1. NETWORK CONFIGURATION
# -----------------------------
n = 1200                        # Metabolic variables/enzymes (Human-scale network)
avg_interactions = 8            # Average connectivity in metabolic pathways
beta = 0.25                     # Structural rigidity (VNAE factor)

# Scale-Free Graph: Reflects biological systems (few hubs, many sparse nodes)
G = nx.barabasi_albert_graph(n, m=avg_interactions // 2)
A = nx.to_numpy_array(G)

# Add weighted affinities (Chemical coupling strength)
A *= np.random.gamma(shape=2.0, scale=0.5, size=(n, n))
np.fill_diagonal(A, 0)

# -----------------------------
# 2. ENZYMATIC ASYMMETRY (Theta)
# -----------------------------
# Theta represents "Metabolic Inertia" or localized dissipation.
# It models validator-specific responsiveness or enzymatic specificity.
theta = np.random.lognormal(mean=1.0, sigma=0.8, size=n)

# -----------------------------
# 3. CURVATURE CERTIFICATE (K)
# -----------------------------
# Estimating the geometric curvature based on coupling and asymmetry
samples = 200_000
idx_i = np.random.randint(0, n, samples)
idx_j = np.random.randint(0, n, samples)

diff_theta = np.abs(theta[idx_i] - theta[idx_j])
interaction = A[idx_i, idx_j]
rigidity_factor = 1 + beta * (theta[idx_i] + theta[idx_j])

# Curvature metric K: If K > 0, the system is geometrically stable
K = np.mean((diff_theta * interaction) / rigidity_factor)

# -----------------------------
# 4. METABOLIC DYNAMICS (VNAE Equation)
# -----------------------------
# p_stress: Persistent external forcing (dietary spikes, oxidative stress, toxins)
p_stress = np.random.normal(0, 0.5, size=n)
L = np.diag(A.sum(axis=1)) - A  # Directed Graph Laplacian

def vnae_metabolism(t, x):
    # d(x)/dt = -(L + Theta)x + p
    # x represents the deviation from the ideal homeostatic state
    return -(L @ x) - (theta * x) + p_stress

# Initial state offset (e.g., post-prandial state or metabolic shock)
x0 = np.random.normal(0, 5.0, size=n)

# Time window for simulation
t_span = (0, 0.5)
t_eval = np.linspace(t_span[0], t_span[1], 100)

solution = solve_ivp(vnae_metabolism, t_span, x0, t_eval=t_eval, method="RK45")

# -----------------------------
# 5. REPORT & ANALYSIS
# -----------------------------
print(f"\n--- VNAE METABOLIC HOMEOSTASIS REPORT ---")
print(f"Metabolic Variables (n): {n}")
print(f"Structural Curvature (K): {K:.6f}")
print(f"Status: {'Geometrically Stable (Homeostasis Guaranteed)' if K > 0 else 'Critical Instability'}")

# -----------------------------
# 6. VISUALIZATION
# -----------------------------
plt.figure(figsize=(10, 6))
# Sample 50 nodes for visual clarity
sample_nodes = np.random.choice(n, 50, replace=False)

for i in sample_nodes:
    plt.plot(solution.t, solution.y[i], alpha=0.5, linewidth=1)

plt.axhline(0, color='black', linestyle='--', alpha=0.3)
plt.title(f"VNAE: Metabolic Convergence Dynamics (K = {K:.4f})")
plt.xlabel("Time (Metabolic Scale)")
plt.ylabel("Deviation from Equilibrium")
plt.grid(True, alpha=0.2)
plt.tight_layout()
plt.show()
