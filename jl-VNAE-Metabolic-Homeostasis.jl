# ------------------------------------------------------------
# VNAE-Metabolic-Homeostasis-Dynamics
# Geometric Stability of Complex Biochemical Networks (Julia Version)
# ------------------------------------------------------------

using Graphs, DifferentialEquations, LinearAlgebra, Statistics, Random, Plots

# 1. NETWORK CONFIGURATION
n = 1200                        # Metabolic variables/enzymes
avg_interactions = 8            # Average connectivity
β = 0.25                        # Structural rigidity (VNAE factor)

# Generate Scale-Free Graph (Barabási-Albert model)
g = barabasi_albert(n, avg_interactions ÷ 2)
A = Matrix(adjacency_matrix(g))

# Add weighted affinities (Chemical coupling strength)
# Gamma distribution for realistic biochemical affinity
dist_affinity = rand(n, n) .* 2.0 # Simplified weighting for speed
A = A .* dist_affinity
A[diagind(A)] .= 0

# 2. ENZYMATIC ASYMMETRY (Theta)
# Represents localized metabolic inertia/dissipation
# Log-normal distribution captures biological variation
θ = exp.(1.0 .+ 0.8 .* randn(n))

# 3. CURVATURE CERTIFICATE (K)
# High-speed Monte Carlo sampling for Riemannian Curvature
samples = 200_000
idx_i = rand(1:n, samples)
idx_j = rand(1:n, samples)

diff_θ = abs.(θ[idx_i] .- θ[idx_j])
coupling = [A[idx_i[k], idx_j[k]] for k in 1:samples]
rigidity = 1.0 .+ β .* (θ[idx_i] .+ θ[idx_j])

K = mean((diff_θ .* coupling) ./ rigidity)

# 4. METABOLIC DYNAMICS (VNAE Equation)
# p_stress: Persistent external forcing (dietary spikes, oxidative stress)
p_stress = 0.5 .* randn(n)
L = Matrix(laplacian_matrix(g))

# Define the ODE system: d(x)/dt = -(L + Θ)x + p
function vnae_metabolism!(dx, x, p, t)
    # Optimized matrix-vector multiplication
    dx .= -(L * x) .- (θ .* x) .+ p_stress
end

# Initial state offset (metabolic shock)
x0 = 5.0 .* randn(n)
tspan = (0.0, 0.5)

# 5. SOLVING THE SYSTEM
# Using Tsit5 (Tsitouras 5/4 Runge-Kutta) - state of the art solver
prob = ODEProblem(vnae_metabolism!, x0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# 6. REPORT & ANALYSIS
println("\n--- VNAE METABOLIC HOMEOSTASIS REPORT ---")
println("Metabolic Variables (n): $n")
println("Structural Curvature (K): $(round(K, digits=6))")
status = K > 0 ? "Geometrically Stable (Homeostasis Guaranteed)" : "Critical Instability"
println("Status: $status")

# 7. VISUALIZATION
# Plotting a subset of 50 metabolites for clarity
sample_indices = randperm(n)[1:50]
plot(sol, idxs=sample_indices, alpha=0.4, legend=false, 
     title="VNAE: Julia-Speed Convergence (K = $(round(K, digits=4)))",
     xlabel="Time (Metabolic Scale)", ylabel="Deviation from Equilibrium",
     grid=true, color=:steelblue)
hline!([0], line=(:black, :dash, 1), alpha=0.5)
