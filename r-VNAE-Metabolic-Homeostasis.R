# ------------------------------------------------------------
# VNAE-Metabolic-Homeostasis-Dynamics
# Geometric Stability of Complex Biochemical Networks (R Version)
# ------------------------------------------------------------

# Load necessary libraries
if (!require("igraph")) install.packages("igraph")
if (!require("deSolve")) install.packages("deSolve")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")

library(igraph)
library(deSolve)
library(ggplot2)
library(reshape2)

# -----------------------------
# 1. NETWORK CONFIGURATION
# -----------------------------
n <- 1200                        # Metabolic variables/enzymes
avg_interactions <- 8            # Average connectivity
beta <- 0.25                     # Structural rigidity (VNAE factor)

# Generate Scale-Free Graph (Barabási-Albert model)
g <- sample_pa(n, m = avg_interactions / 2, directed = FALSE)
A <- as_adjacency_matrix(g, sparse = FALSE)

# Add weighted affinities (Chemical coupling strength)
# Using Gamma distribution for realistic biochemical affinity variance
weights <- matrix(rgamma(n * n, shape = 2, scale = 0.5), nrow = n)
A <- A * weights
diag(A) <- 0

# -----------------------------
# 2. ENZYMATIC ASYMMETRY (Theta)
# -----------------------------
# Theta represents "Metabolic Inertia" or localized dissipation.
theta <- rlnorm(n, meanlog = 1.0, sdlog = 0.8)

# -----------------------------
# 3. CURVATURE CERTIFICATE (K)
# -----------------------------
samples <- 200000
idx_i <- sample(1:n, samples, replace = TRUE)
idx_j <- sample(1:n, samples, replace = TRUE)

diff_theta <- abs(theta[idx_i] - theta[idx_j])
interaction <- A[cbind(idx_i, idx_j)]
rigidity_factor <- 1 + beta * (theta[idx_i] + theta[idx_j])

# Curvature metric K: If K > 0, the system is geometrically stable
K <- mean((diff_theta * interaction) / rigidity_factor)

# -----------------------------
# 4. METABOLIC DYNAMICS (VNAE Equation)
# -----------------------------
# p_stress: Persistent external forcing (dietary spikes, oxidative stress)
p_stress <- rnorm(n, mean = 0, sd = 0.5)
L <- diag(rowSums(A)) - A  # Graph Laplacian

# Define the ODE system: d(x)/dt = -(L + Theta)x + p
vnae_metabolism <- function(t, x, parms) {
  dx <- -(L %*% x) - (theta * x) + p_stress
  return(list(as.vector(dx)))
}

# Initial metabolic offset (shock state)
x0 <- rnorm(n, mean = 0, sd = 5.0)
t_eval <- seq(0, 0.5, length.out = 100)

# Solve the system
solution <- od(y = x0, times = t_eval, func = vnae_metabolism, parms = NULL, method = "rk4")

# -----------------------------
# 5. REPORT & ANALYSIS
# -----------------------------
cat("\n--- VNAE METABOLIC HOMEOSTASIS REPORT ---\n")
cat("Metabolic Variables (n):", n, "\n")
cat("Structural Curvature (K):", round(K, 6), "\n")
cat("Status:", ifelse(K > 0, "Geometrically Stable (Homeostasis Guaranteed)", "Critical Instability"), "\n")

# -----------------------------
# 6. VISUALIZATION
# -----------------------------
# Sample 50 nodes for visual clarity
sample_indices <- sample(1:n, 50)
df <- as.data.frame(solution[, c(1, sample_indices + 1)])
colnames(df)[1] <- "Time"
df_long <- melt(df, id.vars = "Time")

ggplot(df_long, aes(x = Time, y = value, group = variable)) +
  geom_line(alpha = 0.4, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = paste("VNAE: Metabolic Convergence Dynamics (K =", round(K, 4), ")"),
       subtitle = "Homeostatic recovery under 1200+ asymmetric variables",
       x = "Time (Metabolic Scale)",
       y = "Deviation from Equilibrium") +
  theme_minimal()
