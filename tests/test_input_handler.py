# -------------------------------------------
# OPTIMISATION DEMO  ‑‑  min‑cost CPU/GPU allocation
#
# We rebuild the Delta regression on GPU seconds (same
#   raw data as earlier) and then formulate a tiny MILP
#   with PuLP: choose thread count k_i∈{1,…,8} for each
#   protein so that wall‑clock ≤ D and total Service Units
#   (SU) are minimised.
#
# Assumptions
#   • one A100 GPU per job  (M_i = 1, cannot be split)
#   • GPU SU rate  c_gpu = 16 SU / GPU‑h   (NCSA Delta)
#   • CPU SU rate  c_cpu = 1  SU / core‑h
#   • queue wall‑clock limit D = 2 h = 7200 s
# -------------------------------------------
import pandas as pd, numpy as np, statsmodels.api as sm, pulp, itertools, textwrap

print("Starting optimisation demo...")

# 1 ─── raw execution‑time records for Delta
records = [
    # monomer non‑random
    (97,  1999,  661),
    (255, 3005,  867),
    (493, 2252,  755),
    (985, 3815, 1635),
    # monomer random
    (97,  2550, 1484),
    (255, 2263,  967),
    (493, 2499, 1251),
    (985, 2191, 3752),
    # multimer (use total residues)
    (   97+97,   2692,   687),
    (   97+255,  9495,  1122),
    (   97+493,  7972,  1917),
    (   97+985,  7110,  3565),
    ( 255+255,   3821,  1665),
    ( 255+493,   6431,  2534),
    ( 255+985,  11564,  5307),
    ( 493+493,   2255,  3558),
    ( 493+985,   7375,  8822),
    ( 985+985,   3163, 13836),
]
print("Raw records loaded:", records)
df = pd.DataFrame(records, columns=["length","cpu_s","gpu_s"])
print("DataFrame created:\n", df)

# 2 ─── linear regression  GPU_s = a + b * length
X = sm.add_constant(df["length"])
print("Regression input X:\n", X)
gpu_model = sm.OLS(df["gpu_s"], X).fit()
print("Regression model fitted.")
a_g, b_g = gpu_model.params.values
print(f"Regression coefficients: a_g={a_g}, b_g={b_g}")
print(f"Delta A100 GPU seconds ≈ {a_g:,.1f}  +  {b_g:,.3f} × length")

# 3 ─── optimisation data
lengths = [255, 493, 985]                  # example batch
print("Lengths for optimisation:", lengths)
D = 7200                                   # queue limit (s)
print("Queue wall-clock limit D:", D)
c_cpu, c_gpu = 1, 16                       # SU/ h
print("CPU SU rate:", c_cpu, "GPU SU rate:", c_gpu)
k_choices = range(1,9)                     # 1…8 CPU threads
print("Thread choices:", list(k_choices))
a_c, b_c = 0, 0                            # ignore CPU regression; use table
cpu_lookup = { 255:3005, 493:2252, 985:3815 }  # use non‑random table
print("CPU lookup table:", cpu_lookup)
# fallback: linear fit if length not listed
def cpu_seconds(L):
    result = cpu_lookup.get(L, np.interp(L, df["length"], df["cpu_s"]))
    print(f"cpu_seconds({L}) = {result}")
    return result

# 4 ─── build MILP
print("Building MILP problem...")
prob = pulp.LpProblem("MinCostResourceAssignment", pulp.LpMinimize)
print("MILP problem created.")

# decision vars x[i,k] ∈ {0,1}
x = {(i,k): pulp.LpVariable(f"x_{i}_{k}", cat="Binary")
     for i in range(len(lengths)) for k in k_choices}
print("Decision variables created:", x)

# objective  Σ SU
def su_cost(i,k):
    L = lengths[i]
    cpu_sec = cpu_seconds(L)/k
    gpu_sec = a_g + b_g*L
    cost = (c_cpu*k*cpu_sec + c_gpu*gpu_sec) / 3600  # SU
    print(f"su_cost(i={i}, k={k}) = {cost}")
    return cost
prob += pulp.lpSum( su_cost(i,k)*x[i,k] for i,k in x )
print("Objective function added.")

# constraints  ∑_k x_{i,k}=1  and wall‑clock
for i,L in enumerate(lengths):
    print(f"Adding constraints for i={i}, L={L}")
    prob += pulp.lpSum(x[i,k] for k in k_choices) == 1
    gpu_sec = a_g + b_g*L
    print(f"gpu_sec for L={L}: {gpu_sec}")
    for k in k_choices:
        wall = cpu_seconds(L)/k + gpu_sec
        print(f"Wall time for i={i}, k={k}: {wall}")
        prob += wall * x[i,k] <= D    # enforced only on chosen k

print("Solving MILP...")
prob.solve(pulp.PULP_CBC_CMD(msg=False))
print("Status:", pulp.LpStatus[prob.status])
for i,L in enumerate(lengths):
    k_star = next(k for k in k_choices if pulp.value(x[i,k])>0.5)
    print(f"Selected k for i={i}, L={L}: k_star={k_star}")
    wall_time = cpu_seconds(L)/k_star + a_g + b_g*L
    su = su_cost(i,k_star)
    print(f"  L={L:4d} aa  →  k={k_star} threads ;  wall={wall_time:6.0f}s ; SU={su:5.2f}")