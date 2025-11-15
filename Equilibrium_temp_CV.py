import numpy as np
from scipy.optimize import minimize, root_scalar

R = 8.314462618  # J/mol-K
P_in = 10 * 101325  # Pa (10 bar)
T_in= 300
p0 = 101325      # reference pressure Pa
V = ((1 + 18.5) * R * T_in)/P_in



thermo = {

    
    'C12H26': {
        'low': [1.677e+00, 5.499e-02, -2.838e-05, 6.353e-09, -5.298e-13, -3.885e+04, 2.308e+01],
        'high': [3.850e+00, 1.213e-01, -6.236e-05, 1.408e-08, -1.222e-12, -4.414e+04, 1.418e+00],
        'T_switch': 1000.0,
        'Mw': 170.33,
        'Hf298': -291700
    },
    'O2': {
        'low': [3.78245636e+00, -2.99673416e-03, 9.84730200e-06, -9.68129500e-09, 3.24372800e-12, -1.06394356e+03, 3.65767573e+00],
        'high': [3.28253784e+00, 1.48308754e-03, -7.57966669e-07, 2.09470555e-10, -2.16717794e-14, -1.08845772e+03, 5.45323129e+00],
        'T_switch': 1000.0,
        'Mw': 32.00,
        'Hf298': 0.0
    },
    'CO2': {
        'low': [3.85746029e+00, 4.41437026e-03, -2.21481404e-06, 5.23490188e-10, -4.72084164e-14, -4.87591660e+04, 2.27163806e+00],
        'high': [4.63659493e+00, 2.74131991e-03, -9.95828531e-07, 1.60373011e-10, -9.16103468e-15, -4.90249341e+04, -1.93534855e+00],
        'T_switch': 1000.0,
        'Mw': 44.01,
        'Hf298': -393522
    },
    'CO': {
        'low': [3.57953347e+00, -6.10353680e-04, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -1.43440860e+04, 3.50840928e+00],
        'high': [3.04848583e+00, 1.35172818e-03, -4.85794075e-07, 7.88536486e-11, -4.69807489e-15, -1.42661171e+04, 6.01709790e+00],
        'T_switch': 1000.0,
        'Mw': 28.01,
        'Hf298': -110527
    },
    'O': {
        'low': [3.16826710e+00, -3.27931884e-03, 6.64306396e-06, -6.12806624e-09, 2.11265971e-12, 2.91222592e+04, 2.05193346e+00],
        'high': [2.56942078e+00, -8.59741137e-05, 4.19484589e-08, -1.00177799e-11, 1.22833691e-15, 2.92175791e+04, 4.78433864e+00],
        'T_switch': 1000.0,
        'Mw': 16.00,
        'Hf298': 249175
    },
    'H2O': {
        'low': [4.19864056e+00, -2.03643410e-03, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -3.02937267e+04, -8.49032208e-01],
        'high': [3.03399249e+00, 2.17691804e-03, -1.64072518e-07, -9.70419870e-11, 1.68200992e-14, -3.00042971e+04, 4.96677010e+00],
        'T_switch': 1000.0,
        'Mw': 18.02,
        'Hf298': -241826
    },
    'OH': {
        'low': [3.99201543e+00, -2.40131752e-03, 4.61793841e-06, -3.88113333e-09, 1.36411470e-12, 3.61508056e+03, -1.03925458e-01],
        'high': [3.09288767e+00, 5.48429716e-04, 1.26505228e-07, -8.79461556e-11, 1.17412376e-14, 3.85865700e+03, 4.47669610e+00],
        'T_switch': 1000.0,
        'Mw': 17.01,
        'Hf298': 38987
    },
    'H': {
        'low': [2.50000000e+00, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18, -9.27732332e-22, 2.54736599e+04, -4.46682853e-01],
        'high': [2.50000001e+00, -2.30842973e-11, 1.61561948e-14, -4.73515235e-18, 4.98197357e-22, 2.54736599e+04, -4.46682914e-01],
        'T_switch': 1000.0,
        'Mw': 1.01,
        'Hf298': 218000
    },
    'N2': {
        'low': [3.53100528e+00, -1.23660988e-04, -5.02999437e-07, 2.43530612e-09, -1.40881235e-12, -1.04697628e+03, 2.96747468e+00],
        'high': [2.95257626e+00, 1.39690057e-03, -4.92631691e-07, 7.86010367e-11, -4.60755321e-15, -9.23948645e+02, 5.87189252e+00],
        'T_switch': 1000.0,
        'Mw': 28.01,
        'Hf298': 0.0
    }
} 


def nasa_props(sp, T):
    """Return h(T), s(T), g(T) for species [J/mol]."""
    data = thermo[sp]
    if T < 1000:
        a = data["low"]
    else:
        a = data["high"]

    
    # Heat capacity
    cp_R = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4
    h_RT = a[0] + a[1]*T/2 + a[2]*T**2/3 + a[3]*T**3/4 + a[4]*T**4/5 + a[5]/T
    s_R = a[0]*np.log(T) + a[1]*T + a[2]*T**2/2 + a[3]*T**3/3 + a[4]*T**4/4 + a[6]

    # Convert to absolute values
    cp = cp_R * R  # J/mol·K
    h = h_RT * R * T  # J/mol
    s = s_R * R  # J/mol·K
    g = h - T*s
    

    return h, s, g



n1 = 1
n2 = 18.5
# Atom balance matrix A (rows = [C,H,O], cols = species)
A = np.array([
    [12 , 1, 1, 0, 0, 0, 0, 0],   # Carbon
    [26 , 0, 0, 0, 0, 2, 1, 1],   # Hydrogen
    [0 , 2, 1, 1, 2, 1, 1, 0]    # Oxygen
], dtype=float)

# Initial atoms from reactants
b = np.array([12*n1, 26*n1, 2*n2], dtype=float)

# Placeholder species list
species = ["C12H26","CO2","CO","O","O2","H2O","OH","H"]  # Example list


# ----------------------------
# Gibbs free energy function
# ----------------------------
def G_of_n(n, T, V=0.048000888):
    Ntot = np.sum(n)
    y = n / Ntot
    g0 = np.array([nasa_props(sp,T)[2] for sp in species])
    pi = y * Ntot * R * T / V   # partial pressure from EOS
    return np.sum(n*(g0 + R*T*np.log(np.maximum(pi/p0,1e-30))))

def element_balance(n):
    return A.dot(n) - b

# ----------------------------
# Inner minimizer: composition at T
# ----------------------------
def equilibrium_at_T(T):
    n0 = np.ones(len(species))  # initial guess
    cons = [{'type':'eq','fun':lambda n,j=j: element_balance(n)[j]} for j in range(3)]
    bounds = [(1e-20,None)]*len(species)
    res = minimize(lambda n: G_of_n(n,T), n0, method='SLSQP',
                   constraints=cons, bounds=bounds)
    return res.x


# ----------------------------
# Internal energy mixture
# ----------------------------
def internal_energy_mixture(n, T):
    """Return mixture internal energy [J]"""
    return np.sum([ni*(nasa_props(sp,T)[0] - R*T) for ni,sp in zip(n,species)])



# Reactant internal energy (fuel + oxidizer at 300 K)
n_react = np.zeros(len(species))
n_react[species.index("C12H26")] = 1.0
n_react[species.index("O2")] = 18.5
Uref = internal_energy_mixture(n_react, 300.0)

def energy_residual(T):
    n_eq = equilibrium_at_T(T)
    Uprod = internal_energy_mixture(n_eq, T)
    return Uprod - Uref

# ----------------------------
# Solve for constant-volume equilibrium T
# ----------------------------
sol = root_scalar(energy_residual, bracket=[1500,5000], method='bisect')
T_eq = sol.root
n_final = equilibrium_at_T(T_eq)
X_final = n_final/np.sum(n_final)

# ----------------------------
# Print results
# ----------------------------

P = (sum(n_final)* R *T_eq)/V

print(f"Equilibrium Temperature (const-V) = {T_eq:.1f} K at {P/1e5:.1f} bar\n")
print("Equilibrium mole fractions:")
for sp, xi in zip(species, X_final):
    if xi > 1e-6:
        print(f"{sp:6s} : {xi:.6e}")

# ----------------------------
# Energy balance check
# ----------------------------
Uprod_final = internal_energy_mixture(n_final, T_eq)
residual = Uprod_final - Uref
rel_error = residual / abs(Uref) if Uref != 0 else np.nan

print("\n--- Energy Balance Check (const-V) ---")
print(f"Uref  (reactants at 300 K) = {Uref:.6e} J")
print(f"Uprod (products at Teq)    = {Uprod_final:.6e} J")
print(f"Residual (Uprod - Uref)    = {residual:.6e} J")
print(f"Relative error             = {rel_error:.3e}")

# ----------------------------
# Elemental balance check
# ----------------------------
elem_react = A.dot(n_react)
elem_prod  = A.dot(n_final)

print("\n--- Elemental Balance Check ---")
elements = ["C", "H", "O"]
for e, r, p in zip(elements, elem_react, elem_prod):
    print(f"{e:2s}: Reactants = {r:.6f}, Products = {p:.6f}, Difference = {p-r:.3e}")

