# Equilibrium Temperature Estimation

## Project Overview
This Python code implements a thermodynamic model to compute the equilibrium composition and temperature for **constant-volume combustion** of **dodecane (C₁₂H₂₆)** with oxygen. The model uses **Gibbs free energy minimization** combined with **energy conservation** to determine the final equilibrium state.

---

## Problem Statement

- **Fuel:** C₁₂H₂₆ (Dodecane)  
- **Oxidizer:** O₂  
- **Equivalence Ratio:** φ = 1 (Stoichiometric)

### Initial Conditions
- **Temperature:** T₀ = 300 K  
- **Pressure:** P₀ = 10 bar  
- **Process:** Constant-volume adiabatic combustion  

### Equilibrium Species Considered
CO₂, CO, O, O₂, H₂O, OH, H

---

## Methodology

### 1. Two-Step Iterative Approach

#### **A. Gibbs Free Energy Minimization (Inner Loop)**
For a given temperature **T**, minimize the total Gibbs free energy:

\[
G_{\text{mix}} = \sum_i n_i \mu_i(T)
\]

Subject to:
- Carbon, hydrogen, and oxygen elemental balance  
- SLSQP optimization for mole numbers  
- Atom balance matrix ensuring conservation of C, H, O  

#### **B. Energy Balance (Outer Loop)**
Solve temperature \( T_{eq} \) from:

\[
U_{\text{products}}(T_{eq}) = U_{\text{reactants}}(T_0)
\]

Uses root-finding to enforce internal-energy conservation.

---

## Key Results

### **Equilibrium State**
- **Adiabatic Flame Temperature:** 4381.9 K  
- **Final Pressure:** 232.3 bar  

### **Equilibrium Composition**

| Species | Mole Fraction |
|--------|---------------|
| CO₂ | 0.159 |
| CO  | 0.233 |
| O   | 0.033 |
| O₂  | 0.077 |
| H₂O | 0.351 |
| OH  | 0.112 |
| H   | 0.027 |

---

## Verification

### **Energy Balance**
- Reactants Internal Energy: **−3.4791 × 10⁵ J**  
- Products Internal Energy: **−3.4750 × 10⁵ J**  
- Residual: **+4.07 × 10² J**  
- Relative Error: **0.117%**

### **Elemental Balance**

| Element | Reactants | Products | Error |
|---------|-----------|----------|-------|
| Carbon  | 12.000000 | 12.000002 | 1.61 × 10⁻⁶ |
| Hydrogen| 26.000000 | 26.000003 | 3.49 × 10⁻⁶ |
| Oxygen  | 37.000000 | 37.000000 | 0 |

---

## Technical Implementation

### **Core Functions**

#### `nasa_props(sp, T)`
- Computes thermodynamic properties using NASA polynomials  
- Provides enthalpy, entropy, and Gibbs free energy  
- Supports low (<1000 K) and high (≥1000 K) temperature ranges  

#### `internal_energy_mixture(n, T)`
- Computes mixture internal energy  
- Converts enthalpy into internal energy  

#### `energy_residual(T)`
- Computes internal-energy mismatch  
- Used in root-finding for equilibrium temperature  

---

## Thermodynamic Data
- NASA polynomials for 9 species  
- Includes formation enthalpies, molecular weights, and coefficients for two temperature ranges  

---

## Validation Against Cantera

| Species | This Solver | Cantera | Error (%) |
|---------|-------------|---------|-----------|
| CO₂ | 0.159 | 0.161 | 1.24 |
| CO  | 0.233 | 0.222 | 4.95 |
| O   | 0.033 | 0.039 | 15.38 |
| O₂  | 0.077 | 0.092 | 16.30 |
| H₂O | 0.351 | 0.298 | 17.79 |
| OH  | 0.112 | 0.109 | 2.75 |
| H   | 0.027 | 0.028 | 3.57 |

- **Temperature:** 4381.9 K (solver) vs. 4141.68 K (Cantera) → **5.80% error**  
- **Pressure:** 232.3 bar (solver) vs. 229.872 bar (Cantera) → **1.06% error**

---
