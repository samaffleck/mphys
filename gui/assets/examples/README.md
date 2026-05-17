# GUI state examples

Each `.json` file can be loaded via **File → Open** in `mphys_gui`.  
The state is saved with [cereal](https://uscilab.github.io/cereal/) (JSON archive).

---

## darcy_packed_bed.json

1D compressible packed-bed with Darcy's Law + ideal-gas mass balance.

### Equations

| # | Name | Form |
|---|------|------|
| 1 | Darcy's Law (algebraic) | u + (κ/μ) ∂P/∂x = 0 |
| 2 | Mass balance (transient) | ∂P/∂t + ∂(Pu)/∂x = 0 |

### Setup

| Parameter | Value |
|-----------|-------|
| Domain length | 1 m |
| Permeability κ | 1 m² |
| Viscosity μ | 1 Pa·s |
| P (x=0) | 2 Pa |
| P (x=1) | 1 Pa |
| n\_cells | 100 |
| t\_end | 10 s |

κ/μ = 1 and the domain length are chosen for unit-scale validation (not physical values).

### Analytical steady-state

Substituting Darcy into the mass balance gives ∂(P·∂P/∂x)/∂x = 0, so P² is linear in x:

```
P(x) = sqrt(P_L² + (P_R² - P_L²)·x)  =  sqrt(4 - 3x)
u(x) = (κ/μ) · 3 / (2·sqrt(4 - 3x))
```

Mass flux ρu ∝ Pu = (3/2) = const  (conserved, as expected).

Run to t = 10 s; the solution is effectively at steady state by t ≈ 2 s.
