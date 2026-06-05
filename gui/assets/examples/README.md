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

---

## single_particle_model.json

Single Particle Model (PyBaMM SPM): two spherical electrode particles on a
normalised radius r/R ∈ [0,1] with a closed-form terminal voltage.  Load it,
open the **Study** node and press **Run** for a constant-current discharge.

See `include/mphys/models/spm.hpp`.

---

## single_particle_model_electrolyte.json

Single Particle Model with Electrolyte (PyBaMM SPMe).  Adds a macroscopic
electrolyte-concentration field across the negative | separator | positive
sandwich plus the leading-order electrolyte corrections to the terminal
voltage.

### Equations

| # | Name | Form |
|---|------|------|
| 1,2 | Particle diffusion | ∂c_s,k/∂t = D_s,k/r² ∂/∂r(r² ∂c_s,k/∂r),  k∈{n,p} |
| 3 | Electrolyte transport | ε ∂c_e/∂t = ∂/∂x(ε^b D_e ∂c_e/∂x) + S(x) |

with source S = +(1−t⁺)i_app/(F L_n) in the negative electrode, 0 in the
separator, −(1−t⁺)i_app/(F L_p) in the positive electrode, and no-flux
electrolyte boundaries.

### Voltage

```
V = U_p(y) − U_n(x)                              (open-circuit)
  − reaction overpotentials (j0 uses local c_e)
  + 2(1−t⁺)(RT/F)(ln c̄_e,p − ln c̄_e,n)          (concentration overpotential)
  − i_app(L_n/3κ_n + L_s/κ_s + L_p/3κ_p)         (electrolyte ohmic)
  − (i_app/3)(L_n/σ_n + L_p/σ_p)                 (solid ohmic)
```

### Checks

- Negative-particle couliometry: d⟨c_n⟩/dt = −I/(A F ε_n L_n) (exact).
- Electrolyte conservation: porosity-weighted ⟨ε c_e⟩ is constant in time.
- Under constant current the electrolyte relaxes to the analytic SPMe
  quadratic/linear/quadratic steady profile (see `tests/test_spme.cpp`).

The **Results** node plots terminal voltage vs time, the particle
concentrations vs r/R, and the electrolyte concentration vs x.

See `include/mphys/models/spme.hpp`.
