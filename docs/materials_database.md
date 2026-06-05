# Materials database

A library of reusable, geometry-independent physical properties for battery
components, so a user can pick (say) **graphite** for the negative electrode and
**NMC811** for the positive electrode and feed them straight into a model such
as the Single Particle Model.

## Design principles

1. **A material owns only intrinsic properties.** The equilibrium potential of
   NMC811, the maximum lithium concentration of graphite, the conductivity of
   1 M LiPF6 in EC:EMC — these are properties of the *material* and are reused
   across every cell that uses it. They live in the database.

2. **Cell design and operating conditions live elsewhere.** Electrode
   thickness, particle radius, electrode area, active-material volume fraction,
   applied current, temperature — these describe a *cell* or an *experiment*,
   not a material. They are supplied separately (`CellDesign`, current, T).

3. **Materials are grouped by `Category` and tagged with a `Domain`.**
   - `Category::kElectrode` — active materials (graphite, NMC811, LFP, NMC532).
   - `Category::kElectrolyte` — salt + solvent systems (LiPF6 in EC:EMC, …).
   - `Domain` (`kNegativeElectrode`, `kPositiveElectrode`, `kElectrolyte`, …) is
     a soft hint used to filter pickers — e.g. only offer graphite-like
     materials for the negative slot.

4. **A model declares the material inputs it needs.** This is the link between
   "material" and "model". Each model exposes a list of `MaterialSlot`s naming
   the roles to fill (negative electrode, positive electrode, electrolyte) and
   the category/domain that fits each. A UI builds its pickers from this list;
   an adapter turns the chosen materials + cell design into the model's
   parameter struct.

## File layout

```
include/mphys/materials/
  material.hpp        Core types: Category, Domain, property function types,
                      ElectrodeMaterial, ElectrolyteMaterial, MaterialSlot.
  electrodes.hpp      Electrode data + OCP fits (graphite, NMC811, LFP, NMC532).
  electrolytes.hpp    Electrolyte data (LiPF6 EC:EMC Nyman 2008, EC:DMC Marquis).
  database.hpp        Database: enumerate / look up materials, filter by domain.
  spm_materials.hpp   CellDesign + SpmMaterialSlots() + MakeSpmParameters():
                      the bridge from a material selection to SpmParameters.
  spme_materials.hpp  SpmeCellDesign + SpmeMaterialSlots() + MakeSpmeParameters():
                      bridge to SpmeParameters; samples the electrolyte's
                      transport functions at the nominal concentration.
```

Everything is header-only, matching the existing `include/mphys/models/`.

## Property representation

| Kind | Type | Example |
|------|------|---------|
| Electrode OCP, entropic coeff. | `ScalarProperty` = `double(double x)` | `U(x)` vs Li/Li+, `x = c/c_max` |
| Electrolyte κ, D_e, t₊, TDF | `StateProperty` = `double(double c, double T)` | `kappa(c, T)` [S/m] |
| Concentrations, diffusivities, rates | `double` | `c_max`, `D_s`, `k` |

Solid-phase diffusivity and the kinetic rate constant are stored as scalars
today (the SPM uses constants); each material also carries an Arrhenius
activation energy field (default 0 = no temperature dependence) so a future
model can apply a correction without changing the schema. The electrolyte
properties are already concentration/temperature functions.

## Usage

Built-in materials are identified by **strongly-typed enums** (`ElectrodeId`,
`ElectrolyteId`), not strings — editor autocomplete lists every available
material and a typo is a compile error, not a runtime throw. Each material still
carries a human-readable `name` for UI labels.

```cpp
#include "mphys/materials/database.hpp"
#include "mphys/materials/spm_materials.hpp"

using namespace mphys::materials;

CellDesign cell;  // geometry + initial state
auto p = MakeSpmParameters(ElectrodeId::kGraphiteChen2020,   // negative
                           ElectrodeId::kNmc811Chen2020,     // positive
                           ElectrolyteId::kLipf6EcEmcNyman2008,
                           cell, /*current_A=*/5.0, /*T=*/298.15);
// ... build the SpmModel with p as today.

// Or grab a material directly:
const auto& nmc = Database::Electrode(ElectrodeId::kNmc811Chen2020);
double u = nmc.ocp(0.5);   // OCP at 50% lithiation
```

For a UI, enumerate ids (optionally filtered by domain) and label each with its
display name:

```cpp
for (ElectrodeId id : Database::ElectrodeIds(Domain::kPositiveElectrode)) {
  const std::string& label = Database::ElectrodeName(id);
  // add a selectable for `label`, storing `id`
}
```

## Built-in materials

| Material | Category | Domain | Source |
|----------|----------|--------|--------|
| Graphite (Chen 2020) | electrode | negative | Chen et al. 2020 (LG M50) |
| NMC811 (Chen 2020) | electrode | positive | Chen et al. 2020 (LG M50) |
| LFP (Prada 2013) | electrode | positive | Prada et al. 2013 |
| NMC532 (Ecker 2015) | electrode | positive | Ecker et al. 2015 (representative) |
| LiPF6 in EC:EMC (Nyman 2008) | electrolyte | electrolyte | Nyman et al. 2008 |
| LiPF6 in EC:DMC (Marquis 2019) | electrolyte | electrolyte | Marquis et al. 2019 |

## Extending

- **New electrode/electrolyte:** add a builder in `electrodes.hpp` /
  `electrolytes.hpp`, add an enumerator to `ElectrodeId` / `ElectrolyteId`, and
  wire it into the `Database::Electrode(id)` switch and the `ElectrodeIds()` list.
- **New model:** add its own `…MaterialSlots()` and a `Make…Parameters()`
  adapter alongside the model, following `spm_materials.hpp` /
  `spme_materials.hpp`. The SPMe adapter shows how to consume the electrolyte's
  `conductivity` / `diffusivity` / `transference` functions (sampled at the
  nominal concentration for the leading-order constant-property model).
- **Temperature dependence:** populate `activation_energy_D` /
  `activation_energy_k` and apply the Arrhenius factor in the adapter.

## GUI integration

The SPM and SPMe parameter editors (`gui/mphys_gui.cpp`) start with a
**Materials** section: dropdowns for the negative electrode, positive electrode
and electrolyte, populated from `Database::Electrode/ElectrolyteNames()` (the
electrode lists filtered by domain). Choosing a material loads its scalar
properties into the editable fields below and selects its OCP curve — the OCP is
sourced from the database in `ToSpmParameters()` because it is a curve, not an
editable scalar (this is what makes a non-graphite/NMC choice such as LFP behave
correctly). The selection round-trips through the saved-state JSON.

## Status / next steps

- [x] Core types, electrode + electrolyte data, registry, SPM adapter, tests.
- [x] SPMe adapter (`spme_materials.hpp`) consuming electrolyte transport.
- [x] Material-selection pickers wired into the SPM + SPMe GUI editors.
- [ ] Optional: load materials from JSON (Cereal is already a dependency) so
      users can add materials without recompiling.
- [ ] Optional: GUI preview of the selected OCP curve.
