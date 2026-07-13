(omega-dev-frazil)=

# Frazil

This page describes frazil design and implementation details in Omega,
including both the `basic` and `teos` pathways.

## Purpose and coupling points

Frazil computes phase-change-related tendencies that modify:

- pseudo-thickness tendency
- temperature tracer tendency
- salinity tracer tendency

The tendency hook-up is implemented through `FrazilOnCell` in the tracer
tendency phase, where frazil contributions are added to the accumulated
`PseudoThicknessTend` and `TracerTend` arrays.

## Data flow and call sequence

1. `Tendencies::computeTracerTendenciesOnly` checks
   `Tendencies.FrazilTendencyEnable`.
2. If enabled, `FrazilOnCell::operator()` retrieves the default `Frazil`
   object and zeros frazil tendency/accumulator arrays.
3. `FrazilOnCell` extracts `Temperature` and `Salinity` tracer subviews,
   then calls `Frazil::computeFrazil(CT, SA, PressureMid, PseudoThickness)`.
4. `Frazil::computeFrazil` dispatches to `computeFrazilBasicImpl` or
   `computeFrazilTeosImpl` based on `FrazilType`.
5. Returned frazil tendencies are added into `PseudoThicknessTend` and
   `TracerTend` for all active cell layers.

## Configuration coupling

Frazil behavior is configured with:

- `Omega.Tendencies.FrazilTendencyEnable`
  - global switch for applying frazil tendency terms
- `Omega.Frazil.FrazilType`
  - implementation choice (`basic` or `teos`)
- `Omega.Frazil.MassLimit`
  - per-layer mass/thickness limiter used by frazil formation/melt pathways
- `Omega.Frazil.Phi`
  - teos pathway liquid-fraction parameter for new frazil partitioning
- `Omega.Frazil.DepthLimit`
  - optional depth cutoff for frazil activity; negative means no cutoff
- `Omega.Frazil.ConservationCheck`
  - optional post-compute column conservation diagnostic/logging

## Physics/algorithm summary

### Common behavior

- Freezing-point checks are based on conservative temperature and absolute
  salinity with pressure-dependent freezing temperature.
- Vertical accumulation order is bottom-to-top within each active column.
- Frazil tendencies are not time-step scaled inside `Frazil`; they are
  accumulated as tendency contributions.

### Basic pathway (`FrazilType: basic`)

- Uses simplified energetics: the energy of the super-cooled water sets the
amount of solid ice formed (used constant latent heat of fusion of fresh ice).
Salt is added based on a constant bulk salinity `IceRefSal` (default) or a
manual toggle (for now) using the local salinity and the frazil porosity.
Melting of existing frazil is set by the amount of pure ice that can be melted
by the warm layer, and the enthalpy of melt water at the local freezing point
is added to the temperature tendency.
  temperature comparisons and fractional thickness limits.
- Computes local layer tendencies (`HTend`, `TTend`, `STend`) and updates
  accumulated frazil stores.
- Applies surface salt redistribution adjustment and converts accumulators to
  coupler units at the end of the column loop.
  -Warnings: 1) basic frazil formation/melt is does not conserve energy.
  2) Using porosity to set the salt content includes a redistribution of excess
  salt in the surface layer, which can be very significant.

### Teos pathway (`FrazilType: teos`)

- Uses teos-10 Gibbs SeaWater routines for frazil formation and melt state
  transitions.
- Formation uses `Phi` and `MassLimit` to partition and limit newly formed
  frazil contributions.
- Melt computes fraction melted subject to available thermodynamic energy and
  mass-limit constraints.

## Existing ctest coverage

The existing frazil test driver is in
`components/omega/test/ocn/FrazilTest.cpp` and covers:

- teos frazil formation in cold and warm single-layer states
- basic frazil formation in cold and warm single-layer states
- mixed warm/cold column behavior with sign checks for branch switching
- depth-limit behavior ensuring excluded deep layers have zero frazil tendency


## Related pages

- User-facing options: [User Frazil Guide](../userGuide/Frazil.md)
- Tendency container/hook-up: [Tendencies](Tendencies.md)
