# Changelog

All notable changes to this project will be documented in this file.

## 2.0.2 - 2025-05-18 
- Documentation update

## 2.0.1 - 2025-05-14 
- Bugfixes + additions wrt. 2.0.0. 2.0.0 wasn't fully working

## 2.0.0 - 2025-05-14

### Breaking
All code which just used default configurations should continue to work.

- `scheme` paramter of ElectrolyteData replaced by `upwindflux!` parameter which should be one of `act_flux!`, `μex_flux!` and `cent_flux!`. The default is `μex_flux!`, based on the chemical excess potential, formerly selected via `schem=:μex`. This has been shown to be superior to the alternatives (`act_flux!` / `:act`,  `cent_flux!` / `:cent`). Unless experiments comparing these are required, there is no need to keep carrying the selection of schemes in the examples. Code setting `ElectrolyteData.scheme` still works, but the setting is ignored and a warning is emitted.
- `aflux`, `sflux` and `cflux` calculating just one species upwind flux replaced with  `act_flux!`, `μex_flux!` and `cflux!`, respectively, which  calculate all ionic fluxes at once.
- PNPSystem and PBSystem are now distinct types in order to be able to implement more specific methods
- `pnpunknowns` has been removed in favor of a method of `VoronoiFVM.unknowns`
- `rexp`, `rlog` are now part of `ElectrolyteData`

### Features
- ElectrolyteData now has a field  `actcoeff!` for a user definable function which calculates activity coefficients
- Fixed implementation of PBSystem, now cases with different `κ`  etc. for different species should work
- Pre-calculate more data in ElectrolyteData, use `update_derived!` to update these fields.

## 1.1.0 - 2025-03-16

### Features
- Add cvsweep 
- Reorgnization of source files
- Document which entries of ElectrolyteData are modified by algorithms
- Use ReTest, tighten CI with Aqua etc

### Breaking
- `molarity` kwarg is removed from `dlcapsweep`.
- Removed undocumented EquilibriumData stuff, moved to tests

## 1.0.0 - 2025-02-26

### Features
- Add rotating disk code from RotatingElectrodesProject
- ExtendableFEM extension & coupling with Stokes
- add `chemical_potentials` and `electrochemical_potentials` methods for  calculation from stationary solutions
- Remove bikerman special case in pnpflux
- forward/backward progress for dlcap
- Scale pressure equation by 1/pscale
- Added pre-factor `1/v0` to activity coefficients (this does not influence excess chemcial potential gradients
  but makes things more compatible to literature)
- Checked computation with Dual64

### Breaking
- charge() renamed to chargedensity()
- remove export of iA, iC
- Define `rexp(::Any)` and `rlog(::Any)` as `Base.log` and `Base.exp()`, remove epsreg from ElectrolyteData.
  Use these functions for all `exp` and `log` calculations in the package.
  By defining `rexp(::Number)` and `rlog(::Number)`, users can overwrite these with e.g. regularized  versions.
- Provide `RExp` and `RLog` functor structs providing regularized alternative to exp and log, minimize use of these regularization in examples

### Bugfixes
- Re-checked and fixed all formulas again
- Fix calculation of chemical_potentials! (take into account solvation)
- Use @view in solventoncentration
- Proper initialization of unknowns for PBSystem in examples, we don't need regularization here

## [0.2.4] - 2024-10-02

- CI test on apple silicon
- Require LessUnitful v1
- Allow for VoronoiFVM v2

## [0.2.3] - 2024-02-23

- Use ExampleJuggler v2
- Bugfix in currents(::IVSweepResult): Faraday constant was wrong by chance.

## [0.2.2] - 2024-02-03

- improve methods to obtain currents from return results
- improved visualizations
- Documentation + tests with ExampleJuggler.jl
