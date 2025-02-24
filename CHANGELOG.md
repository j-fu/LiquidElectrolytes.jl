# Changelog

All notable changes to this project will be documented in this file.

## 0.3.0

### Features
- Add rotating disk code from RotatingElectrodesProject
- ExtendableFEM extension & coupling with stokes
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
- Introduce `exp` and `log` (`Base.exp`, `Base.log` by default= as part of electrolyte data, remove epsreg
  - ElectrolyteData is now parametrized with the types of these functions
  - This allows for a more transparent handling of log and exp regularization
  - Provide `RExp` and `RLog` functor structs providing regularized alternative to exp and log
  - Minimize use of regularization in examples
  - Remove rlog and rexp. Projects can define instead `rlog=RLog()` and `rexp=RExp()`

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
