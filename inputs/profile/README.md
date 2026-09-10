# Profile likelihood points

One file here is one point of a profile scan. It says nothing about which
parameter is being scanned other than by naming it, so the same three-parameter
set of files works unchanged if `configParSet.yaml` is reordered or if
systematics sets are added or removed.

Run them as SLURM arrays with `profile_dcp.sh`, `profile_dm32.sh` and
`profile_sin23.sh` from the gundam directory, or one at a time:

    gundamFitter -c config.yaml -of ./inputs/profile/dcp_m0900.yaml \
      -O /fitterEngineConfig/generateSamplePlots=false --out-dir outputs/profile/

## What a point file does

Two blocks, both merged into the already-unfolded config by `-of`.

`isFrozen: true` on the scanned parameter. Minuit still declares it and then
calls `FixVariable` on it at the top of every fit
(`RootMinimizer.cpp:242-247`), so it holds still while the other five
oscillation parameters and all 70 systematics minimise around it. `isFixed`
would be wrong here: `MinimizerBase.cpp:62` drops fixed parameters from
`_minimizerParameterPtrList_` before Minuit ever sees them, so there would be no
postfit row for the parameter and no way to release it again from a config
override.

`parameterInjection` puts it at the scan point. This is the propagator's own
config key (`Propagator.cpp:70`), the same thing `--inject-parameters` sets from
the command line, so a point file needs no second file and no extra flag. It is
applied to the model only, after the Asimov data has been built
(`LikelihoodInterface.cpp:131-138`), and `moveParametersToPrior` skips frozen
parameters (`ParameterSet.cpp:476`), so the injected value is still there when
the fit starts.

## Why by name and not by index

`-of` files are merged by `applyOverrides`, which matches list entries on their
`name` field (`GenericToolbox.Json.h:472-511`) and adds keys that are not
already present (`:533-540`). So `isFrozen` does not have to exist in
`configParSet.yaml` beforehand, and nothing here depends on the order of
`parameterDefinitions`.

The command-line form, `-O /..../parameterSetList/0/parameterDefinitions/5/isFrozen=true`,
is shorter but indexes positionally and requires the key to exist already
(`ConfigUtils.cpp:192-216`). Reorder `configParSet.yaml` and it silently freezes
a different parameter. That is why these are files rather than command lines.

## Which parameters, and why only three

`PMNS_SIN_SQUARED_23`, `PMNS_DELTA_MASS_SQUARED_32` and `PMNS_DELTA_CP` - the
three with `priorType: Flat`, which are the three HK measures. Freezing a flat
parameter away from its prior adds no penalty term, so the curve is the
statistical likelihood plus the penalties the *other* parameters pay to absorb
the shift, which is what a profile is meant to show.

The other three (`PMNS_SIN_SQUARED_12`, `PMNS_SIN_SQUARED_13`,
`PMNS_DELTA_MASS_SQUARED_21`) are Gaussian-constrained from reactor and solar
measurements. Profiling them would work with exactly the same file shape - copy
one and change the two names and the value - but the curve would mostly retrace
the external constraint that was put in, since nothing in these five samples
measures them.

`PMNS_SIGN_MASS_SQUARED_32` is already `isFixed: true` and has two values, not a
range.

## The grids

148 points, spanning the bounds already declared in `configParSet.yaml`. The
last column is the step in units of that parameter's 1 sigma from
`MakeCovariances.C`, which is what says whether the grid can resolve the curve:

| parameter | bound in config | grid | prior | sigma | step |
| --- | --- | --- | --- | --- | --- |
| `PMNS_SIN_SQUARED_23` | `parameterLimits: [0.3, 0.7]` | 0.300 to 0.692 by 0.008 | 0.528 | 0.021 | 0.38 sigma |
| `PMNS_DELTA_MASS_SQUARED_32` | `parameterLimits: [2.0E-3, 3.0E-3]` | 2.000E-3 to 2.980E-3 by 2.0E-5 | 2.509E-3 | 3.4E-5 | 0.59 sigma |
| `PMNS_DELTA_CP` | `cyclicRange: [-pi, pi]` | -180 to +172.5 deg by 7.5 deg, 48 points | -1.601 (-91.7 deg) | 1.17 (67 deg) | 0.11 sigma |

Steps are round in the parameter, which is why the top of the first two grids
stops short of the upper limit - 0.692 rather than 0.700, and 2.980E-3 rather
than 3.000E-3. Nothing lives up there; both ends are already about 5 sigma out.

delta-CP is the exception twice over. It gets the full range, because it is
cyclic and its interval is quoted over all of it rather than as a symmetric
error. And it gets 48 points rather than 50: a 7.5 degree step puts 0 and both
maximal-CP-violation points on the grid, where 360/50 = 7.2 degrees would miss
+-90. Forty-eight steps of 7.5 close 360 exactly, so +180 is not a separate
point - it is `dcp_m1800`.

File names carry the value in a form that survives GUNDAM's output-name
cleaning: sin23 in thousandths, dm32 in units of 1E-6 eV^2, delta-CP in tenths
of a degree with `m`/`p` for the sign. `generateCleanFileName`
(`GundamUtils.cpp:79`) turns a `-` or a `.` into `_`, which would make
`dcp_-0.785` and `dcp_+0.785` collide.

No grid is centred on the prior, so the vertex of each curve falls between two
sampled points. That costs nothing: the Asimov minimum is the prior by
construction, with chi2 = 0, and `asimov.sh` already reports it.

## Reading the output

One ROOT file per point in `outputs/profile/`, named from the override file, e.g.
`gundamFitter_config_With_dcp_m0900.root` - the `With_` comes from the
`overrideFiles` appendix (`gundamFitter.cxx:214`). The number to plot against
the scan value is the postfit total likelihood. `compare3way.C` reads these
next to the 1D and 3D grid outputs.
