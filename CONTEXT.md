# myFuelCell: current state (context for a new chat)

Written 2026-11-12 from the code at `e1ad218`; updated 2026-11-18 for `master` = `0cc4a1e` (species source, zone-normalised rescale, rate constant from Y, in-run stop rule). It replaces the untracked `report.md` in this repo, which was written at `9b28b91` and is now out of date (it still says alpha = 0.5 only, and it still lists the stale volumetric flux as an open item).

**`master` is the main version for all simulations (fit and sweep), locally and on ARIS.**

## 0. Quick facts

| | |
|---|---|
| Repo | `~/OpenFOAM/cfdlab-v2312/applications/myFuelCell` (local git, no remote) |
| What it is | a single-phase, pruned copy of openFuelCell2 (branch `coupled-phi-parallel` @ `7a4480b`), about 19–21k lines instead of 83.8k, plus a restructured electric solve |
| OpenFOAM | **v2312**. The login shell loads v2406, so always run `source /usr/lib/openfoam/openfoam2312/etc/bashrc` first. `WM_PROJECT_VERSION` is `v2312` (not `2312`). |
| Binaries | app `myFuelCell`; libs `libmyFuelCellSystem`, `libmyRegionCourantNo`, `libmyFuelCellTurbulenceModels`; tool `regionProcAddressing`. All live in `~/OpenFOAM/cfdlab-v2312/platforms/linux64GccDPInt32Opt`, next to the openFuelCell2 libs, which is why the names differ. |
| Not built here | `decomposeParID`, `topoSet`, `libmasterRegionToCell`. These are still taken from openFuelCell2. |
| Build | `src/Allwmake` (libs, app, tool). If you add a new header, run `wmakeLnInclude -u src/libSrc/fuelCellSystems` first. |
| Reference case | `~/OpenFOAM/cfdlab-v2312/run/H2/HTPEMFC/warmStart/base`: HT-PEMFC straight channel, 144,840 cells, 3×4 = 12 ranks, steady state, coupled algorithm, implicit galvanostatic control |
| HPC production case | `~/coarse_08` (serpentine). It is huge, so **never run it locally**. Use it only as the dictionary template. |

### Git history
```
0cc4a1e electricControl: publish ibarErr = |ibar/target - 1| as a function-object result
b20a1e9 electroChemicalReaction: default sinkYmin 1e-10
1f8b096 ButlerVolmer: skip the zone rescale when a zone sum is <= VSMALL
cc0a1a0 ButlerVolmer: rate constant from the solved mass fraction (one solve less lag)
c0a9924 ButlerVolmer: zone-normalised rescale of the carried current (coupled)   (tag stageA)
2a9493e electroChemicalReaction: log the cells where the sink cap acts
8f301fd electroChemicalReaction: explicit water production
b2ca59e electroChemicalReaction: bound the implicit sink of consumed species
e1ad218 ButlerVolmer: any transfer coefficient 0 < alpha < 1 in the coupled solve
feb6458 alphaRhoPhi: start from zero when no file is present
c3a8f6e Restart: read the mass flux, eta and X; keep the volumetric flux current
9b28b91 electricControl: implicit galvanostatic control by superposition
37195d0 electricSystem: algorithm and operating control in constant/electricProperties
610042b Baseline: single-phase prune of openFuelCell2
```
Other branches (not merged): `nernstY` (3d35485, Nernst X from Y: no gain), `change3-unnormalised` (052cc6e, failed), `stageB-Ykinetics` (= master up to b20a1e9).

**Uncommitted:**
- `src/Allwmake` and `src/Allwclean` (modified to build the tool);
- `src/tools/regionProcAddressing/` (untracked);
- `report.md` (untracked, stale);
- this file.

The user has not yet said whether the tool should be committed. If it is, commit it on its own.

---

## 1. Architecture

### 1.1 Parent mesh and regions
- The solver runs on one **parent mesh** (the whole cell) and several **child regions** split from it.
- The addressing files `cellRegionAddressing`, `faceRegionAddressing` and `boundaryRegionAddressing` link each child to the parent. `regionType` reads them (`cellMapIO_` maps child to parent, `cellMap_` maps parent to child).
- Two regions only meet through parent ids. For example, an ion cell `i` maps to the electron cell `eR.cellMap_[iR.cellMapIO_[i]]`.

`constant/regionProperties`:
```
regions ( fluid (air fuel)  solid (electrolyte interconnect)  electric (phiE phiI) );
```
`regionTypeList::reset` builds the regions in the order fluid → solid → electric. It then builds one **`electricSystem`**, registered on `Time` as `"electricProperties"`, which owns both potential solves and the operating control. The kinetics models find it through `activationOverpotentialModel::electricSys()`.

| object | what it does |
|---|---|
| fluid region (`singlePhaseSystem`) | p–U (SIMPLE/SIMPLEC), species, electrochemical sources |
| solid region | solves nothing; supplies ρCp, k and heat sources to the parent |
| electric region (phiE, phiI) | holds φ, **i**, J, σ, T and the dissolved model; **solves nothing itself** |
| `electricSystem` (+ `coupledMatrix`, `electricControl`) | φE/φI solve (segregated or coupled) and the control |
| parent mesh (`fuelCellSystem` + `appSrc/EEqns.H`) | the single energy equation for T |

### 1.2 Main loop (`src/appSrc/myFuelCell.C`, unchanged from openFuelCell)
```cpp
fuelCell.mapToCell();
while (runTime.run()) {
    runTime++;
    fuelCell.mapFromCell();   // T: parent -> regions
    fuelCell.correct();       // properties, Nernst + kinetics (inside fluid correctThermo), control beforeSolve
    fuelCell.solve();         // fluid Y -> U -> p, then the electric system
    fuelCell.mapToCell();     // Qdot, rhoCp, k, rhoCpPhi -> parent
    #include "EEqns.H"        // T on the parent mesh
    runTime.write();
}
```
- With steadyState, one time step is one **outer iteration**. The physics lag each other by one iteration (block Gauss–Seidel).
- The coupled electric algorithm makes φE–φI–kinetics implicit within an iteration. The implicit control also makes V implicit.
- T, species (through `coeff`, built from the solved Y, so one species solve behind), the Nernst potential (X, two solves behind) and σ remain lagged.

### 1.3 Energy (parent mesh)
`ddt(ρcpT) + div(ρcp u T) − (div ρcp u) T − lap(k T) = Q`
- k is porosity-weighted in porous zones.
- Q = reaction heat `j(sign η − TΔS/(|n|F))` + Joule heat `|i|²/σ` + gravity work.
- There is no pressure work. The `dpdt` field and the parent `Cp` field were removed.

---

## 2. Fluid solve

- **Phase model.** The only phase model is `reactingPhaseModel` = Anisothermal< MultiComponent< Reacting< Moving< Thermo<phaseModel, rhoReactionThermo> >, CombustionModel > > >.
  - The thermo is heRhoThermo, perfect gas, JANAF, Sutherland.
  - The phase name is empty, so the fields are `U`, `alphaRhoPhi`, `p_rgh`, `X.<sp>`, `iDmdt<region>.<sp>`.
  - The phase fraction is still a field ≡ 1.
- **Ignored dictionaries.** The two-phase dictionaries in the case `regionProperties` (`water`, `blending`, `drag` …) are ignored, but they must stay for compatibility.
- **`fluid::correct()`.** It calls correctEnergyTransport, then **correctThermo**, then correct. `electroChemicalReaction::correct()` runs inside correctThermo: it computes Nernst, then Butler–Volmer, then the water activity and sources. **The kinetics use X from before it is updated from Y** (lagged one iteration).
- **`singlePhaseSystem::solve()`.** The order is `YEqns.H` → `UEqn.H` (+ porous Darcy/Forchheimer) → `pEqn.H` (SIMPLE) or `pcEqn.H` (SIMPLEC).
  - The flow is compressible and variable-density.
  - The flux in the equations is the mass flux `alphaRhoPhi` (kg/s).
- **Species.** Every species except `inertSpecie` (N₂ in air, H₂O in fuel) is solved:
  - equation: `ddt + div(alphaRhoPhi, Yi) − lap(ρD) = R_i`;
  - after the solve, Y is clipped to [SMALL, 1] and the inert species is set to 1 − ΣY.
  - D comes from `binaryFSG` in the channels and `porousFSG` (Knudsen + tortuosity) in the porous layers.
- **Electrochemical source** (`electroChemicalReaction::R()`, `splitSource(w, Y)`, b2ca59e, 8f301fd, b20a1e9).
  - The rate is `w_i = j W_i ν_i/(|ν_e| F)` (kg/m³/s); `iDmdt<region>.<sp>` stays the explicit w.
  - **Sign split:** production (w > 0) is explicit, `fvm::Su(w)`; consumption (w < 0) is implicit, `fvm::Sp(w/max(Y + SMALL, sinkYmin), Y)`, so the coefficient never exceeds |w|/sinkYmin.
  - Water uses the same split (8f301fd): production is explicit, and it gets the bounded implicit form only where the dissolved exchange wins (w < 0). The old `fvm::Sp(+w/Y)` that weakened the diagonal is gone.
  - **`sinkYmin`**: optional key in `constant/<region>/combustionProperties` (same dictionary as `dissolved`), **default 1e-10** (b20a1e9; it was 1e-6), printed once. Where it acts (w < 0 and Y + SMALL < sinkYmin) the sink is w·Y/sinkYmin, i.e. smaller than w; the cells are logged as `R(<region>/<specie>): N cells below sinkYmin` (twice per iteration: `dmdt()` and the species equation both call R()).
  - With 1e-10 no cell is capped up to 3.0 A/cm² on the straight channel; 1e-15 gives identical runs.
  - The total mass source `dmdt = ΣR_i` enters the pressure equation.
- **Pressure.** After the solve:
  - `phi = phiHbyA + flux`;
  - p_rgh is relaxed **before** p is updated (openFuelCell2 fix 7a4480b);
  - then correctRho;
  - **since c3a8f6e the volumetric flux is updated, `phase_.phiRef() = phi/interpolate(rho)`**. The written `phi`, the Courant number and the inletOutlet BCs (`phi phi;`) are no longer stale.

### Restart behaviour (c3a8f6e, feb6458)
| field | read on restart | cold start |
|---|---|---|
| `alphaRhoPhi` | READ_IF_PRESENT, AUTO_WRITE | starts at **zero** until the first pressure solve (feb6458). The log says "…not found: starting from zero". |
| `eta` (activation overpotential) | READ_IF_PRESENT | grows from 0 |
| `X.<sp>` | READ_IF_PRESENT (preserves the X lag across a restart) | computed from Y |

- A restart now follows the continuous run (V within ~2e-5). openFuelCell2's segregated restart collapses instead (ibar → 0).
- Use writePrecision ≥ 11–12. At 7, p ≈ 2e5 Pa is resolved to only 0.1 Pa, which gives a U/p_rgh residual kick on restart.

---

## 3. Electric system

### 3.1 Settings: `constant/electricProperties`
```
electronRegion  phiE;
ionRegion       phiI;
algorithm       $electricAlgorithm;   // segregated | coupled  (inputData key is electricAlgorithm, NOT algorithm)
segregatedCoeffs { ionReferenceRelax $r_i; }
coupledCoeffs    { omega 1; }         // Newton damping (0,1]
control
{
    type    galvanostatic;            // galvanostatic | potentiostatic
    patch   interconnect0;            // controlled Dirichlet patch of phiE
    zone    ccl;                      // phiE holds both electrodes -> integrate over ONE CL
    ibar    constant $ibar;           // potentiostatic: voltage <Function1>
    galvanostatic
    {
        mode    implicit;             // explicit | implicit (implicit = coupled only)
        rE      $r_E;                 // explicit gain; implicit: fallback
        relax   1;  maxStep 0.1;  Vmin 0;  Vmax 1.2;  // gMin optional
    }
}
```
- The coupled linear solver is `phiCoupled` in `system/phiI/fvSolution` (PCG + DIC, tol 1e-11, relTol 0).
- The old keys (`control`, `galvanostatic`, `relax`, `coupled` in `phiE/phiI regionProperties`) are **fatal errors**, and the message says where each key moved.
- Other configuration errors are also fatal at startup:
  - an unknown algorithm or mode;
  - a missing patch or zone;
  - implicit mode with the segregated algorithm;
  - Tafel/Constant kinetics with the coupled algorithm;
  - `omega` outside (0, 1];
  - in coupled mode, `alpha` outside (0, 1).
- The BV parameters sit in `constant/{air,fuel}/combustionProperties` → `activationOverpotentialModel` (`type ButlerVolmer; zoneName ccl|acl; j0; relax; alpha; gamma; species{… cRef}; nernstModel{standard; RxnList}`).
  - `relax` is used only in segregated mode. Coupled mode prints a one-time note that it is ignored.
  - The optional coupled keys are `inverseTolerance` (1e-12) and `inverseMaxIter` (100).

### 3.2 Equations
- **Charge conservation** in each region: `−∇·(σ∇φ) − J = 0`, `i = −σ∇φ`.
- **phiE** covers the BPP, GDL, MPL and CL of both electrodes. σ is `constantSigma` in the plates and `porousSigma` σ(1−ε) in the porous layers. The Dirichlet patches are `interconnect1` (0 V) and `interconnect0` (V).
- **phiI** covers both CLs and the membrane, all zero-gradient.
- **Kinetics** in the CL pairs: `η = φE − φI − E`, `J_I = G(η)`, `J_E = −J_I`. E comes from `nernstModels/standard` (ΔH, ΔS, partial pressures, `coeff = Π (xρ/(W cRef))^ξ`). The sign of n (negative at the cathode) fixes the source signs.

### 3.3 Segregated algorithm (original openFuelCell scheme)
1. **Explicit BV** during fluid correct:
   - `eta = (1−relax)eta + relax(φE−φI−E)`;
   - `j = max(j0·coeff·s^γ (exp(nαFη/RT) − exp(−n(1−α)Fη/RT)), 0)`;
   - the result goes into `J`.
2. **φE solve** with J frozen.
3. **φI solve** with J frozen. Because φI is all-Neumann, `setReference` is applied on every cell with the value `φ − iDot·ionReferenceRelax`. This doubles the diagonal (implicit under-relaxation 0.5) and drifts the level until ∫J_I = 0. Charge is conserved only at convergence.

### 3.4 Coupled algorithm (`electricSystem::solveCoupled`, `coupledMatrix`)
**Base point.** ButlerVolmer, coupled branch, in `correct()`. The state carried between iterations is the realised current `Jr` (the ion `J`), not η:
```
kin = max(j0·coeff·s^γ, SMALL·j0),  y = nFη/(RT)
solve h(y) = exp(αy) − exp(−(1−α)y) = sign(n)·Jr/kin        (butlerVolmerInverse.H)
eta* = y/nf,  dJdEta = kin·|nf|·h'(y),  etaE = eta* + E,  j = max(sign·Jr, 0)
```
- **Rate constant from Y (cc0a1a0, both algorithms).** `coeff = Π (Y_i ρ / W_i / cRef)^ξ` with the solved mass fraction Y_i and W_i = W(kg/kmol)/1000; this is the same molar concentration as X ρ/W, but one species solve fresher (X is updated in correctThermo after the kinetics). Nernst and the water activity still use X (branch `nernstY` tried Y there: no gain).
- **Zone-normalised rescale of the carried current (c0a9924).** Before the inverse, per catalyst zone (sums over all ranks):
  ```
  rho = kin/kinOld            (1 where kinOld = 0)
  a   = Σ V J / Σ V rho J,    J = sign(n)·SI (realised current of the last iteration)
  Jr  = a·rho·SI              written back to SI and SE (= −Jr)
  ```
  The zone total stays the realised total, so the species solve consumes exactly the cell current; only the distribution follows the new rate constants. `kinOld.<region>` is a volScalarField, READ_IF_PRESENT and AUTO_WRITE, so a restart continues the rescale. The `Total current … (coupled)` line ends with `, a = <a>`. (The unnormalised kin/kinOld rescale, branch `change3-unnormalised`, closed an oscillating O2-inventory loop and failed.)
- **Fallback (1f8b096).** If Σ V J or Σ V rho J of a zone is ≤ VSMALL, the rescale is skipped for that zone in that iteration: rho = 1, a = 1, Jr = SI; kinOld is still updated; one line `ButlerVolmer <zone>: rescale skipped (sum V J = …, sum V rho J = …)`. It appears only at iteration 1 of a cold start (no current yet).
- **Inverse (e1ad218).** The general-α inverse uses Newton–Raphson inside a bracket [0, Tafel limit `ln(1+r)/α`] (or the mirrored bracket for r < 0).
  - It falls back to bisection when a step leaves the bracket.
  - It uses `expm1` and stops on a relative step below tol.
  - Non-converged cells produce one warning per zone.
- **α convention.** α multiplies the dominant branch on both electrodes (on the cathode n < 0 and η < 0).
- **Regression.** α = 0.5 reproduces the old asinh closed form to ~1e-9 V.

**Linearised system** (SPD), with `r = ω V G'(η*)` per cell pair:
```
electron row: (aP_E + r)φE − Σ aN φE,N − r φI = −V J* + r(η*+E)
ion row:      (aP_I + r)φI − Σ aN φI,N − r φE = +V J* − r(η*+E)
```
- φI needs no reference, because the coupling ties it to φE.
- The run prints `Coupled: pairs = …, max|eta_new − eta*| = …, I+ = …, I- = …, rel. imbalance = …`.

**`coupledMatrix`.** Built once:
- the unknown numbering `[φE | φI]`;
- the cell pairs;
- the faces and upper-triangular order;
- the processor interfaces (offset faceCells, reusing the φ halo exchange);
- the `lduPrimitiveMesh` and `lduMatrix`.

Every iteration runs `assemble(dJdEta, etaE)` → `solve(psi, rhs, controls)` (repeatable on the same matrix) → `scatter` → `realise(etaE)`.

**Realise.** `Jreal = J* + (r/V)(φE − φI − etaE)`. It is written to `iR.J` and `eR.J = −Jreal`, so ∫J_I = 0 to solver precision. Then i = −σ∇φ is computed and the dissolved model is solved.

**Parallel requirement.** Every pair must be on one processor. `decomposeParID` guarantees this. Coupled patches other than processor patches are rejected.

**Why it is robust.** The iteration is Newton on the *inverse* kinetics, so a big η step gives only a linear current step. ω = 1 is safe, and ω does not change the fixed point.

### 3.5 One outer iteration (coupled algorithm, implicit control)
```
mapFromCell   T -> regions
correct       air/fuel correctThermo: Nernst, BV inverse (Jr -> eta*, dJdEta, etaE, j), water/dissolved
              phiE/phiI sigma, dissolved correct
              electricControl::beforeSolve   -> log "ibar: <realised> voltage: <V_old>"
solve         air/fuel: Yi -> U -> p (dmdt) -> phi = alphaRhoPhi/rho_f
              electric: assemble M,b0; solve1 M phi1 = b0 (patch at V_old)
                        afterSolve: unitSource c, solve2 M phic = c, ibar1, g, dV;
                                    psi = phi1 + dV phic; patch = V_old + dV
                        scatter; realise; i_E, i_I; dissolved
mapToCell     Qdot, rhoCp, k, fluxes;  EEqns: T
```

---

## 4. Operating control (`electricControl`)

The measured quantity is `ī = Σ_{zone} J_E V / A_patch`.
- Use the zone `ccl`: over the whole of phiE, the anode and cathode currents cancel.
- Every mode prints `ibar: … voltage: …` once per iteration (ibar = realised zone current, V = patch value before this iteration's update). `post_VI.py` and `compare.py` parse this line.
- **`ibarErr` (0cc4a1e).** At that line, galvanostatic control stores `ibarErr = |ibar/target − 1|` as the function-object result `ibarErr` of object `electricControl` (Time's function-object properties, the store runTimeControl/convTests read). Potentiostatic runs have no target and publish nothing.

| mode | where | update |
|---|---|---|
| potentiostatic | beforeSolve | V = V(t) (Function1); both algorithms |
| explicit galvanostatic | beforeSolve | `V += rE (ī − ī*)` with the previous J. This is an integral controller with loop gain `rE·|dī/dV|` (≈0.34 at rE = 1e-5). Unstable near gain 2 (the openFuelCell2 `full_thick` failure). Works with both algorithms. |
| **implicit galvanostatic** | afterSolve, coupled only | superposition, below |

**Implicit control.** Within one iteration `Mφ = b0 + V c` is linear in V, so `φ(V) = φ1 + (V − V_old)φc`.
1. `c = unitSource(patch)` is the boundaryCoeffs of the same `−laplacian` with the patch set to 1. Nothing about the patch is hand-coded.
2. Solve 2, `M φc = c`, is logged as `phiUnit` and warm-started from the previous φc.
3. `ī1` and `g = dī/dV` are taken from the tangent current with the same zone sum. Expect g < 0; `−1/g` is the fast area-specific resistance (≈2.9e-5 Ω m²).
4. `ΔV = relax·(ī* − ī1)/g`, clamped to `|ΔV| ≤ maxStep` and `V ∈ [Vmin, Vmax]`.
5. `psi += ΔV·φc`, then the patch is set and scatter/realise run as usual.
6. **Fallback.** If `g ≥ −gMin`, the step is skipped and the explicit update with rE is used, with a warning.

**Log line:**
```
Implicit control: ibar1 = …, g = … A/m2/V, -1/g = … Ohm m2, dV = …, V = …, ibar = … [(limited)]
```

**Judging convergence.** With an attainable target, ī = ī* holds by construction from about iteration 9, so the current error alone says little; the species balances (errFar) are the slow quantity. Above the limiting current, however, V sits at the Vmin clamp and the realised current stays short of the target (tortuosity 3 at 2.5 A/cm²: 1.8 %), while errFar can still pass, because it uses the realised current. The stop rule therefore needs both (§5).

**Start-up.** A cold start has J = 0, so |g| is tiny. V dives towards about 0.1–0.3 V, gets clamped by maxStep, and overshoots ī for a few iterations. Warm starts avoid this (§6).

---

## 5. Case setup used now (`warmStart/base`, coarse_08 dictionaries)

- **inputData:**
  - operating point: p 2e5, T 473.15, ibar (A/m²);
  - inlet composition: yO2 0.23, yN2 0.769, yH2O_ca 0.001, yH2 0.99, yH2O_an 0.01;
  - kinetics and conductivity: j0_ca 47000, j0_an 1e8, σ (bpp 83300, gdl/cl 150, mem 22);
  - electric settings: `electricAlgorithm coupled`, r_E 0.5e-5, r_i −4e-12, r_o/r_f 0.001 (BV relax);
  - geometry and flow: A_act 142e-6, A_an = A_ca 1e-6, l_an = l_ca 2, nx 3, ny 4;
  - derived with `#eval`: NPROCS, U_an/U_ca, cRef_an/cRef_ca.
- **BV:** cathode zone `ccl`, `alpha 0.5` literal, gamma 1.5, `cRef $cRef_ca`. The old case bug, where air used `$cRef_an`, is fixed in this base.
  - The `alphaTest` cases add the inputData keys `alpha_ca`/`alpha_an`.
  - The HPC calibration "set B" uses j0_ca 3300, alpha 0.6, mem 20.
- **Species scheme (2026-11-06), air and fuel:** `div(alphaRhoPhi,Yi) Gauss linearUpwind limitedGrad;` with `limitedGrad cellLimited Gauss linear 1;`.
  - The old `Gauss limitedLinear01 1` gave a period-2 flip-flop of O₂/H₂O (errFar ±1e-3 forever).
  - On a cold start it also gave an O₂ slug at iteration 1.
  - U keeps `Gauss limitedLinearV 1`.
- **Stop rule (2026-11-18, the same everywhere: local, sweep, fit):** decided inside the run by the coded function object **`convTests`** (`system/convTests`, `#include`d in `controlDict.run`, listed after airBalance and fuelBalance).
  - It reads `O2_errFar`, `H2O_errFar` (airBalance), `H2_errFar` (fuelBalance) and `ibarErr` (electricControl, §4). An iteration counts if all four |values| < `tol` 1e-3; any failure, and any missing result, resets the count to 0.
  - At `nConsecutive` 10 it prints `convTests: converged at iteration <n>`, writes the fields of that iteration and ends the run (`Time::writeAndEnd`). Every iteration it logs `convTests: <count> of 10` (plus `missing <fo>/<result>` if a name is wrong).
  - It is the `nBelowTol` counter of `run/.../convTests/system/speciesBalance`, moved into its own function object, generalised to any result, with missing = not satisfied and its own stop. No runTimeControl is used any more (in v2312 runTimeControl counts a missing result as satisfied).
  - The coded `speciesBalance` (one per fluid region) is the convTests version: the `inertSpecie` gets its patch flux as `alphaRhoPhi` minus the other species (it has no diffusivity). It publishes `<sp>_errFar/_errTra/_errSrc/_res`.
  - `O2ccl` (coded): min and volume average of O2 over zone `ccl` of air every iteration (`O2ccl: min … mean …`, `postProcessing/air/O2ccl/<first it>/O2ccl.dat`). `libfieldFunctionObjects` (volFieldValue, fieldMinMax) cannot be loaded into myFuelCell: its `phaseSystem`/`saturationModel` classes clash with the solver's ("Duplicate entry function1", then a failed `phaseProperties` lookup).
  - On the straight channel the new rule stops w25 at it 43 and c08 at it 34 (the errFar-only rule gave the same iterations; ibarErr at the w25 stop 8.9e-16). At that stop V is still 2–5e-4 V above the 1e-4-converged V at 2.5–2.9 A/cm².
- **Pre-processing (`base/Allpre`, one log `log.pre`):** ofc2mesh mesh (`ssc_thick.yaml`) → `decomposeParID` → `parallel.sh` (regions split **in parallel after decomposition**) → setFields → controlDict.run.
  - `Allsolve` runs the solver to `log.run`. `Allrun` does both.

---

## 6. Workflows that exist around the solver

- **Warm start / polarization curves.** Directory `run/H2/HTPEMFC/warmStart/`; read its `CLAUDE.md` first.
  - `warmStart.sh <src> <dst> <ibar> [time]` copies the decomposed `processorsN/<time>` to `processorsN/0`, sets ibar, and resets only the inlet U with parallel collated `changeDictionary -subDict dictionaryReplacement -enableFunctionEntries`.
  - Inlet Y, T and p are **not** reset.
  - `polCurve.sh <dir> <i0> <di> <i1>` runs a chain up or down, with a watchdog.
  - `compare.py` summarises a run.
  - **sinkFix** (`run/H2/HTPEMFC/sinkFix/`, `RESULTS.md`): the O2-starvation work behind b2ca59e…b20a1e9. Local chains: `hc/polCurve.sh`, `hc/single.sh` (watchdog only; convergence from the convTests line). With master the straight channel converges warm up to 3.0 A/cm² (period-2 V mode, decaying ×1.8–2.4 per 10 its); tortuosity 3 at 2.5 A/cm² is above the limiting current (V at Vmin from it 47, never converges).
  - **Results (straight channel):** cold 0.8 A/cm² converges in 22 iterations. Warm steps of 0.1 A/cm² converge from 0.1 to 2.5 A/cm² in 12–48 iterations. 2.6 fails from O₂ starvation (min Y_O2 in the CL under the rib at the outlet end wall hits the clip).
  - Cold starts fail at 0.1, 0.2 and 3.0.
- **Changing the processor count / reconstructing.** Directory `run/H2/HTPEMFC/redecompose/`; read its `CLAUDE.md` first.
  - Plain `reconstructPar/decomposePar -region` is wrong here, because the processor region meshes come from parallel.sh, not from decomposePar.
  - Run `regionProcAddressing -allRegions` (serial), then `reconstructPar -allRegions` (`-withZero` for time 0), then `decomposePar -allRegions -fields`.
  - Script: `reDecompose.sh <src> <dst> <nx> <ny> [time]`.
  - A cold start on 2×4 crashes at iteration 7 (PBiCGStab FPE on H2O); this was not investigated.
- **α tests.** `run/H2/HTPEMFC/alphaTest`, with `fieldcmp.py` and `physcheck.py`.
  - Cold 0.8 A/cm²: V = 0.64495 at α = 0.5, 0.72068 at α = 0.6 (+75.7 mV), and 0.53125 at α = 0.4.
- **Regression (A/B).** `run/H2/HTPEMFC/myFuelCellTest/` with `run.sh` and `cmp.sh A B`.
  - Cold-start references for the current code are `z_{seg,cpl}`.
  - Re-run the A/B check after every refactor.
- **HPC.** Polarization runs of coarse_08 (AnodeAnode) and polCurve_ac1 (AnodeCathode) on ARIS use the same scripts with `MPIRUN=srun`. Their state is kept outside this repo (`~/i04000/OBSERVATIONS.md` for the 2026-11-12 crash at 0.4 A/cm² with set B).
  - Sweep: `~/polCurve_coarse08` (`POLCURVE_SWEEP.md`): two SLURM chains, up 8000 → 30000 then down 7000 → 1000 (`afterany`); a point is converged only if log.run has the convTests line.
  - Fit: `~/fit_coarse08` (5 parameter sets × 9 currents), same function object and criteria once its controlDicts are switched to convTests.

---

## 7. Known limitations / open items

1. **One electron region and one ion region.** Separate phiEA/phiEC must be merged into one phiE. The control zone is then a single CL.
2. **Coupled algorithm:**
   - Butler–Volmer only, with 0 < α < 1;
   - processor patches only;
   - pairs must be co-located (`decomposeParID`).
3. **Implicit control:**
   - coupled only;
   - a cold start from J = 0 gives a tiny |g| and clamped steps;
   - at low current the implicit controller can bang-bang between ±maxStep (BV is steep at high V) and blow up T;
   - a lagged φc refresh is not implemented.
4. **H₂O source sign.** Done (8f301fd): production is explicit.
5. **X lag.** The BV rate constant now uses Y (cc0a1a0); Nernst and the water activity still use X from the previous iteration. Nernst from Y (branch `nernstY`) gave no gain.
6. **Forced dissolved model.** `electroChemicalReaction` requires a `dissolvedModel` (HT-PEM: `dissolveOnOff true` with `type none`).
7. **Phase fraction.** The phase fraction is a field ≡ 1. It multiplies Y and U and enters BV as s^γ.
8. **Slow species mode.** The H2 errFar decays slowly when species relaxation is < 1.
9. **CL edge effect.** The last x-column of the CL next to the end wall carries about 2× the j of its neighbours; not investigated.
11. **High current.** Above ~2.6 A/cm² (straight channel) a period-2 V/O2 mode decays only ×2 per 10 iterations; at the convTests stop V can be a few 1e-4 V from the fully converged value. Points above the limiting current end at the Vmin clamp and do not converge (by design of the stop rule).
10. **Stale doc.** `report.md` §3.4/§6 are stale (alpha 0.5, stale phi). This file and the Claude doc are current.

**Planned "step 3"** (numerics, each to be verified on its own): forced dissolvedModel, and alpha ≡ 1. The stale phi item was done in c3a8f6e, the H2O `fvm::Sp` sign in 8f301fd.

---

## 8. Where to look (paths relative to `src/libSrc/fuelCellSystems/`)

| topic | file |
|---|---|
| main loop, energy | `src/appSrc/myFuelCell.C`, `src/appSrc/EEqns.H` |
| parent system, maps | `fuelCellSystem/fuelCellSystem.C`, `regions/region/regionType*.C`, `regionTypeList.C` |
| fluid, SIMPLE | `regions/fluid/fluid.C`, `solvers/singlePhaseSystem/{singlePhaseSystem.C,YEqns.H,UEqn.H,pEqn.H,pcEqn.H}` |
| phase chain | `phaseModel/` (`MovingPhaseModel` for alphaRhoPhi/phiRef, `MultiComponentPhaseModel` for YiEqn/X/dmdt) |
| electrochemical sources | `electroChemicalModel/electroChemicalReaction.C` (`splitSource`, `sinkYmin`) |
| kinetics | `activationOverpotentialModels/ButlerVolmer/{ButlerVolmer.C,butlerVolmerInverse.H}` (coeff from Y, rescale, fallback), heat in `activationOverpotentialModel.C` |
| Nernst | `nernstModels/standard/standard.C` |
| electric regions | `regions/electric/electric.C` |
| electric system | `electricSystem/electricSystem.C`, `coupledMatrix/coupledMatrix.C`, `electricControl/electricControl.C` (ibarErr) |
| stop rule, balances | case files `system/convTests`, `system/speciesBalance`, `system/controlDict.run` (`run/H2/HTPEMFC/warmStart/base`) |
| σ, D, porosity, dissolved, saturation | `sigmaModels/`, `diffusivityModels/`, `porosityModel/`, `dissolvedModels/`, `saturationModels/` |
| re-decompose tool | `src/tools/regionProcAddressing/` (uncommitted) |

The longer written report on the implicit galvanostatic control math and the BV Newton inverse is a Claude doc: https://claude.ai/artifact/VGczbkm7ZChvoKJNJt9fZ8

## 9. Working rules from the user
- Ask before major changes: solver code, base dictionaries, long runs.
- Work step by step, and keep the case dictionaries compatible unless the user agrees otherwise.
- Test locally on straight channels only.
- Don't `pkill` (shared machine). Stop a run with `foamDictionary -entry stopAt -set writeNow <case>/system/controlDict`.
- Use `LC_ALL=C` for awk/printf (the Greek locale prints decimal commas).
- Source the OpenFOAM bashrc before `set -e`.
- Keep few log files per case (`log.pre`, `log.decompose`, `log.run`).
