# myFuelCell: how the solver works

**Commit:** `9b28b91` ("electricControl: implicit galvanostatic control by superposition").
Earlier commits: `37195d0` (electricSystem, settings in `constant/electricProperties`) and `610042b` (baseline, a single-phase prune of openFuelCell2).
**Forked from:** openFuelCell2, branch `coupled-phi-parallel`, commit `7a4480b`. That code is described in `openFuelCell2/report.md`.
**OpenFOAM:** v2312 (COM). The default shell on the workstation loads v2406, so `source /usr/lib/openfoam/openfoam2312/etc/bashrc` first.
**Reference case used for the examples:** `run/H2/HTPEMFC/myFuelCellTest/cpl8000_implicit`. It is the HT-PEMFC `ssc_1` geometry on 24 ranks, steady state, coupled algorithm, implicit galvanostatic control at 8000 A/m².

All paths below are relative to `src/libSrc/fuelCellSystems/` unless they say otherwise.

---

## 0. What myFuelCell is

myFuelCell is the part of openFuelCell2 that the HT-PEMFC robustness cases actually use, plus a restructured electric solve:

- **One phase only.** Every fluid region is a `singlePhaseSystem` with a `reactingPhaseModel`.
- **Two electric algorithms**, chosen in `constant/electricProperties`: `segregated` (the original openFuelCell scheme) or `coupled` (a single φE–φI solve with Newton-linearised kinetics).
- **Three control modes:** potentiostatic, explicit galvanostatic, and implicit galvanostatic (coupled only).
- **All material and transport models are kept:** diffusivity, conductivity, activation overpotential, Nernst, dissolved water, porosity and saturation.

| | openFuelCell2 | myFuelCell |
|---|---|---|
| Application | `openFuelCell` | `myFuelCell` |
| Libraries | `libfuelCellSystem` + 7 helper libraries | `libmyFuelCellSystem`, `libmyRegionCourantNo`, `libmyFuelCellTurbulenceModels` |
| Source size | 83.8k lines | ≈21k lines |
| Fluid systems | single-phase, two-phase, drift-flux | single-phase |
| Electric settings | `phiE/regionProperties` (`control`, `galvanostatic`), `phiI/regionProperties` (`relax`, `coupled {}`) | `constant/electricProperties` |
| Galvanostatic update | explicit integral controller | explicit, or implicit (superposition of two linear solves) |
| Meshing tools | `decomposeParID`, `topoSet`, `libmasterRegionToCell` | not built here; the openFuelCell2 binaries are used |

**Removed:**
- the two-phase and drift-flux systems and the `PhaseSystems` templates;
- population balance, diameter models, phase pairs, blending, and all interfacial models (drag, lift, …);
- the interfacial composition models except saturation;
- the custom boundary-condition libraries (`derivedFvPatchFields`, `thermoTools`);
- the Inert/Pure/Stationary/Isothermal phase models;
- the `phaseChange` branches, and the Nernst instantiations whose "other" phase is a liquid;
- `ChangJaffe` and `generalSigma`, which were never compiled in openFuelCell2 either.

None of the removals changes the numerics. Segregated, coupled and potentiostatic runs reproduce openFuelCell2 bit for bit (§5).

---

## 1. Overall architecture

### 1.1 Parent mesh and child regions

The region structure is the same as in openFuelCell2. The solver works on one **parent mesh** (the whole cell) and several **child regions** cut out of it with `splitMeshRegions`/`topoSet`. The addressing files `cellRegionAddressing`, `faceRegionAddressing` and `boundaryRegionAddressing` link each child to the parent. They are read by `regionType` ([regions/region/regionType.C](src/libSrc/fuelCellSystems/regions/region/regionType.C)):

| member | meaning |
|---|---|
| `cellMapIO_[childCell]` | parent cell id |
| `cellMap_` (hash) | parent cell id → child cell id |
| `faceMap_`, `faceMask_` | parent face id and orientation sign for every child face |
| `patchesMapIO_` | parent patch id for each child patch |

Two regions only meet through parent ids. For example, ion cell `i` → electron cell is `eR.cellMap_[iR.cellMapIO_[i]]`.

`constant/regionProperties`:

```
regions
(
    fluid    (air fuel)
    solid    (electrolyte interconnect)
    electric (phiE phiI)
);
```

`regionTypeList::reset` ([regions/region/regionTypeList.C](src/libSrc/fuelCellSystems/regions/region/regionTypeList.C)) builds the regions in the order **fluid → solid → electric**, whatever order the dictionary uses. If an `electric` entry exists, it then builds one **`electricSystem`**, which owns the potential equations and the operating control of the two electric regions (§3):

| object | class | what it does |
|---|---|---|
| `fluid` region | `regionTypes::fluid` → `singlePhaseSystem` | p–U (SIMPLE/SIMPLEC), species, electrochemical sources |
| `solid` region | `regionTypes::solid` | nothing is solved; supplies ρCp, k and heat sources to the parent |
| `electric` region | `regionTypes::electric` | holds φ, **i**, J, σ, T and the dissolved-water model. It **no longer solves anything** |
| electric system | `electricSystem` (+ `coupledMatrix`, `electricControl`) | solves φE and φI (segregated or coupled) and applies the control |
| parent mesh | `fuelCellSystem` + `appSrc/EEqns.H` | the single energy equation (T) |

`electricSystem` is registered on the `Time` database under the name `electricProperties`. The activation-overpotential models look it up there (`activationOverpotentialModel::electricSys()`).

### 1.2 Main loop

[src/appSrc/myFuelCell.C](src/appSrc/myFuelCell.C) is unchanged from openFuelCell:

```cpp
fuelCell.mapToCell();                    // initial properties -> parent
while (runTime.run())
{
    runTime++;
    (Courant / LTS time-step control)
    fuelCell.mapFromCell();              // 1. T: parent -> every region
    fuelCell.correct();                  // 2. properties, kinetics, operating control
    fuelCell.solve();                    // 3. region PDEs, then the electric system
    fuelCell.mapToCell();                // 4. Qdot, rhoCp, k, rhoCpPhi -> parent
    #include "EEqns.H"                   // 5. energy equation on the parent mesh
    runTime.write();
}
```

`regionTypeList::correct()` and `solve()` first loop over the regions, then call `electricSystem::correct()` and `electricSystem::solve()`. The order of operations is therefore exactly the openFuelCell order.

With `ddtSchemes default steadyState;` one time step is one **outer iteration**. The physics are coupled by lagging each other by one outer iteration (block Gauss–Seidel). The coupled electric algorithm (§3.4) makes the φE–φI–kinetics part implicit. The implicit control (§4.3) also makes the cell voltage implicit within the iteration.

| phase | fluid regions (air, fuel) | electric regions + electric system |
|---|---|---|
| `mapFromCell` | T → `he(p,T)` | T copied into each electric region |
| `correct` | transport, `correctThermo` → **Nernst, Butler–Volmer, water sources** | regions: σ models, dissolved `correct`; system: **`electricControl::beforeSolve()`** |
| `solve` | Y → U → p (§2) | regions: nothing; system: φE, φI segregated (§3.3) **or** one coupled solve (§3.4), with the implicit control inside it |
| `mapToCell` | reaction heat, ρCp, k_eff, ṁCp fluxes | Joule heat $\mathbf i\cdot\mathbf i/\sigma$ |

### 1.3 Energy equation (parent mesh)

The energy equation is unchanged ([src/appSrc/EEqns.H](src/appSrc/EEqns.H)):

$$
\partial_t(\rho c_p T) + \nabla\cdot(\rho c_p \mathbf u\,T) - (\nabla\cdot \rho c_p\mathbf u)\,T - \nabla\cdot(k\nabla T) = \dot Q
$$

- `rhoCpPhi` comes from the fluid regions' `alphaRhoPhi·Cp` and is zero in the solids.
- k is porosity-weighted in the porous zones, $k_s(1-\varepsilon)+k_f\varepsilon$.
- $\dot Q$ is the sum of the electrochemical heat $j\,(\mathrm{sign}\,\eta - T\Delta S/(|n|F))$ in the catalyst layers, the Joule heat $|\mathbf i|^2/\sigma$, and the gravity work $\rho\,\mathbf u\cdot\mathbf g$.

There is no pressure-work term, so the `dpdt` field of openFuelCell had no consumer and was removed. The parent-mesh `Cp` field, which was always zero, was removed as well.

---

## 2. The fluid solve

### 2.1 Construction

`regionTypes::fluid` ([regions/fluid/fluid.C](src/libSrc/fuelCellSystems/regions/fluid/fluid.C)) owns a `phaseSystem`. The only concrete system is `singlePhaseSystem` (`type singlePhaseSystem;` in `constant/<region>/regionProperties`). Its phase is the only registered phase model, `reactingPhaseModel`:

```
AnisothermalPhaseModel<
  MultiComponentPhaseModel<
    ReactingPhaseModel<
      MovingPhaseModel< ThermoPhaseModel<phaseModel, rhoReactionThermo> >,
      CombustionModel<rhoReactionThermo> > > >
```

What each layer provides, after the prune:

| layer | provides |
|---|---|
| `phaseModel` | the phase fraction field (a `volScalarField`, ≡ 1), name, index |
| `ThermoPhaseModel` | `rhoReactionThermo` (here `heRhoThermo`, perfect gas, `multiComponentMixture`, JANAF, Sutherland) |
| `MovingPhaseModel` | U, the mass flux `alphaRhoPhi`, laminar turbulence model, `UEqn()`, continuity error |
| `ReactingPhaseModel` | the combustion model `electroChemicalReaction` (§2.4): `R(Y)`, `Qdot` |
| `MultiComponentPhaseModel` | `YiEqn`, mole fractions X, `iDmdt`, the diffusivity models, `dmdt()` |
| `AnisothermalPhaseModel` | `thermo.correct()` |

Notes:
- In single-phase mode the phase name is empty. The fields are therefore `U`, `alphaRhoPhi`, `p_rgh`, `iDmdt<region>.<species>`, `X.<species>`, …, without a phase suffix.
- The phase-level enthalpy equation (`heEqn`) and its pressure-work source (`heQdot`) were removed. Energy is only solved on the parent mesh.
- The two-phase dictionaries still found in the case `regionProperties` files (`water`, `blending`, `drag`, …) are **ignored**. They are kept only so that the robustness cases stay compatible.

Porous zones (GDL, MPL, CL) are `porousZone` objects (`porosityModel/`: DarcyForchheimer, powerLaw, fixedCoeff). They carry the porosity, the solid k and Cp, and the momentum resistance.

### 2.2 `fluid::correct()`

```cpp
phases_->correctEnergyTransport();
phases_->correctThermo();   // -> ReactingPhaseModel::correctThermo -> reaction_->correct()
phases_->correct();
```

**The electrochemistry is evaluated here**, before any electric solve. `electroChemicalReaction::correct()` calls the activation-overpotential model (§3.2), which writes the current source into the electric regions' `J` fields (segregated), or the base point of the linearisation (coupled). This is why the fluid regions come first.

### 2.3 `singlePhaseSystem::solve()`

[solvers/singlePhaseSystem/singlePhaseSystem.C](src/libSrc/fuelCellSystems/solvers/singlePhaseSystem/singlePhaseSystem.C):

```cpp
while (pimple_.loop())                    // nOuterCorrectors (1)
{
    correct();                            // phase_.correct(): continuity error
    #include "YEqns.H"                    // species
    #include "UEqn.H"                     // momentum (+ porous resistance)
    while (pimple_.correct())             // nCorrectors (1)
        consistent ? "pcEqn.H" : "pEqn.H" // SIMPLEC / SIMPLE pressure
    if (turbCorr) correctTurbulence();
}
```

The reference case uses `momentumPredictor yes`, `consistent false` (so SIMPLE, `pEqn.H`), relaxation `p_rgh 0.3`, `U 0.7`, `"(O2|H2O).*" 0.99` (air) and `"H2.*" 0.99` (fuel). The flow is **compressible and variable-density**. The flux `phi` in the equations is `phase_.alphaRhoPhiRef()` in kg/s.

**Species** (`MultiComponentPhaseModel::YiEqn`). Each active species, i.e. every species except `inertSpecie` (N₂ in air, H₂O in fuel), is solved in turn:

$$
\partial_t(\alpha\rho Y_i) + \nabla\cdot(\alpha\rho\mathbf u Y_i) - \nabla\cdot\big(\overline{\alpha\rho}\,\overline{D_i}\,\nabla Y_i\big) = \alpha R_i
$$

The equation is relaxed, solved with the `Yi` solver, and $Y_i$ is clipped to $[\text{SMALL},1]$. The inert species is set to $1-\sum Y$ in `correctThermo`, all mass fractions are normalised, and X is updated. The diffusivities $D_i$ come from `diffusivityModels` (`binaryFSG` in the channels, `porousFSG` with Knudsen and tortuosity in the porous layers).

**Momentum** (`UEqn.H`): `phase_.UEqn()`, which is ddt + div + continuity-error SuSp + turbulence. `porousZone_->addResistance(UEqn)` adds the Darcy/Forchheimer resistance, then the equation is relaxed and, with `momentumPredictor yes`, solved against $-\nabla p_{rgh}$.

**Pressure** (`pEqn.H`):

$$
\partial_t\rho + \psi\,\partial_t p'_{rgh} + \nabla\cdot\phi_{HbyA} - \nabla\cdot\big(\overline{\rho\,r_{AU}}\,\nabla p_{rgh}\big) = \dot m
$$

- **$\dot m=\sum_i R_i$** (`MultiComponentPhaseModel::dmdt()`) is the net mass the electrochemistry adds to the gas (O₂ consumption and H₂O production at the cathode, H₂ consumption at the anode). It couples p–U to the current and the species.
- After the solve, `phi = phiHbyA + pEqn.flux()`, and then p_rgh is **relaxed before p is updated**: `p_rgh.relax(); p = max(p_rgh+ρgh, pMin); p_rgh = p−ρgh`. This is the openFuelCell2 fix `7a4480b`.
- Then `thermo.correctRho(ψΔp)` and the velocity correction follow.
- SIMPLEC (`pcEqn.H`, `consistent true`) is still available.

### 2.4 Electrochemical sources in the species equations

`electroChemicalReaction` ([electroChemicalModel/electroChemicalReaction.C](src/libSrc/fuelCellSystems/electroChemicalModel/electroChemicalReaction.C)):

$$
w_i = \frac{j\, W_i\, \nu_i}{|\nu_e|\,F},\qquad R_i = \texttt{fvm::Sp}\!\left(\frac{w_i}{Y_i+\text{SMALL}},\,Y_i\right)
$$

- $j\ge 0$ is the volumetric reaction rate (A/m³) from the activation model, and $\nu_i$ are the `RxnList` coefficients of the Nernst model.
- For a consumed species the implicit coefficient strengthens the diagonal. For a produced species (H₂O) it **weakens** it, most of all where $Y_{H_2O}$ is small. This is one of the open numerical items (§6).
- Water is handled in `correct()`:
  - the water activity $a = x_{H_2O}\,p/p_{sat}(T) + 2(1-\alpha)$ is computed in the catalyst-layer cells (saturation model, here ArdenBuck) and passed to the ion region's `dissolvedModel`;
  - the water source is the production minus the uptake into the dissolved phase, and is stored in `iDmdt(H2O)`;
  - with `dissolvedModel none` the uptake is zero.

In myFuelCell the water always stays in the same (single) phase. The openFuelCell2 branch that moved it to a second, liquid phase is gone.

### 2.5 Mapping back

`fluid::mapToCell` sends the following to the parent mesh:
- the heat source (reaction heat and gravity work);
- ρCp and the continuity error·Cp;
- `alphaRhoPhi` and `alphaRhoPhi·Cp` on all faces, including boundary faces;
- the porous k_eff.

`mapFromCell` recomputes `he` from the parent T at fixed p.

---

## 3. The electric system

### 3.1 Settings: `constant/electricProperties`

```
electronRegion  phiE;
ionRegion       phiI;

algorithm       coupled;            // segregated | coupled

segregatedCoeffs
{
    ionReferenceRelax  -4e-11;      // ion potential reference (segregated)
}

coupledCoeffs
{
    omega   1;                      // Newton damping, (0, 1]
}

control
{
    type    galvanostatic;          // galvanostatic | potentiostatic
    patch   interconnect0;          // controlled Dirichlet patch of phiE
    zone    ccl;                    // current collected over this zone (optional)
    ibar    constant 8000;          // potentiostatic: voltage <Function1>

    galvanostatic
    {
        mode    implicit;           // explicit | implicit
        rE      1e-5;               // explicit gain; implicit: fallback
        relax   1;                  // implicit only
        maxStep 0.1;                // [V]
        Vmin    0;                  // [V]
        Vmax    1.2;                // [V]
        gMin    1e-8;               // optional [A/m2/V]
    }
}
```

- Only the coefficient dictionary of the selected algorithm is read.
- The coupled linear solver is read from `system/<ionRegion>/fvSolution`, entry `phiCoupled` (PCG + DIC, tolerance 1e‑11, relTol 0 in the reference case).
- In the segregated algorithm each region is solved with its own `phi` solver.
- In the robustness inputData the algorithm variable is called `electricAlgorithm`. It is deliberately not called `algorithm`, because inputData is `#include`d in many dictionaries.

The electric regions keep only region data in `constant/<region>/regionProperties`: `sigma {}`, `dissolveOnOff`, `dissolved {}` and `T0`. The old keys `control`, `galvanostatic`, `relax` and `coupled` are **fatal errors**, and the message says where each one moved. Other configuration errors also stop the run at startup with a message:
- an unknown `algorithm` or `mode`;
- a missing control patch or zone;
- implicit mode with the segregated algorithm;
- Tafel/Constant kinetics with the coupled algorithm;
- `omega` outside (0, 1].

### 3.2 Governing equations and kinetics

Each electric region carries a potential φ with steady charge conservation:

$$
-\nabla\cdot(\sigma\nabla\phi) - J = 0,\qquad \mathbf i = -\sigma\nabla\phi
$$

- **phiE (electron region):** BPP, GDL, MPL and CL of both electrodes. σ is `constantSigma` in the plates and `porousSigma` = σ(1−ε) in the porous layers. It has Dirichlet patches `interconnect1` (0 V) and `interconnect0` (the cell voltage V, set by the control).
- **phiI (ion region):** the two CLs and the membrane, with `constantSigma` (or `lambdaSigma`). All its boundaries are zero-gradient.

The kinetics link the two regions cell by cell inside each CL:

$$
\eta = \phi_E - \phi_I - E,\qquad J_I = G(\eta),\qquad J_E = -J_I
$$

The Nernst potential $E$ is computed by `nernstModels/standard` from ΔH, ΔS and the local partial pressures. The concentration factor is $\text{coeff}=\prod_k (x_k\rho/(W c_{ref,k}))^{\xi_k}$. The sign of the electron stoichiometry $n$ (negative at the cathode) fixes the sign of the source: at the cathode $J_E=+j$, $J_I=-j$, and the reverse at the anode.

### 3.3 Segregated algorithm

**Butler–Volmer (segregated branch).** `ButlerVolmer::correct()` ([activationOverpotentialModels/ButlerVolmer/ButlerVolmer.C](src/libSrc/fuelCellSystems/activationOverpotentialModels/ButlerVolmer/ButlerVolmer.C)) runs during `fluid::correct()`. For every CL cell (electron cell `e`, ion cell `i`):

```cpp
eta = (1-relax)*eta + relax*(phiE[e] - phiI[i] - nernst);     // BV 'relax' (combustionProperties)
j   = max(j0*coeff*s^gamma*(exp(n*alpha*F*eta/RT) - exp(-n*(1-alpha)*F*eta/RT)), 0);
SE[e] = -sign(n)*j;   SI[i] = +sign(n)*j;
```

The current is an explicit function of the potentials of the previous outer iteration, damped only through η. Tafel and Constant follow the same pattern.

**Potential solves.** `electricSystem::solveSegregated()` ([electricSystem/electricSystem.C](src/libSrc/fuelCellSystems/electricSystem/electricSystem.C)) solves the electron region and then the ion region, each with J frozen:

```cpp
phiEqn = -laplacian(sigma, phi) - J;
if (phi.needReference())                    // phiI: all-Neumann
{
    iDot = gSum(J*V)/gSum(V);               // net charge source
    for every cell: phiEqn.setReference(id, phi[id] - iDot*ionReferenceRelax, true);
}
phiEqn.solve();
i = -sigma*grad(phi);
dissolved->solve();
```

Calling `setReference(…, true)` on every cell doubles the diagonal. That is an implicit under-relaxation of φI with factor 0.5, and it shifts the φI level by $-\text{ionReferenceRelax}\cdot\overline{J}$. This drives the net charge $\int J_I\,dV$ towards zero over the outer iterations. Charge is only conserved at convergence. If the electron region ever needs a reference (no fixed-value patch), the run stops with an error.

### 3.4 Coupled algorithm

The coupled solve removes three lags of the segregated scheme in one linear system: φE vs φI, kinetics vs potentials, and reference level vs charge balance.

**Butler–Volmer, coupled branch (base point).** The branch is selected by `electricSystem::coupled()`. It requires α = 0.5 and that the model's ion region is the system's `ionRegion`. With

$$
G(\eta) = 2k\,\sinh(f\eta),\qquad k = \max\!\big(j_0\,\text{coeff}\,s^\gamma,\ \text{SMALL}\,j_0\big),\qquad f = \frac{|n|F}{2RT}
$$

the base point is rebuilt from the **realised current** $J_r$ of the previous coupled solve (the ion region's `J`):

```cpp
eta[fluidId]   = asinh(Jr/(2*kin))/fh;            // eta* = G^-1(Jr)
dJdEta[ionId]  = fh*sqrt(sqr(2*kin) + sqr(Jr));   // G'(eta*)
etaE[ionId]    = eta[fluidId] + nernst[fluidId];  // eta* + E
j[fluidId]     = max(sign*Jr, 0);                 // species and heat sources
```

- `dJdEta` and `etaE` are fields owned by `electricSystem` on the ion mesh (not written).
- The BV `relax` is not used. The run prints a one-time note if it is set.
- Tafel and Constant stop with an error in coupled mode.

**Linearisation.** In each CL cell pair (e, i) with cell volume $V$:

$$
J_I \approx J^* + G'(\eta^*)\,\big(\phi_E-\phi_I-(\eta^*+E)\big),\qquad r = \omega\, V\, G'(\eta^*)
$$

$$
\begin{aligned}
\text{electron row } e:&\quad (a_P^E + r)\,\phi_E - \sum a_N^E\phi_{E,N} - r\,\phi_I = -V J^* + r\,(\eta^*+E)\\
\text{ion row } i:&\quad (a_P^I + r)\,\phi_I - \sum a_N^I\phi_{I,N} - r\,\phi_E = +V J^* - r\,(\eta^*+E)
\end{aligned}
$$

The matrix is **symmetric positive definite**: two symmetric Laplacian blocks, a positive semi-definite coupling block $r\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$, and Dirichlet patches on φE. φI needs no reference because the coupling ties it to φE.

**`coupledMatrix`** ([electricSystem/coupledMatrix/coupledMatrix.C](src/libSrc/fuelCellSystems/electricSystem/coupledMatrix/coupledMatrix.C)) holds this system and splits the work into a part done once and a part done every iteration.

*Built once (constructor):*
- the unknown numbering, `[0, nE)` for φE and `[nE, nE+nI)` for φI;
- the cell pairs (ion cell → electron cell with the same parent id), 17 040 in the reference case;
- the face list (φE faces, φI faces shifted by nE, one link per pair) and its upper-triangular order;
- the processor interfaces, whose faceCells are offset like the cells and which reuse the φ processor patch fields for the halo exchange;
- the `lduPrimitiveMesh` with a non-blocking schedule, and the `lduMatrix`.

*Each iteration:*
- `assemble(dJdEta, etaE)` builds `-laplacian(σ,φ) - J` for both regions, folds the non-coupled boundary coefficients into the diagonal and the source b₀, stores the processor-interface coefficients, and adds the $r$ terms;
- `solve(psi, rhs, controls, name)` runs `lduMatrix::solver::New(...)->solve`. It can be called several times on the same matrix, which the implicit control does;
- `scatter(psi)` writes the internal potentials back and corrects the boundary conditions;
- `realise(etaE)` sets the new sources from the solved potentials (below);
- `unitSource(patch)` returns the right-hand side of a unit patch value (§4.3).

**Parallel requirement:** every pair must be on one processor. `decomposeParID` guarantees this by decomposing the parent mesh and carrying the decomposition to the children. Coupled patches other than processor patches are rejected.

**Realised current.** After the solve (and after the implicit control, §4.3):

```cpp
delta = phiE[e] - phiI[i] - etaE[i];        // eta_new - eta*
Jreal = J*[i] + (r/V)*delta;                // tangent current
iR.J[i] = Jreal;   eR.J[e] = -Jreal;
```

These are exactly the sources the linear system balanced, so $\int J_I\,dV$ over both CLs is zero to solver precision:

```
Coupled: pairs = 17040, max|eta_new - eta*| = … V, I+ (A) = …, I- (A) = …, rel. imbalance = …
```

Then $\mathbf i=-\sigma\nabla\phi$ is computed in both regions and the dissolved-water models are solved.

**Why it converges.** The iteration is Newton's method on the *inverse* kinetics. The state carried between iterations is the current, mapped back onto the curve by $\eta^*=G^{-1}(J)$ (asinh, nearly logarithmic in the Tafel branch). A large step in η therefore only produces a linearly extrapolated current, not an exponential overshoot, so full steps (ω = 1) are safe. At the fixed point $J=G(\eta)$ for any ω, so ω damps the iteration without changing the answer. T, the species (through `coeff`), the Nernst potential and σ are still lagged by one outer iteration.

### 3.5 Sequence of one outer iteration (coupled algorithm, implicit control)

```
mapFromCell            T -> air, fuel, phiE, phiI
correct
  air/fuel   correctThermo -> electroChemicalReaction::correct
               Nernst E(T, x, p)
               ButlerVolmer (coupled): Jr -> eta*, dJdEta, etaE, j
               water source / activity -> dissolved model
  phiE, phiI sigma; dissolved correct
  electric   electricControl::beforeSolve: log "ibar: <realised> voltage: <V_old>"
solve
  air/fuel   SIMPLE: Y_i (sources from j) -> U -> p (dmdt)
  electric   coupledMatrix::assemble            M, b0
             solve 1   M phi1 = b0              (patch at V_old)
             afterSolve (implicit control):
               unitSource, solve 2  M phic = c
               ibar1, g, dV;  psi = phi1 + dV phic;  patch = V_old + dV
             scatter; realise -> J_E, J_I; i_E, i_I; dissolved
mapToCell    Qdot (reaction + Joule), rhoCp, k, fluxes
EEqns        T on the parent mesh
```

In the segregated algorithm the electric part is: `beforeSolve` (explicit update) in `correct`, then the φE solve and the φI solve in `solve`.

---

## 4. Operating control (`electricControl`)

[electricSystem/electricControl/electricControl.C](src/libSrc/fuelCellSystems/electricSystem/electricControl/electricControl.C). The controlled quantity is the Dirichlet value V of `patch` in the electron region. The measured quantity is the mean current density through that patch:

$$
\bar i = \frac{\sum_{c\,\in\,\text{zone}} J_{E,c}\,V_c}{A_{\rm patch}}
$$

This sum runs over the electron-region cells of `zone`, or over the whole region if no zone is given. Because phiE holds both electrodes, the zone must be one CL (`ccl`): over the whole region the anode and cathode currents cancel. The controller has two hooks:
- `beforeSolve()`, called from `electricSystem::correct()`;
- `afterSolve(matrix, psi, etaE, controls)`, called inside the coupled solve between solve 1 and the scatter.

Every mode prints `ibar: <ibar from the current J> voltage: <patch value before the update>` once per iteration. That line is what `post_VI.py` reads.

### 4.1 Potentiostatic

`type potentiostatic; voltage <Function1>;` sets the patch to V(t) in `beforeSolve()`. It works with both algorithms.

### 4.2 Explicit galvanostatic

`mode explicit` is the openFuelCell controller, applied in `beforeSolve()` with the J of the previous iteration:

$$
V^{k+1} = V^k + r_E\,(\bar i - \bar i^*)
$$

This is an integral controller with loop gain $r_E\,|\partial\bar i/\partial V|$. With $|\partial\bar i/\partial V|\approx 3.4\cdot10^4$ A/m²/V (the reference case) and $r_E = 10^{-5}$, the gain is ≈ 0.34, which is stable. The remaining error decays by a factor of about $1-0.34$ per iteration, but only while the rest of the system holds still. It works with both algorithms. At a gain close to 2 or above, V oscillates or diverges. That was the `full_thick` failure in openFuelCell2.

### 4.3 Implicit galvanostatic (coupled only)

**Idea.** Within one outer iteration the coupled system is linear, and V enters only through the Dirichlet value of `patch`:

$$
M\,\phi = b_0 + V\,c
$$

where $c$ is zero except in the electron cells next to the patch. So

$$
\phi(V) = \phi_1 + (V - V_{\rm old})\,\phi_c,\qquad M\phi_1 = b_0 + V_{\rm old}\,c,\qquad M\phi_c = c
$$

The realised current is affine in V too, so the V that gives $\bar i = \bar i^*$ follows from one scalar equation. No boundary condition is ever an unknown.

**Implementation, step by step (`afterSolve`):**

1. **`c` comes from the assembly.** `coupledMatrix::unitSource(patch)` copies φE (same patch types), sets the patch value to 1 and assembles `-laplacian(σ, φ, "laplacian(sigma,phi)")`, the same operator as `assemble()`. The patch's `boundaryCoeffs` are accumulated into the face cells. The J source and the explicit non-orthogonal correction stay in b₀, so they are absent from `c`. Nothing about the patch discretisation or its sign is written out by hand.
2. **Solve 2.** $M\phi_c = c$ on the same `lduMatrix`, interfaces and solver controls (`phiCoupled`), logged as `Solving for phiUnit`. It is warm-started from the previous iteration's φ_c, which `electricControl` stores.
3. **Current and sensitivity.** With $s_p = r_p/V_p$, for every pair:
   $$
   J_{E,1} = -\big(J^*_p + s_p\,(\phi_{1,E}-\phi_{1,I}-\eta^*_p-E_p)\big),\qquad
   \frac{dJ_E}{dV} = -s_p\,(\phi_{c,E}-\phi_{c,I})
   $$
   $\bar i_1$ and $g = d\bar i/dV$ are formed with the **same zone sum** as the measured $\bar i$. Expect $g<0$, because raising V lowers the current. $-1/g$ is the fast, frozen-species area-specific resistance (Ω m²).
4. **Step.**
   $$
   \Delta V = \theta\,\frac{\bar i^* - \bar i_1}{g},\quad |\Delta V|\le \texttt{maxStep},\quad V_{\rm old}+\Delta V \in [\texttt{Vmin},\texttt{Vmax}]
   $$
   (ΔV is recomputed if a limit is active), with θ = `relax`.
5. **Superposition.** `psi += ΔV·φ_c` and patch = $V_{\rm old}+\Delta V$. Then `scatter` and `realise` run as usual, so the realised current, the next base point and all downstream code see the combined solution.
6. **Fallback.** If $g \ge -\texttt{gMin}$ the step is skipped and φ₁ is kept. The patch gets the explicit update $V_{\rm old} + r_E(\bar i_1 - \bar i^*)$ instead, and a warning is printed.

The log line is:

```
DICPCG:  Solving for phiUnit, Initial residual = …, Final residual = …, No Iterations …
Implicit control: ibar1 = …, g = … A/m2/V, -1/g = … Ohm m2, dV = …, V = …, ibar = … [(limited)]
```

**Check.** With θ = 1 and no limit active, the predicted `ibar` equals $\bar i^*$. The `ibar:` line of the next iteration, which is the realised current, shows the same value to the printed digits. In the reference run this holds for every unlimited step.

**Convergence criterion.** In implicit mode $\bar i=\bar i^*$ holds by construction, so the current error says nothing about convergence. Judge convergence by |ΔV| per iteration, together with `max|eta_new - eta*|` and the flow and species residuals.

### 4.4 Behaviour: explicit vs implicit at 8000 A/m²

Both runs use the coupled algorithm, 2000 iterations, start at V = 0.4 V with J = 0, and have species relaxation 0.99 (`myFuelCellTest/cpl8000_{explicit,implicit}`, figures in `cpl8000_explicit_vs_implicit/`):

| | explicit (rE = 1e‑5) | implicit (maxStep 0.1 V) |
|---|---|---|
| ī within 1 % from iteration | 17 | 9 |
| ī within 0.1 % / 0.01 % from iteration | 349 / 1134 | 9 / 9 |
| start-up extremes | V ≥ 0.087 V, ī ≤ 23 000 A/m² | V ≥ 0.100 V, ī ≤ 16 700 A/m² |
| V at iteration 500 / 1000 / 2000 | 0.42011 / 0.40595 / 0.40228 | 0.42004 / 0.40591 / 0.40228 |
| PCG iterations per outer iteration | 265 | 145 (solve 1) + 85 (solve 2) |

- **Start-up.** With J = 0 at the start, G′ is tiny and g ≈ −1.4 A/m²/V. The implicit step is clamped at −0.1 V three times and overshoots for two iterations before it settles. The explicit update dives in the same way, since $r_E(0-8000) = -0.08$ V per iteration.
- **After ~50 iterations** the voltage histories coincide. The slow drift that remains (−14, −3, −0.6 mV per 500 iterations) is the species mode caused by the Y relaxation, not the controller.
- **Cost.** Solve 1 becomes cheaper because the superposed field already matches the new patch value. In total the implicit run did ≈13 % fewer PCG iterations.
- **ASR estimate.** At the end, −1/g ≈ 2.9·10⁻⁵ Ω m².

---

## 5. Verification

Every restructuring step was checked against openFuelCell2 (`openFuelCell` built from `7a4480b`) by running the same decomposed case for 50 iterations (`myFuelCellTest/{ref,my,s2}_{seg,cpl}`, `run.sh`, `cmp.sh`).

| step | segregated | coupled | potentiostatic (10 its, both algorithms) |
|---|---|---|---|
| prune, whole files (`610042b`, part 1) | logs and fields identical | identical | – |
| prune, class members (`610042b`) | identical apart from the `dpdt` message and the `Cp` file | identical, same | – |
| electricSystem (`37195d0`) | identical apart from the startup summary | identical, plus the one-time BV `relax` note | identical |
| implicit control (`9b28b91`), explicit mode | identical | identical | – |

"Fields identical" means a byte-for-byte comparison of every field file at t = 50.

---

## 6. Limitations and open items

- **One electron region, one ion region.** `electronRegion`/`ionRegion` are single names, so separate `phiEA`/`phiEC` regions must be merged into one `phiE` (as in the robustness cases). The control zone must then be a single CL.
- **Coupled algorithm:** Butler–Volmer with α = 0.5 only; processor patches only; pairs must be co-located (`decomposeParID`).
- **Implicit control:** coupled algorithm only. A start from J = 0 gives a small |g| and clamped steps for the first few iterations. Possible remedies: a realistic initial J*, or a smaller initial `maxStep`. `lagged` mode (refreshing φ_c only every N iterations) is not implemented.
- **Out-of-date volumetric flux.** `MovingPhaseModel::phi_` (field `phi` of the fluid regions) is computed from U once at construction and never updated. The pressure equation updates the mass flux `alphaRhoPhi`. The written `phi`, the Courant number (`phirMag`) and the flux handed to the turbulence model are therefore stale.
- **Water source sign.** `R(H2O) = fvm::Sp(+w/Y, Y)` weakens the diagonal of the species equation where $Y_{H_2O}$ is small.
- **Forced dissolved model.** `electroChemicalReaction` requires a `dissolvedModel` in the ion region (`dissolveOnOff true` with `type none` for HT-PEM).
- **Phase fraction.** The phase fraction is still a field ≡ 1 that multiplies the Y and U equations and enters Butler–Volmer as $s^\gamma$.
- **Ignored dictionaries.** The two-phase dictionaries in the fluid `regionProperties` are ignored but still required by the case files for compatibility.
- **Case bug (not code).** The air `combustionProperties` uses `cRef $cRef_an` for O₂.

---

## 7. Where to look

| topic | file |
|---|---|
| main loop, energy equation | [src/appSrc/myFuelCell.C](src/appSrc/myFuelCell.C), [EEqns.H](src/appSrc/EEqns.H) |
| parent system, maps | [fuelCellSystem/fuelCellSystem.C](src/libSrc/fuelCellSystems/fuelCellSystem/fuelCellSystem.C), [regions/region/](src/libSrc/fuelCellSystems/regions/region/) (`regionType`, `regionTypeList`) |
| fluid region | [regions/fluid/fluid.C](src/libSrc/fuelCellSystems/regions/fluid/fluid.C) |
| single-phase SIMPLE/PIMPLE | [solvers/singlePhaseSystem/](src/libSrc/fuelCellSystems/solvers/singlePhaseSystem/) (`singlePhaseSystem.C`, `YEqns.H`, `UEqn.H`, `pEqn.H`, `pcEqn.H`, `createFields.H`), [phaseSystem/](src/libSrc/fuelCellSystems/phaseSystem/) |
| phase model chain | [phaseModel/](src/libSrc/fuelCellSystems/phaseModel/) (`MultiComponentPhaseModel::YiEqn`, `dmdt`; `MovingPhaseModel::UEqn`) |
| electrochemical sources | [electroChemicalModel/electroChemicalReaction.C](src/libSrc/fuelCellSystems/electroChemicalModel/electroChemicalReaction.C) |
| kinetics | [activationOverpotentialModels/](src/libSrc/fuelCellSystems/activationOverpotentialModels/) (`ButlerVolmer.C`; heat in `activationOverpotentialModel.C::Qdot`) |
| Nernst | [nernstModels/standard/standard.C](src/libSrc/fuelCellSystems/nernstModels/standard/standard.C) |
| electric regions | [regions/electric/electric.C](src/libSrc/fuelCellSystems/regions/electric/electric.C) |
| electric system, algorithms | [electricSystem/electricSystem.C](src/libSrc/fuelCellSystems/electricSystem/electricSystem.C) |
| coupled matrix | [electricSystem/coupledMatrix/coupledMatrix.C](src/libSrc/fuelCellSystems/electricSystem/coupledMatrix/coupledMatrix.C) |
| operating control | [electricSystem/electricControl/electricControl.C](src/libSrc/fuelCellSystems/electricSystem/electricControl/electricControl.C) |
| conductivity, diffusivity, porosity, dissolved water, saturation | [sigmaModels/](src/libSrc/fuelCellSystems/sigmaModels/), [diffusivityModels/](src/libSrc/fuelCellSystems/diffusivityModels/), [porosityModel/](src/libSrc/fuelCellSystems/porosityModel/), [dissolvedModels/](src/libSrc/fuelCellSystems/dissolvedModels/), [saturationModels/](src/libSrc/fuelCellSystems/saturationModels/) |
| regression and comparison cases | `run/H2/HTPEMFC/myFuelCellTest/` (`run.sh`, `cmp.sh`) |
