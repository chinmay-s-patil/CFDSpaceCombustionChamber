
# Fluent → OpenFOAM Migration Guide

## TC1 Single-Element GOX/GCH4 Combustion Chamber

Solver: **reactingFoam** (compressible, multi-species, PIMPLE pressure-based reacting flow).
This is the right call — reactingFoam natively has kOmegaSST, native least-squares gradients,
and a CHEMKIN reader, so a surprising amount of your setup ports over almost directly.

---

## 0. Before you do anything: try to recover the Fluent case

Check your working directory for autosaved files like `FFF-*.cas.h5` / `FFF-*.dat.h5`
(you had autosave writing these every 10 iterations per your solver log). If one exists,
you may still have your mesh + most BCs intact — worth a quick look before rebuilding
from scratch in either tool.

---

## 1. Mesh: Fluent → OpenFOAM

**fluentMeshToFoam only reads legacy ASCII/binary `.msh` files, not `.cas.h5` (CFF format).**

1. In Fluent (or Fluent Meshing), export the mesh only:
   `File → Export → Mesh...` → save as Fluent Format `.msh` (not `.cas.h5`).
2. Convert:
   ```
   fluentMeshToFoam case.msh
   ```
3. Fix patch types — Fluent zone names carry over as patch names, but everything
   imports as generic `patch`. Edit `constant/polyMesh/boundary`:
   - `gox_inlet`, `gch4_inlet` → type `patch`
   - `outlet` → type `patch`
   - `chamber_top_wall`, `chamber_side_wall`, etc. → type `wall`
4. Sanity check:
   ```
   checkMesh -allTopology -allGeometry
   renumberMesh -overwrite
   ```

   Fluent poly-hexcore meshes with your BOI/inflation regions convert fine topologically,
   but check non-orthogonality/skewness numbers — thresholds Fluent tolerates aren't
   always what OpenFOAM's linear solvers like.

---

## 2. Chemistry: CHEMKIN → OpenFOAM

You have `gri16_reduced.inp` (16 species, 60 reactions, `!THERMO` commented out — no
inline thermo data, and no transport file). You'll need companion files:

- **Missing THERMO block** → grab standard **GRI-Mech 3.0 `therm.dat`** — all 16 of
  your species (H2, H, O, O2, OH, H2O, HO2, H2O2, CH4, CH3, CO, CO2, HCO, CH2O, N2, AR)
  exist in it, since your mechanism is a GRI-derived reduced set.
- **No transport file** → grab GRI-Mech 3.0's companion `tran.dat` (Lennard-Jones
  transport params) — same reasoning, all species are covered.

Two ways to bring it in — either works, pick one:

**A) Native CHEMKIN reader (least manual work)** — in `constant/thermophysicalProperties`:

```
chemistryReader chemkinReader;
foamChemistryFile        "$FOAM_CASE/constant/reactions";
foamChemistryThermoFile  "$FOAM_CASE/constant/thermo";
CHEMKINFile          "$FOAM_CASE/chemkin/gri16_reduced.inp";
CHEMKINThermoFile    "$FOAM_CASE/chemkin/therm.dat";
CHEMKINTransportFile "$FOAM_CASE/chemkin/tran.dat";
newFormat yes;
```

**B) Pre-convert with `chemkinToFoam` (lets you inspect/edit the parsed mechanism)**:

```
chemkinToFoam gri16_reduced.inp therm.dat tran.dat constant/reactions constant/thermo
```

**Check the log output confirms 60 reactions parsed** — the same number I counted
directly from your `.inp` earlier. If it reports something else (e.g. dropping a
DUPLICATE pair), that's a parser mismatch worth chasing down before running anything.

`thermoType` block (standard for reactingFoam + CHEMKIN mixtures):

```
thermoType
{
    type            hePsiThermo;
    mixture         reactingMixture;
    transport       sutherland;
    thermo          janaf;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}
```

---

## 3. Combustion + turbulence models

| Fluent                            | OpenFOAM equivalent                                                                                                                                                                                                                                                                                |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Finite-Rate/Eddy-Dissipation      | `combustionModel PaSR;` (closest 1:1 analog) — or step up to `EDC` once stable, which is a more rigorous version of what FR/ED approximates                                                                                                                                                   |
| Finite-rate chemistry integration | `constant/chemistryProperties`: `chemistrySolver ode; odeCoeffs { solver seulex; }` — stiff ODE integrator, handles your 16-species mechanism's stiffness directly                                                                                                                            |
| SST k-ω, Low-Re                  | **No porting needed.** `constant/turbulenceProperties`: `simulationType RAS;` + `RASModel kOmegaSST;` in `constant/momentumTransport`. Set `lowReCorrection yes;` if using an OpenFOAM version with that switch — same caveat as before: verify actual y+ on your (coarser) mesh. |

---

## 4. Boundary conditions

| Fluent BC                                       | OpenFOAM equivalent                                                                                                                                                                                                                                                                                                                                                      |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `gox_inlet` mass-flow inlet (T=278K)          | `U`: `flowRateInletVelocity` (massFlowRate); `T`: `fixedValue 278`; `O2`: `fixedValue 1`, all other Yi `fixedValue 0`; `p`: `zeroGradient`                                                                                                                                                                                                             |
| `gch4_inlet` mass-flow inlet (T=269K)         | Same pattern:`CH4`: `fixedValue 1`, rest 0; `T`: `fixedValue 269`                                                                                                                                                                                                                                                                                                |
| Outlet, gauge P=1.0 bar                         | **Check where your outlet patch actually sits relative to the throat.** If it's downstream of where flow goes sonic/supersonic (per your own earlier note — negligible backpressure influence), a hard `fixedValue` p is wrong; use `waveTransmissive` or `inletOutlet`/`zeroGradient` so it doesn't force a backpressure the flow can't physically feel. |
| Walls,`chamber_top_T`/`chamber_side_T` UDFs | Port near-verbatim as a`codedFixedValue` BC in `0/T` — see below                                                                                                                                                                                                                                                                                                    |
| Walls, velocity                                 | `noSlip`                                                                                                                                                                                                                                                                                                                                                               |

**Porting your wall-temperature UDFs** — this is the easiest part of the whole migration,
your interpolation logic drops in almost unchanged:

```cpp
// 0/T, boundaryField for e.g. chamber_top_wall
chamber_top_wall
{
    type            codedFixedValue;
    value           uniform 300;
    name            chamberTopT;
    code
    #{
        const vectorField& Cf = patch().Cf();
        scalarField& field = *this;

        const scalar zOffset = 0.029; // same Z_OFFSET_M as your UDF

        static const scalar zTab[34] = {0.001,0.011, /* ...same as chamber_top_T.c... */ 0.331};
        static const scalar TupTab[34] = {294.9,296.1, /* ... */ 617.2};

        forAll(Cf, i)
        {
            scalar z = Cf[i].z() + zOffset;
            // same clamped piecewise-linear interpolate() logic as your .c file
            field[i] = interpolateTable(z, zTab, TupTab);
        }
    #};
}
```

Just paste your existing `z_tab`/`Tup_tab`/`Tside_tab` arrays in and reuse the same
clamped-linear-interpolation logic — it's a copy-paste port, not a rewrite.

---

## 5. Discretization — your Fluent choices map directly

| Fluent                     | `system/fvSchemes`                                                                                                                                       |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| First-order upwind (all)   | `div(phi,U) Gauss upwind; div(phi,Yi_h) Gauss upwind; div(phi,K) Gauss upwind;` etc.                                                                     |
| Least-squares gradient     | `grad(U) leastSquares; default leastSquares;` — exact match, no approximation needed                                                                    |
| Implicit, global time step | `ddtSchemes { default Euler; }` in `system/controlDict`: `adjustTimeStep yes; maxCo 0.3;` (start conservative given the stiffness you were fighting) |

Note: reactingFoam has no literal "AUSM" flux option — it's PIMPLE (pressure-based), not
density-based, so that setting doesn't carry over. Not a loss: PIMPLE's implicit
pressure-velocity coupling tends to be considerably more forgiving of the mesh-size
disparity that was driving your Fluent divergence in the first place.

---

## 6. Initial patches — your cell registers become `topoSet` + `setFields`

Direct equivalents of what you built in Fluent:

```
// system/topoSetDict
actions
(
    {
        name    ignitionKernel;
        type    cellSet;
        action  new;
        source  sphereToCell;
        origin  (0 0 0.005);
        radius  0.005;
    }
    {
        name    goxCore;
        type    cellSet;
        action  new;
        source  cylinderToCell;
        p1      (0 0 -0.05); p2 (0 0 0); radius 0.002;
    }
    {
        name    gch4Annulus;
        type    cellSet;
        action  new;
        source  cylinderToCell;
        p1      (0 0 -0.05); p2 (0 0 0); radius 0.003;
    }
    // then intersect gch4Annulus with NOT goxCore if needed for the true annulus
);
```

Then `system/setFieldsDict` patches T and every Yi in each region exactly like your
Fluent registers — same composition table we built for the ignition sphere (H2O 0.35,
CO2 0.15, CO 0.30, H2 0.15, OH 0.05, CH4/O2/N2/AR = 0), and pure O2=1 / CH4=1 in the
core/annulus registers respectively. Run in this order:

```
topoSet
setFields
```

---

## 7. Execution order, start to finish

```
fluentMeshToFoam case.msh
checkMesh -allTopology -allGeometry
renumberMesh -overwrite

chemkinToFoam gri16_reduced.inp therm.dat tran.dat constant/reactions constant/thermo
# (or configure native chemkinReader instead — see §2)

# build 0/, constant/, system/ per tables above

topoSet
setFields

reactingFoam
```

Watch `Co` (Courant number) and species boundedness the first several hundred timesteps
the same way you were watching T-limiter cell counts in Fluent — if PaSR/EDC + your
patched kernel doesn't ignite within the first few hundred iterations, the same root
cause list applies (mixing resolution, patch normalization, mesh near the shear layer).
