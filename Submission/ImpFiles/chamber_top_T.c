/*-----------------------------------------------------------------
  chamber_wall_temperatures.c

  DEFINE_PROFILE UDFs for the TUM SFB/TR-40 single-element GCH4/GOX
  combustion chamber (Test Case 1) wall temperature boundary
  conditions.

  Source data: Test Case 1 report, Annex 3, Table 9 (values
  reconstructed from the in-house tool "Thermtest").

  Two profiles are provided:
    chamber_top_T   -> T_up   (thin walls: top AND bottom)
    chamber_side_T  -> T_side (thick walls: left AND right sides)

  WHY TWO WALL GROUPS, NOT FOUR:
  Internal channel is 12x12 mm, but the block is NOT square outside:
  external width 85 mm, external height 50 mm (Table 1). So:
    - top/bottom wall thickness  = (50-12)/2 = 19 mm   -> use T_up
    - left/right wall thickness  = (85-12)/2 = 36.5 mm -> use T_side
  This matches Figure 7 and the fact that TCLU (top) and TCLD
  (bottom) thermocouple readings track closely together, distinct
  from the side-plane reconstruction.

  *** CHECK Z_OFFSET_M BEFORE USE ***
  Table 9's own z-origin is the START of the first chamber segment,
  29 mm UPSTREAM of the faceplate (see note under Table 9). Other
  tables in the report (feed sensors, Table 7/8) use z=0 AT the
  faceplate instead.
    - If your mesh has z=0 at the FACEPLATE -> Z_OFFSET_M = 0.029
    - If your mesh has z=0 at the start of the first segment
      (i.e. already matches Annex 3's own origin) -> Z_OFFSET_M = 0.0
  Get this wrong and every wall temperature is shifted by 29 mm.

  EXTRAPOLATION NOTE:
  Table 9 only covers z = 1..331 mm in its own frame, i.e. up to
  roughly the throat. Outside that range (e.g. into the nozzle
  segment) this UDF clamps to the nearest tabulated value rather
  than extrapolating -- treat that region's wall temperature as a
  rough approximation, not measured data.

  Compile as a COMPILED UDF (Build, then Load) -- the helper
  function and arrays are not supported in interpreted UDFs.
-------------------------------------------------------------------*/

#include "udf.h"

#define N_POINTS 34
#define Z_OFFSET_M 0.029   /* <-- SET THIS: 0.0 or 0.029, see note above */

/* Table 9: axial coordinate, converted mm -> m */
static const real z_tab[N_POINTS] = {
    0.001, 0.011, 0.021, 0.031, 0.041, 0.051, 0.061, 0.071, 0.081, 0.091,
    0.101, 0.111, 0.121, 0.131, 0.141, 0.151, 0.161, 0.171, 0.181, 0.191,
    0.201, 0.211, 0.221, 0.231, 0.241, 0.251, 0.261, 0.271, 0.281, 0.291,
    0.301, 0.311, 0.321, 0.331
};

/* Table 9: T_up [K] -- top and bottom (thin) walls */
static const real Tup_tab[N_POINTS] = {
    294.9, 296.1, 300.9, 323.4, 337.8, 347.3, 355.3, 362.8, 369.9, 376.9,
    383.6, 390.1, 396.5, 402.7, 408.6, 414.4, 420.0, 425.4, 430.6, 435.6,
    440.5, 445.2, 449.7, 454.0, 458.1, 462.1, 466.0, 469.7, 473.4, 477.6,
    483.3, 494.2, 527.3, 617.2
};

/* Table 9: T_side [K] -- left and right (thick) side walls */
static const real Tside_tab[N_POINTS] = {
    294.9, 296.1, 300.8, 323.3, 337.8, 347.2, 355.2, 362.7, 369.8, 376.8,
    383.5, 390.0, 396.3, 402.5, 408.5, 414.2, 419.8, 425.2, 430.4, 435.4,
    440.3, 444.9, 449.4, 453.8, 457.9, 461.9, 465.7, 469.4, 473.2, 477.2,
    482.8, 492.5, 517.5, 582.7
};

/* Piecewise-linear interpolation, clamped at both ends. */
static real interpolate(real z, const real *zt, const real *Tt)
{
    int k;
    real frac;

    if (z <= zt[0])
        return Tt[0];

    if (z >= zt[N_POINTS - 1])
        return Tt[N_POINTS - 1];

    for (k = 0; k < N_POINTS - 1; k++)
    {
        if (z >= zt[k] && z <= zt[k + 1])
        {
            frac = (z - zt[k]) / (zt[k + 1] - zt[k]);
            return Tt[k] + frac * (Tt[k + 1] - Tt[k]);
        }
    }

    return Tt[N_POINTS - 1]; /* unreachable fallback */
}

DEFINE_PROFILE(chamber_top_T, t, i)
{
    face_t f;
    real x[ND_ND];
    real z_query;

    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        z_query = x[2] + Z_OFFSET_M;   /* x[2] assumed to be the flow axis */
        F_PROFILE(f, t, i) = interpolate(z_query, z_tab, Tup_tab);
    }
    end_f_loop(f, t)
}
