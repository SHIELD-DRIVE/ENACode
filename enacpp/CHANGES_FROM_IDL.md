# Changes from the IDL code (version 2.0)

Version 2.0 of the C++ code has the fixes below built in. They were evaluated one at a time against the IDL results in step 3b and accepted on 2026-10-09; everything except the Q22 shock-shift option.
- To reproduce the IDL results, use the **version 1.0 archive** (`enacpp_v1.0_faithful.zip`). Its faithful mode matches IDL (validate.py PASS 537/537), and every fix below is a switch there.
- The quirk numbers (Qnn) follow `CPP_PORT_PLAN.md`.

**Effect sizes:**
- They are from the step 3b evaluation (full inputs, band 13, unless noted).
- Each is the ratio to the IDL maps for the **total** flux over the 1800 sky pixels: median, 5-95% range, and largest change.
- The solar-wind component lies far out on its Maxwellian tail at IBEX-Hi energies, so its ratio swings by orders of magnitude for small changes in temperature or sampled energy. It is negligible in the band 13 total but matters at IBEX-Lo energies.

## Summary

| ID | Change | Total flux, band 13 | Status |
|---|---|---|---|
| Q29 | ds: r in cm in the angular terms | -7.5% (-16% to -1%; down to -20%) | built in |
| Q01 | ds: phi step from the phi difference | with Q29: -16% (-22% to -7%; down to -32%) | built in |
| Q03 | order-independent gap fill that fills every cell | median unchanged; up to -12% near the pole | built in |
| Q04 | E_sw divided in the gap fill | none with constant E_sw | built in |
| Q10 | survival: E_rel in keV | +0.7% (up to +2.9%) | built in |
| Q11 | get_limits: true radius | none (diagnostics) | built in |
| Q12 | binning to the nearest grid point | +0.3% (-6% to +7%; up to +14%) | built in |
| Q14 | double precision throughout (also Q02) | < 2e-4 | built in |
| Q17 | Moscow region code 0 -> 5 | Moscow only | built in |
| Q18 | sub-band centres | -3.5% (-9% to -2%) | built in |
| Q19 | region 6 flux included | none in this model | built in |
| Q23 | pole fill | replaced by the Q03 fill | built in |
| Q28 | region codes and B in double precision | none | built in |
| Q30 | streamline radial index from r_inner | none while r_lower = r_inner | built in |
| Q31 | integration stops at the first exit from the heliosheath | 2e-5 at one pixel | built in |
| Q32 | accelerated PUIs: explicit cutoff below E_l (no NaN at negative shock-frame speed) | none where IDL runs; IDL stops at IBEX-Lo bands 0-5 | built in |
| **all of the above** | | **-18% (-29% to -11%; down to -40%); accelerated PUIs -33%** | |
| Q21 | fitted E_sw | +0.03% (-0.6% to +0.1%) | option, on by default (`physics.fitted_esw`) |
| Q24 | TS from the regions | +1.1% (0 to +4%); accelerated +7% | option, on by default (`physics.ts_detection = "regions"`) |
| r_TS,V2 | reference TS distance computed along the V2 trajectory | pending (first full-input run) | option, on when a V2 path is set (the reference configuration uses `configs/v2_trajectory_sim.txt`) |
| Q22 | shock-frame shift for the Maxwellians | x2.9 scalar, x2.2 vector | option, **off** (`physics.shock_shift = "none"`) |

## Built-in fixes

### Q29: streamline path length units (`read_stream_g`, line 70) - new
- **IDL:** `ds = sqrt(dr0^2 + r^2 dtheta0^2 + r^2 sin^2(theta) dphi0^2)`. `dr0` is converted to cm (line 65), but r in the two angular terms is still in AU. The angular terms are therefore about 2e26 times too small, and ds is effectively |dr|.
- **Effect of the bug:** in the heliosheath the flow turns tailward and streamlines move largely in angle, so IDL underestimates the path length, the optical depth (tau = n_H/v ds) and therefore the extinction.
- **Fix:** r in cm in the angular terms.
- **Evidence (full inputs, band 13):**

  | Case | Total flux, median (5-95%) | Lowest | Accelerated PUIs |
  |---|---|---|---|
  | Q29 alone | -7.5% (-16% to -1%) | -20% | -17% |
  | Q01 + Q29 | -16% (-22% to -7%) | -32% | -31% |

  - Q01 matters once Q29 is fixed: with Q29 alone, the phi term still uses the theta step (Q01).
  - The largest decreases (about -35% with all fixes) are on the flanks near the ecliptic, where streamlines turn tailward; the smallest are near phi = 0 just north of the ecliptic, where they run nearly radially.
- **Check against the streamline coordinates:** the distance between consecutive points from their x, y, z (slice inputs, subsonic points):
  - With Q01 + Q29, ds / |dxyz| has a median of 1.002 (5-95%: 0.985-1.020), and the total heliosheath path matches to 0.2%. The residual comes from evaluating r and theta at the end of each step.
  - IDL's ds has a median of 0.94 (5-95%: 0.64-1.00), and the total path is 12% short on these polar streamlines.
- Found through the full-input evaluation: Q01 alone gave identical maps.
- **Status:** built in.

### Q01: streamline path length, phi step (`read_stream_g`, line 69)
- **IDL:** `dphi0 = dtheta0*(!pi/180.)`. The phi step is replaced by the theta step, converted to radians twice.
- **Fix:** dphi from the phi difference between consecutive points (wrapped to +/-180 deg).
- **Evidence:** identical maps on its own, because Q29 makes the angular terms negligible. Together with Q29: -16% median (see Q29).
- **Status:** built in.

### Q03: gap filling (`read_stream_g`)
- **IDL:**
  - Ten passes update the arrays in place, so the result depends on the visiting order.
  - The extinction wraps in phi; the fractions only partly.
  - The pole is filled separately.
  - A ring with no streamline data at some radius is never reached, since the extinction only uses neighbours at the same radius.
  - Full inputs: 5 region 3/6 cells left empty (no extinction there).
- **Fix:**
  - Each pass reads the previous pass (order-independent), and every quantity wraps in phi.
  - The pole is one point whose neighbours are the whole of row 1.
  - Where no same-radius neighbour has data, the extinction uses neighbours at r -/+ 1.
  - Passes repeat until nothing changes.
  - All ten fractions are divided by the neighbour count (no Q04).
- **Evidence:**
  - Full inputs: 0 empty cells after 9 passes.
  - Total flux: median unchanged, 5-95% -0.9% to +0.6%, 8% of pixels change by more than 1%; largest -12% at theta 6 deg (next to the pole).
  - Accelerated PUIs: -4% to +4% (5-95%), up to +39%.
  - On IDL's own binned data, cells with streamline data are unchanged; in filled cells the extinction changes by up to 15% at the 90th percentile.
- **Status:** built in.

### Q04: E_sw in the first gap-fill pass
- **IDL:** E_sw is summed over the neighbours but not divided by their number, so filled cells get E_sw x n (up to 26).
- **Fix:** divide.
- **Evidence:** identical maps while the constant E_sw = 0.04 is used. It is required by the Q21 option.
- **Status:** built in.

### Q10: survival probability units (`ionization_rate`)
- **IDL:** E_rel = 0.5 m_p v^2 in erg is passed to `xsection`, which expects keV.
- **Fix:** convert to keV.
- **Evidence:** total +0.7% median (0 to +2.1%, up to +2.9%); accelerated +1.5%.
- **Status:** built in.

### Q11: get_limits radius checks
- **IDL:** uses `ri_opher + k*dr` (30 + 2k) instead of the true radius 2 + 2k.
- **Fix:** true radius.
- **Evidence:** identical maps (r_limits diagnostics only).
- **Status:** built in.

### Q12: streamline binning
- **IDL:** `round((theta - 2.9)/dtheta)` (same for phi), from the cell-centred grid of the streamline tracer, while the ENA grid points are at 0, 6, 12 ... deg.
- **Fix:** nearest grid point. The streamline positions are true coordinates (phi is recomputed from x, y).
- **Evidence:**
  - Total +0.3% median, -6% to +7% (5-95%), up to +14%.
  - The pattern is a band just south of the ecliptic plane, as expected from a half-cell shift of a steep latitude gradient.
- **Status:** built in.

### Q14: precision
- **IDL:** single-precision literals, !pi, energies and float arrays.
- **Fix:** double precision throughout (covers Q02).
- **Evidence:** median 4e-8, largest 2e-4.
- **Status:** built in.

### Q17: Moscow region code 0
- **IDL:** a Moscow region code 0 leaves the cell at 0.
- **Fix:** code 0 becomes 5 (undetermined), so the region-5 check reports it.
- **Evidence:** Moscow inputs only.
- **Status:** built in.

### Q18: band integration
- **IDL:** samples the left edge of each of the 10 sub-bands, while the IBEX-Hi response weights describe the sub-band centres.
- **Fix:** sample the centres.
- **Evidence:** total -3.5% median (-9% to -2%), on every pixel.
- **Status:** built in.

### Q19: region 6 flux
- **IDL:** builds the populations in region 6 (open heliosheath) cells but adds only region 3 flux.
- **Fix:** add region 6 too, as the MI code does.
- **Evidence:** identical maps: this model has no region 6 flux cells in the window.
- **Status:** built in.

### Q23: pole fill subscripts
- **IDL:** `ntr[0,0,rad] = ntr[0,*,rad+1]` writes a row into theta rows 1-29.
- **Fix:** `ntr[0,*,rad] = ntr[0,0,rad+1]`.
- **Evidence:** identical: the case never arises with these inputs. Q03 replaces this code when fixed.
- **Status:** built in.

### Q28: get_limits precision
- **IDL:** region codes and B are read through single-precision scalars.
- **Fix:** double precision.
- **Evidence:** identical maps (B is not used).
- **Status:** built in.

### Q30: streamline radial index (`read_stream_g`, line 157) - new
- **IDL:** `r_in = round((r - r_lower)/dr)` measures the grid index from the start of the integration window. Everywhere else, grid index k means radius r_inner + k*dr (`rlower = (r_lower - r_inner)/dr`, line 161).
- **Effect of the bug:** none in the reference run, where r_lower = r_inner = 2 AU. With, for example, r_lower = 30 AU, every streamline point would be binned 14 cells (28 AU) too far in.
- **Fix:** measure from r_inner.
- **Status:** built in.

### Q31: stop at the first exit from the heliosheath (`make_flux`, line 609) - requested
- **IDL:** a region 3/6 cell is integrated unless the cell just inside it is region 2 or 1. Along a line of sight going 3 -> 2 -> 3 (e.g. a pocket near 600 AU), only the first re-entry cell is skipped; every later region 3/6 cell is integrated again.
  - The survival losses, summed over region 3 cells only, also restart from zero after the gap.
- **Fix:** integration stops at the first cell, after entering the heliosheath, that is beyond the heliopause (region 1/2) or undetermined (region 5, with a warning giving the number of such lines of sight and the first location). Region 6 counts as heliosheath, so 3 <-> 6 does not stop it.
- **Region-4 pockets (version 2.3, user decision 2026-10-13):** a supersonic (region 4) stretch between heliosheath cells, where a wavy termination shock is crossed again, is skipped (no ENAs) and integration continues behind it. `ena_run` warns with the number of pockets, lists the first ten (theta, phi, radii) and writes all of them to `<output>/heliosheath_pockets.txt`, so the region file can be checked. Until version 2.3 these pockets stopped the integration: in the multi-ion model 46 flank lines of sight (theta 66-102, phi 276-318) exit into region 4 at 140-184 AU and lost 54-113 heliosheath cells each (pixels x0.02 of IDL). The single-ion model has none (slice test: identical maps).
- `ena_run` reports, in either mode, how many lines of sight re-enter the heliosheath inside the integration window, how many cells lie beyond the first exit, and the smallest radius at which that happens.
- **Evidence (full inputs, band 13):** one pixel changes, by 2e-5 (theta 72, phi 18): re-entry inside 800 AU is rare in this model.
- **Test:** a planted pocket (region 2 at k = 60-62 on one line of sight; slice inputs):
  - Faithful mode reproduces the IDL rule.
  - With Q31, that line of sight sums only the cells before the exit, and no other pixel changes.
- **Status:** built in.

### Q32: accelerated PUIs at negative shock-frame speed (`calc_diff_flux_zank`, `zank_fn`) - new
- **IDL:** `v_plasma = plasma_frame(...) - vshock` is negative when |v_ENA + u| < v_shock, which happens at IBEX-Lo energies. `(v/vl)^(-del)` of a negative number is NaN, so the flux is NaN and `if j ne j then stop` halts the run. Found when testing all bands on the slice inputs: bands 0-3 stop in IDL (faithful) mode, and bands 0-5 with the region-based TS (which takes v_shock where the flow enters region 3/6).
- **Physics:** the DSA (Zank-Wang) spectrum of the accelerated PUIs has a lower energy limit E_l (`een0` from the TS fits, about 0.79 keV in the shock frame here); below it there are no accelerated PUIs. In IDL this cutoff comes only from the sign of the formula: below v_l it is negative (for delta > 3; delta is 8.3-16.2 here) and is clamped to 0. A negative shock-frame speed is below E_l too, but the formula gives NaN there.
- **Fix:** zank_fn = 0 for a shock-frame speed v <= v_l, applying the E_l cutoff explicitly; this covers negative speeds.
- **Effect:** none where IDL runs (bands 6-17 identical to the bit). Bands 0-5 now run; as expected from E_l, their accelerated component is zero in bands 0-4 and below 4e-5 of the total in band 5.
- **Status:** built in; the E_l interpretation was confirmed by the user.

## Options

### Q21: fitted E_sw (`physics.fitted_esw = true`)
- **IDL:** the fit gives E_sw = 0.04 (V2, flank) to 0.06 (tail), but `make_flux` uses 0.04.
- **Evidence:** total +0.03% median, -0.6% to +0.1%, at most -1.4%; accelerated unchanged.
- **Status:** built in.

### Q22: shock-frame shift for the Maxwellians (`physics.shock_shift`)
- **IDL:** only the accelerated (Zank) PUIs are evaluated in the shock frame.
- **Option:** `"scalar"` or `"vector"` shift for the transmitted, reflected and SW populations too.
- **Evidence:** total x2.9 (scalar) or x2.2 (vector) median, up to x4.3.
- **Status:** built in.

### Q24: TS from the regions (`physics.ts_detection = "regions"`)
- **IDL:** the TS is the first point with M^2 < 0.64 after a supersonic point, at r < 130 AU. Extinction starts where M^2 < 0.81.
- **Option:**
  - The TS is where a streamline enters a region 3/6 cell from a region 4 cell.
  - Extinction starts on leaving region 4.
  - v_shock (`physics.vshock`, version 2.2): `"first_below"` (default) is the flow speed at the first point from the region crossing on with M^2 < 0.64, i.e. IDL's criterion started from the crossing (the 2-cell speed if never reached). `"downstream"` (version 2.1) is the flow speed `vshock_offset_cells` grid cells past the crossing.
  - Until 2026-10-13 v_shock was the speed at the first region 3/6 point itself. That point can still be upstream of the velocity jump (nearest-node binning, MHD shock spread over a few cells). In the MI test this gave the upstream speed in most cells: the accelerated flux vanished at 0.5-0.9 keV and rose x1.5 (median) at 3-8 keV.
  - In the MI model the SW fluid slows at the region boundary (median 2.7 AU past it), but pu3 lags (median 4.9 AU, 90th percentile 51 AU). `ena_run` prints the same diagnostic for the single-ion flow (path to M^2 < 0.64, v_shock under each rule).
  - The fractions and v_shock are carried along the rest of the streamline.
- **Evidence:**
  - The TS is found on the same 1501 of 1582 streamlines as with the Mach test (1513 crossings: 12 streamlines cross twice). On the slices, the radii agree to within one cell.
  - Total +1.1% median (0 to +4%, up to +11%); accelerated +7% median, up to x2.
  - On top of all fixes (Q01-Q31): accelerated x0.90 of faithful instead of x0.67 (+36%); total x0.83 instead of x0.82.
- **Why the accelerated PUIs react:** the fits use r' = ln(r_TS / r_TS,V2) with r_TS,V2 fixed at 78.5 AU. Near V2, a one-cell (2 AU) shift of r_TS changes E_m by about 20% and delta by about 15%, and the accelerated tail at IBEX-Hi energies depends steeply on both. A systematic one-cell offset between the two TS definitions therefore biases every r'.
- **Status:** built in.

### Reference TS distance (`physics.ts_reference = "computed"`)
- **IDL:** r' = ln(r_TS / 78.5 AU) (82 AU Moscow) in the partition fits.
- **Option:** r_TS,V2 from the model itself: the first TS crossing along the Voyager 2 path, with the same method as the TS detection.
  - **Region method:** midway across the 4 -> 3/6 region change, at the nearest line of sight to the path at each grid radius.
  - **Mach method:** the M^2 = 0.64 crossing, interpolated in r; M^2 is bilinear over the four nearest lines of sight.
- **The path:** `physics.v2_trajectory_file`, with columns r, theta, phi. The trajectory is the spacecraft's, valid for any model on this grid, including the multi-ion model (user, 2026-10-13). Beyond the file's last row (97 AU), theta and phi are extrapolated linearly in r from a least-squares line through the last 5 AU of the file (0.156 deg/AU in theta, 0.063 deg/AU in phi; user decision 2026-10-13). The multi-ion TS along V2 (111 AU) lies there: theta 121.9, phi 218.1 on the extrapolated path, on the same nearest line of sight as before. `configs/v2_trajectory_sim.txt` is V2's trajectory in the BU model's simulation frame, digitized from the user's plot (36-97 AU, about +/-0.1 deg; theta 117.5, phi 216.3 at 84 AU). A fixed `physics.v2_direction = [theta, phi]` is the alternative.
- **Test (slice inputs, a test path through the covered rows):** the region method gave 93.0 AU, identical to an independent check from the region array; the Mach method gave 93.9 AU.
- **Full inputs (2026-10-09):** r_TS,V2 = 69.0 AU with the region method (theta 114.1, phi 214.9 on the path); the Mach method gives 69.93 AU (theta 114.3, phi 215.0), within half a grid cell. The two runs still differ in the accelerated PUIs (x0.74 Mach/regions, median) because each streamline's TS is placed differently, not because of r_TS,V2. Band 13 total vs. IDL: 0.82 median (all fixes and options).
- **Decided (user, 2026-10-13):** r' = ln(r_TS / r_TS,V2) with the paper's cap at ln(1.62), unchanged; r_TS,V2 and every TS from the regions, so the whole code uses one TS definition.
- **Status:** on in the reference configuration (`ts_detection = "regions"`, V2 trajectory file).

## Open points
- **81 of 1582 streamlines have no TS** under either method, so they carry no fractions and are left out of the binning. It is worth checking what they are, e.g. streamlines that end before the TS.
- **The V2 direction** in the simulation frame is needed for the computed r_TS,V2.

## Not switches
- **Q20:** streamlines are always processed. `extinction = false` only sets ext = 1, where IDL would stop on zero densities.
- **Q02:** covered by Q14.
- **Q05, Q13, Q25:** concern diagnostic output only.
- **Q06:** the IDL constants remain the defaults; Q24 and the computed r_TS,V2 replace the TS-related ones.
- **Q07:** radial_range 0/1 report "not implemented".
- **Q16, Q26:** intended behaviour.
- **Q27:** milestone 4.
- **Q08, Q09, Q15:** resolved or harmless.

## Other changes in version 2.0
- **Precision:** besides Q14, the survival-exponent array (a float array in IDL) and the binned flow velocities are double, and degrees per radian is 180/pi instead of 57.2958.
- **boundaries.txt:** rows are placed by their i, j, k columns in any order. IDL's get_limits relied on theta -> r -> phi order, which version 1.0 enforced.
- **Configuration:** `run.mode` and `[quirks]` are gone (an error names the version 1.0 archive); `[legacy] ri_opher` is ignored.
- **Output:** `flux_bandNN.h5` records the code version, the settings and the geometric factor. The legacy count-rate `.dat` files are unchanged in format.
- **Checkpoint C4:** the Kornbleuth+2021 alpha(r_TS) diagnostics (`np_ts`, `rad_ts`, `nh_avg`, `rts`, `a`) and the `[15,33]`-referenced ratios `np_rat`, `tp_rat` (Q05) are no longer written; they were not used by the model.
- **Warnings:** `ena_run` reports region-5 cells, lines of sight that leave and re-enter the heliosheath or end at a region-5 cell, gap-fill cells left empty, and a constant r_TS,V2 combined with the region-based TS.

## Multi-ion model (version 2.3)

`model.populations = "multi_ion_inca"` (configuration `configs/run_mi.toml`) is the multi-ion code `make_countrate_files_inca.pro` (Kornbleuth+2026) with the decisions below. The IDL behaviour is reproduced by the version 1.1 archive (`enacpp_v1.1_mi.zip`). It was validated against the hooked IDL run on 2026-10-13: ENA maps agree to < 1e-6, and validate.py passes 710 of 711 arrays.

**Shared with the single-ion version 2:** double precision; path length and grid indices (Q01, Q29, Q12, Q30); the gap fill (Q03, Q04); stop at the first heliosheath exit or a region-5 cell (Q31); E_rel in keV in ionization_rate (Q10); sub-band centres (Q18); boundaries read in double precision (Q28).

**Multi-ion decisions (user, 2026-10-12/13):**

| Item | IDL | Version 2.3 | Effect (full MI inputs, bands 9/13/14, median) |
|---|---|---|---|
| Q33 `read_stream_inca` | q = 7.65 in tail belt cells (left over from a retired population) | dropped: q = 5 everywhere; `read_stream_mi` dropped (unused) | IDL test: band 14 total up to +8.7% (145 tail pixels) |
| Q34 injection energy | E_l from the SW speed at 10 AU on the cell's own line of sight | kept as the default (`multi_ion.injection = "line_of_sight"`; E_l about 1-2 keV, consistent with Kornbleuth+2023, about 0.8-1.2 keV; user decision 2026-10-13). Option `"streamline"`: SW speed just upstream of the region-based TS, carried along the streamline (E_l about 0.5-0.8 keV) | streamline vs line of sight: accelerated x379/0.49/0.44, total x1.016/0.83/0.57 |
| Q35 streamline file | `read_ascii(data_start=32)` on a 31-line header drops the first point | all points read | none on the maps |
| TS fits | every supersonic point (pu3 >= 215 km/s, M^2 > 4); x = ((r - 106)·38/90 + 106)/106 | at the 4 -> 3/6 region crossing; r_TS,V2 from the regions along V2 (111 AU); `ts_scaling = "compressed"`: x = ((r - r_TS,V2)·38/90 + r_TS,V2)/r_TS,V2, tailored to the hybrid fits; `"ratio"` (x = r/r_TS,V2) as an option | with v_shock "first_below": total x0.999/1.008/1.021; reflected x1.00/1.07/1.16; accelerated 0 at band 9 (Q36), x1.03 at bands 13/14 |
| Polar blend | toward fixed polar values beyond 45 deg; at the TS point effectively only at theta <= 25 or >= 155 deg | applied at theta <= 25 or >= 155 deg | included above |
| Q36 v_shock | pu3 speed at the first point below 215 km/s; 0 before it (8.6% of the heliosheath cells, 79-90% of the accelerated flux at 0.5-0.9 keV) | `physics.vshock = "first_below"`: the same criterion started from the region crossing, applied from the crossing on | accelerated flux at 0.5-0.9 keV removed (artifact) |
| Fits carried | also across the header rows into the next streamline | along the same streamline only | included above |
| kappa of the pu2 slot | carried from the previous cell into region 6 (8.6 after a belt cell) | reset in every cell (kappa_oh; 8.6 in belt cells) | not measured separately |
| E_SW (Q21) | constant 0.04 | fitted value (`physics.fitted_esw`), as the single-ion model (accepted by the user, 2026-10-13) | SW component only; totals identical |
| Ne4, Neu maps | always 0 | not written | none |
| Reconnection region | file coordinates, r <= 550 AU, phi < 30 or > 330 deg; at the pole the file gives theta 0.57, phi 0 for every phi index, so the whole polar row qualifies | grid coordinates: the pole (theta 0, x = 0) is excluded (user decision 2026-10-13); limits are settings (`[multi_ion] belt_*`) | polar row only: total x0.29-0.49 at the pole pixel |

The MI `get_limits` relabelling (a region-3 cell directly beyond a region-1/2 cell takes that region) is kept: it is part of how the MI regions are defined.

## Extended single-ion model (version 2.4)

`model.populations = "extended_si"` (configuration `configs/run_esi.toml`) is the extended single-ion code: `make_countrate_files.pro` with `nhpops = 1`, plus the MI versions of `ena_prod_trans.pro` and `distribution_functions.pro`. The decisions below are applied. The IDL behaviour is reproduced by the version 1.2 archive (`enacpp_v1.2_esi.zip`). It was validated against the hooked IDL run on 2026-10-14:
- validate.py passes all 737 arrays;
- the ENA maps agree to < 1e-6;
- the physics functions agree to 3e-13.

**What the model adds to the single-ion version:**
- **Neutral populations.** Four neutral populations, from the 41-column neutral file.
- **Cross section.** It is evaluated at the parent ion's thermally averaged speed relative to the neutral bulk flow (vrel), for every population except the accelerated PUIs.
- **Heliosheath-created PUIs (`heliosheath_pui`, IDL PUI4).**
  - Density: n_p − n1 − n2 − n3, which includes the SW share (Q26, intended).
  - Temperature: ½ m_p |u_p − u_H|²/k.
  - kappa: kappa_oh.
- **Neutral components (`ne4`, `neu`).** Populations 1 and 2 at the inner boundary, with density × 0.639307. Like IDL, they are part of the total (CombPUIs).
- **Injection energy.** The DSA injection energy E_l is ½ m_p v², with v the flow speed at the last streamline point inside 12 AU. It replaces the E_l fit; on the user's inputs it is 1.0–2.2 keV.

**Shared with the single-ion version 2:**
- double precision;
- path length and grid indices (Q01, Q29, Q12, Q30);
- the complete gap fill (Q03, Q04);
- stop at the first heliosheath exit (Q31; region-4 pockets skipped with a warning);
- E_rel in keV in ionization_rate (Q10);
- sub-band centres (Q18);
- boundaries read in double precision (Q28);
- zank_fn = 0 at v <= 0 (Q32);
- TS from the regions, with r_TS,V2 computed along the V2 trajectory (Q24);
- v_shock "first_below";
- fitted E_SW (Q21).

**Extended single-ion decisions (user, 2026-10-14).** The effects below were measured with version 1.2 on the user's full inputs, each change on its own against the faithful run.

| Item | IDL | Version 2.4 | Effect |
|---|---|---|---|
| Q37 region-6 test lines | n1 x 0.8 with 0.2 n1 in the SW slot, T_SW = (1 + (B/0.15)²) T1, no heliosheath PUIs, everything 1e-10 in theta rows 0–1 (left from the reconnection tests) | removed: region 6 gets the region-3 populations, the SW with kappa_oh | measured with Q03 in both runs: total ×1.08 median at band 3 (to ×1.37), ×1.01 at bands 6/9, ×0.97 at band 13 (to ×0.63), ×0.93 at band 14 (to ×0.41). The hot test SW of region 6 disappears; region 6 is on 1679 of 1800 lines of sight |
| Q38 region-2 losses | ionization_rate also accumulates losses in region 2, with T × 7 (a physics test) | region 3 only, as the SI and MI codes | 8 re-entering lines of sight: up to +7.5% at band 3, ≤ 0.5% at bands 6–14 |
| Q27 neutral velocity units | populations 1–2 scaled × 1e2 (m/s), total H × 1e5 | × 1e5 for all: the file is in km/s for every population (user's neutral slice: v_r ≈ −2.4 to −3.1, v_t 19–25 km/s at 2 AU) | Ne4/Neu only, ≤ 1e-40 of the other components at bands 3 and above |
| `inca_belt = 1` path | belt test (κ 1.9, 1090 km/s shift, T × 6.5) | not ported (user) | none (off in IDL) |
| `radial_range` 1/2 block | flux from the ISM and OHS neutrals into slots 0/1 | not ported | none (radial_range = 0) |
| Gap fill | 10 passes left 120 heliosheath cells without fractions (hidden by the Q37 test lines in polar region 6) | complete (Q03): no cell left empty | total within ±7% (band 3 down to ×0.82), median ×1 |
| Q31 | integration continues after the first heliosheath exit | stops there | 8 re-entering lines of sight: band 3 down to ×0.73, ≤ 4% at the other bands |
| v10 for E_l | kept from the previous streamline if a streamline has no point inside 12 AU | taken on the same streamline (first point if none inside 12 AU; counted in the log) | none on the user's inputs (every streamline starts at 2 AU) |

Also: `get_limits` is used without the MI relabelling, as in the extended SI IDL code.


## Version 2.5: survival in region 6 (Q39) and region-6 ions

**Q39, survival probability in region 6 (all three models; user, 2026-10-14).**
- **IDL:** `ionization_rate` accumulates the loss exponent only in region-3 cells, and each cell adds to the cell just before it. Region-6 cells therefore got S = 1, and a region-3 cell after region 6 started again from 0. This was written before region 6 produced ENAs.
- **Now:** losses are counted in regions 3 and 6 from `survival_observer_r` (100 AU) outward. The accumulated exponent is carried along the line of sight through every other cell, so a region-6 cell gets the losses of all the heliosheath in front of it.
- **Where the maps change:** only on lines of sight that cross region 6 inside the integration window. The single-ion slice test and the synthetic multi-ion run are unchanged. Version 1.3 has the switch `quirks.Q39`, to measure the effect on the full inputs.

**Region-6 ions (extended single-ion model; user, 2026-10-14).** IDL used the same ion populations in regions 3 and 6, apart from the SW kappa (kappa_oh). Region-6 ions should generally be more accelerated. `[region6]` sets:
- **The kappas** of the transmitted PUIs, reflected PUIs, heliosheath PUIs and SW. At fixed temperature a kappa changes the shape of the distribution, not the pressure.
- **An accelerated fraction of those ions.**
  - Inputs: a fraction f of their density takes a fraction e of their pressure, with `accelerated_kappa`.
  - Remaining ions: (1 − f) n_i at T_i (1 − e)/(1 − f).
  - `accelerated_from = "each_population"` (default): every population keeps its own accelerated part, f n_i at (e/f) T_i.
  - `accelerated_from = "pooled"`: one population at (e/f) times the density-weighted mean temperature.
  - Density and pressure are conserved cell by cell (checked to 1e-15).
  - Energy conservation means that giving the tail more energy cools the thermal part. On a synthetic case with f = 0.1, e = 0.3 and kappa 1.8, band 14 rose ×1.8–2.5 on the region-6 pixels while band 13 fell to ×0.8.
- **What is not changed:** the DSA-accelerated PUIs keep their settings. There is no magnetic-field dependence; the B-dependent heating was part of the removed reconnection test.
- **`physics.kappa_heliosheath_pui`** sets the region-3 heliosheath PUIs separately. The default is `kappa_oh`, as in IDL, where `kappaIHS = kappaOH`.

All defaults reproduce the version 2.4 populations.


## Version 2.6: threads (milestone 5)

There are no changes to the physics. The maps and dumps are bit-identical to version 2.5 for any thread count. This was checked on the single-ion slice test and on the synthetic multi-ion and extended single-ion runs (all checkpoint arrays), and on the full-size synthetic run (all maps and legacy files).

## Version 2.7: time-dependent maps (milestone 8, stage A)

This version adds a new capability that the IDL code does not have: maps per energy band and per output time window, made from a sequence of snapshots (`[time_dependent]`, README "Time-dependent maps"). Steady-state runs are unchanged: their maps are bit-identical to version 2.6 for the single-ion slice test and the synthetic multi-ion and extended single-ion runs.

Choices, as agreed with the user on 2026-10-14:
- **Time cells.** Snapshot s holds for [T0 + s·dt, T0 + (s+1)·dt), with no interpolation in time. This is a decision, not a simplification to be replaced (user, 2026-10-16): interpolating between snapshots could smooth or distort time-dependent features such as instabilities. Interpolation may be revisited later.
- **Travel.** Straight radial paths at the ENA speed; the observer is at 1 AU.
- **Survival along the ENA's path.** Each cell's loss term comes from the snapshot covering the time the ENA passes it. The emission time is the centre of the time cell, and the last snapshot is used after it.
- **Arrival.** Arrivals are binned by overlap into the output windows (window averages).
- **Coverage** (user, 2026-10-15). A line of sight is a region with missing flux where ENAs made before the first snapshot are missing and the extinction of the missing heliosheath cells is ≥ 0.3 (`coverage_extinction`, about the cooling length).
  - An optional second test flags it when the missing share of the window's band flux reaches `coverage_missing_flux` (off by default). The share includes every component, so the heliosheath-created PUIs (not reduced by the extinction) count with their contribution at each energy.
  - Flagged lines are set to 0 by default, or kept with `include_incomplete_lines = true`, and are always flagged in the output, with a clear warning.

Checks:
- **Steady-state limit.** Repeating one snapshot reproduces the steady-state maps in every complete window to 2e-15: single-ion (slice test), multi-ion and extended single-ion (synthetic), all components.
- **Step change.** A density step in the middle of the sequence is checked against an independent Python recomputation (cell-by-cell path sums and explicit overlaps) to 1e-15.
- **Thread count and energy grouping.** 1 and 4 threads, and 12 energy groups, give bit-identical results.
- **Missing share.** With identical snapshots, the estimated missing share equals 1 − (time-dependent map / steady map) to 7e-16 in every band, window and line of sight.

## Version 2.8: MPI and HDF5 snapshots for every model (milestone 6)

No changes to the physics.
- **MPI.** Time-dependent runs can use several MPI ranks: production over snapshots, an exchange by line of sight, and integration over lines of sight. The maps are bit-identical to version 2.7 for any number of ranks and threads. This was checked with 1, 2, 3, 4, 5 (with energy groups) and 7 ranks on the extended single-ion step-change and steady-limit tests, and with 6 ranks on the multi-ion test.
- **HDF5 snapshots for every model.** `ena_convert` stores all ion fluids and neutral populations (format 2), so the multi-ion and extended single-ion models can read HDF5 snapshots. The maps are identical to those from the text inputs. Version-1 single-ion files are still read.
- **Converter options.** `ena_convert` gains `--first/--last`, `--skip-existing` and `--compress`, and splits its snapshots over MPI ranks.
- **Steady-state runs** are unchanged (bit-identical).
- **HDF5 file locking** is switched off at program start (unless `HDF5_USE_FILE_LOCKING` is set). On the user's NFS workstation it made every open re-read the snapshot over the network (4.3 s instead of 0.3 s).
- **Master configuration `configs/all_settings.toml`:** every setting, with defaults and IDL values. To make one file serve every model and both run modes:
  - `[region6]`, `[multi_ion]` and `physics.kappa_heliosheath_pui` are skipped, with a note, when another model is selected (before, this was an error);
  - `time_dependent.enabled = false` keeps the section in the file without using it;
  - with `[time_dependent.files]`, `[[snapshot]]` entries may also be present: a time-dependent run uses the files, a steady-state run the `[[snapshot]]` entries.

## Milestone 7: maps in ecliptic J2000 (`tools/ena_plot.py`)

This replaces `plotter_eps.pro`.
- **Rotation to ecliptic J2000:** the maps are rotated exactly, fixed by the solar rotation axis (J2000) and the interstellar flow known in both frames (user, 2026-10-15). IDL only shifted the map in longitude: `round(281/dphi)` columns, about +78° in the data, with the latitude correction computed but not applied.
- **Effect:** relative to the IDL figures, latitudes on the flanks move by up to ±7.2°, longitudes by about 4°, and the poles by 7.3°.
- **Display:** cell-centred cells, Mollweide projection, and boxes over the regions with missing flux.
- **`--idl`** reproduces the old placement for comparison.
- **IDL plot settings kept** (in `configs/plot_settings.toml`): centring on the nose, the tail or ecliptic (0, 0), with the tilted centre of IDL `map_set`; flux scaling (`scale`, IDL `plot_scl`, shown as "(Flux x S)" in the title); colour-bar limits for all maps or per band (IDL `manual_minset`/`manual_maxset`); linear scale by default; UPWIND/DOWNWIND/V1/V2 markers; −180..180 labels and edge meridians.
- **Colour bar** (user, 2026-10-15): default limits per band from the IDL `manual_maxset` for bands 9-17 (265, 140, 60, 30, 12, 0.90, 0.08, 0.025, 0.009), minima 0. Numbers as integers, floats to the nearest hundredth, or exponents with the mantissa to the nearest hundredth (`legend_style`). A range taken from the data is rounded outwards so the labels are exact; IDL printed the unrounded tick values.
- **Not carried over** (user, 2026-10-15): the remaining `plot_flags` options.

