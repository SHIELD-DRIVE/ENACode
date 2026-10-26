#pragma once
// Streamline processing: termination-shock partition fits along streamlines
// (Kornbleuth+2023), streamline extinction (Kornbleuth+2021 Eq. 4), binning to
// the grid and filling of empty cells.

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "ena/config.hpp"
#include "ena/grid.hpp"
#include "ena/snapshot.hpp"
#include "ena/streamlines.hpp"

namespace ena {

// Per-point arrays: the streamlines concatenated, with the
// Tecplot zone-header lines between them present as NaN rows (num rows in
// total). Arrays computed between consecutive points have num-1 entries.
struct StreamPoints {
  std::size_t num = 0;
  // read from the file
  std::vector<double> x, y, z, ux, uy, uz, nh, np, tp, vr, vt, theta, phi, r, region;
  std::vector<double> npu3, tpu3, vrpu3, vtpu3;   // multi-ion: pu3 fluid
  std::vector<double> tssw;                       // multi-ion: SW + pu3 temperature (sound speed)
  // derived
  std::vector<double> v, cs, m2;                  // num
  std::vector<double> ds, tau, ext;               // num-1
  std::vector<double> ntr0, nref0, nen0, etr0, eref0, een0, em0, kappa_ref0, vshock0, esw0;  // num-1
  std::vector<std::int32_t> r_in, t_in, p_in;     // num
  // multi-ion, injection = "streamline": SW speed (cm/s) at the last
  // region-4 point before the streamline enters region 3/6, carried along the
  // rest of that streamline (0 before the crossing or without one).
  std::vector<double> u_inj0;                     // num
  // summary of the termination-shock detection
  std::size_t streamlines = 0, streamlines_with_ts = 0, ts_crossings = 0;
  std::size_t esi_no_inner_point = 0;  // extended single-ion: TS crossings on streamlines with no point inside 12 AU
  std::vector<std::size_t> ts_rows;               // rows where the fits were set (regions method)
  double r_ts_reference = 0;                      // AU, normalization of the fits
};

// Cell fields produced by the streamline processing: every region 3/6 cell inside the
// integration window that is connected to streamline data is filled; other
// cells are 0.
struct StreamFields {
  Field3D<double> cool, vx, vy, vz;
  Field3D<double> ntr, nref, nen, etr, eref, een, em, kappa_ref, vshock, esw;
  Field3D<double> ptr_ct, ptrn_ct;
  Field3D<double> u_inj;          // multi-ion: carried upstream SW speed per cell (cm/s)
  bool has_u_inj = false;
  long unfilled_u_inj = 0;        // region 3/6 cells left without u_inj
  // gap-fill summary: region 3/6 cells inside the integration window
  long cells = 0, unfilled_cool = 0, unfilled_fractions = 0;
  int passes = 0;
};

// Called before each gap-fill pass (pass = 0: the binned state, 1: after the
// first pass, ...) and with pass = -1 for the final state, so checkpoints C3
// can be written.
using StreamPassCallback = std::function<void(int pass, const StreamFields&)>;

// Reference TS distance r_TS,V2 (AU): the configured constant, or the TS
// radius along physics.v2_direction (bilinear over the four nearest lines of
// sight) found with the ts_detection method. `note` receives a description.
double reference_ts_distance(const RunConfig& cfg, const Snapshot& snap, std::string* note = nullptr);

// `region` is needed for ts_detection = "regions"; `r_ref` is r_TS,V2 in AU.
// ts_detection = "regions": how far the flow is from the region-based TS
// (path length from the crossing to the first point with M^2 < 0.64) and the
// v_shock each rule gives. Printed by ena_run.
std::string crossing_diagnostics(const RunConfig& cfg, const StreamPoints& pts);

StreamPoints compute_stream_points(const RunConfig& cfg, const StreamlineSet& sl, double r_ref,
                                   const Field3D<std::int16_t>* region = nullptr);
// The two steps of bin_and_fill: average the points in each cell, then fill
// the empty region 3/6 cells from their neighbours (order-independent passes,
// until nothing changes).
StreamFields bin_stream_points(const RunConfig& cfg, const StreamPoints& pts, const Field3D<std::int16_t>& region);
void fill_stream_fields(const RunConfig& cfg, StreamFields& f, const Field3D<std::int16_t>& region,
                        const StreamPassCallback& on_pass = nullptr);
StreamFields bin_and_fill(const RunConfig& cfg, const StreamPoints& pts, const Field3D<std::int16_t>& region,
                          const StreamPassCallback& on_pass = nullptr);

}  // namespace ena
