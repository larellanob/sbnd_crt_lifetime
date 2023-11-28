# SBND nearline electron lifetime measurement

The goal of this project is to make a quick selection of crossing muon
tracks and make an electron lifetime measurement. It uses ROOT to read
data from the SBND *Commissioning Trees*, and makes selections and
figures.

Quick summary of the macros. Numbered ones are necessary to run in
order to make the electron lifetime measurement:

1. `remake_ct.cxx`:
   - (optional) quick and dirty remake of CRT tracks from CRT hits in
     the commissioning trees.
   - transforms the angles of the TPC tracks to something (which I
     believe is) more consistent
   - matches CRT tracks to TPC tracks using angle matching
   - corrects time (x-position) of TPC tracks, maintaining their
     orientation
   - saves "CORRECTED_" tree

2. `Calculate_elifetime`:
   - Uses "CORRECTED_" tree to select suitable TPC crossing tracks
   - Separates TPC tracks in bins in the x direction
   - For each x-bin makes histogram of dQ/dx per hit
   - Performs fits to the dQ/dx distributions
   - Collects fit parameters to make electron lifetime measurement

Other main macros:

- `Plot_matched_tracks_zxy`
  - Uses "CORRECTED_" tree to make 3D-plots of the matched tracks in
    toy SBND detector geometry
  - 3 different views of the detector, includes LaTeX display of matched angles

- `muonhits.py`
  - python macro to draw event displays (wire,time tick) of muonhits in
    the commissioning trees
  - also draws scaled xz plane using collection plane information
  - has option to draw auxiliary lines in xz plane to draw track
    theta_xz angle
 

