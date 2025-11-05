Neural dynamics for working memory and evidence integration during olfactory navigation in Drosophila
Code and simulation materials for Kathman et al., 2025 (bioRxv)

--------------------------------------------------------------------
OVERVIEW
--------------------------------------------------------------------
This repository provides the MATLAB code used to simulate odor-guided navigation
in a turbulent plume environment. The model includes neural dynamics for three
behavioral states (baseline, goal-directed, and search), incorporating adaptive
odor thresholds and probabilistic state transitions.

The code replicates example simulations described in the manuscript:
Kathman, N., Lanz, A., Freed, J., Nagel, K. (2025). "Neural dynamics for working memory and evidence integration
during olfactory navigation in Drosophila." Nature Communications (in review).

--------------------------------------------------------------------
CONTENTS
--------------------------------------------------------------------
File / Folder                Description
--------------------------------------------------------------------
run_nav_simulation.m          MATLAB function to run a navigation simulation
LICENSE.txt                   Open-source license (GNU)
README.txt                    This documentation file

--------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------
 - MATLAB R2021a or later
 - No special toolboxes required
 - Local copy of the plume dataset (10302017_10cms_bounded_2.h5)

--------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------
1. Download this repository as a ZIP file or clone using:
   git clone https://github.com/nagellab/Kathmanetal2025.git

2. Open MATLAB and set the repository folder as your working directory:
   cd path/to/Kathmanetal2025

--------------------------------------------------------------------
DATA SOURCE AND PATH SETUP
--------------------------------------------------------------------
This simulation uses a publicly available plume dataset:

   From Álvarez-Salvado et al. 2019 (Plume dataset originally from Crimaldi lab)
   https://dx.doi.org/10.5061/dryad.g27mq71

1. Download the file:
   10302017_10cms_bounded_2.h5

2. Place it in a convenient directory, for example:
   /Users/yourname/Data/PlumeData/10302017_10cms_bounded_2.h5

3. In the code (run_nav_simulation.m), update the line near the top:
      plume_path = '/path/to/10302017_10cms_bounded_2.h5';
   Replace it with the full path to the .h5 file on your system.

--------------------------------------------------------------------
HOW TO RUN
--------------------------------------------------------------------
In MATLAB:
   out = run_nav_simulation(triallength, tau, plot_figs);

Arguments:
   triallength  - number of samples to simulate (e.g. 1800 = 120 seconds)
   tau           - decay constant (seconds) controlling transition probability
   plot_figs     - 1 to display figures, 0 to suppress plotting

Example:
   out = run_nav_simulation(1500, 3, 1);

--------------------------------------------------------------------
OUTPUTS
--------------------------------------------------------------------
The function returns a structure 'out' containing:
   out.U    - 5 x T matrix of DN activity
   out.v    - forward velocity (mm/s)
   out.a    - angular velocity (rad/s)
   out.x    - x positions (plume pixel coordinates)
   out.y    - y positions
   out.odor - sampled odor concentration at each timestep
   out.C    - compressed odor signal
   out.s    - behavioral state index:
                1 = baseline
                2 = goal-directed
                3 = search

If plotting is enabled, two figures are generated:
   Figure 1 – Time series of odor, state, and velocities
   Figure 2 – Trajectory of the simulated agent in plume coordinates

--------------------------------------------------------------------
EXPECTED RESULTS
--------------------------------------------------------------------
Simulations produce trajectories that reflect state transitions between
baseline, goal-directed, and search modes, similar to those observed in
Drosophila walking experiments.

--------------------------------------------------------------------
LICENSE
--------------------------------------------------------------------
This code is distributed under the GNU License (see LICENSE.txt).

--------------------------------------------------------------------
CITATION
--------------------------------------------------------------------
If you use or adapt this code, please cite:

Kathman, N., Lanz, A., Freed, J., Nagel, K. (2025). "Neural dynamics for working memory and evidence integration during olfactory navigation in Drosophila." (bioRxv).

--------------------------------------------------------------------
CONTACT
--------------------------------------------------------------------
Katherine Nagel
katherine.nagel@nyulangone.org

Nicholas Kathman
nicholas.kathman@nyulangone.org


