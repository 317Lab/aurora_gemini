# Aurora Gemini
Boundary driver development and postprocessing tools for Gemini3D simulations of auroral arc systems.

## Requirements
Please ensure your environment has the following:
- MATLAB $\geq$ r2024a
- GEMINI requirements found here:\
https://github.com/317Lab/gemini3d/blob/main/Readme.md
For HPC environments, avoid using hdf5 modules as they tend to interfere with hdf5 tools installed by the gemini3d build.

## Quick BASH Install
Run the following:
```sh
git clone https://github.com/317Lab/aurora_gemini.git; bash aurora_gemini/install.sh
```
Ensure that ctest outputs '100% tests passed'.

## Manual Install
From the directory `aurora_gemini/..`, run the following:
```sh
mkdir sims
git clone https://github.com/317Lab/gemini3d.git
git clone https://github.com/317Lab/mat_gemini.git
git clone https://github.com/317Lab/mat_gemini-scripts.git
cd mat_gemini
git clone https://github.com/geospace-code/matlab-stdlib.git

cd ../gemini3d
git checkout jvi_save_sig
cmake -B build
cmake --build build --parallel
ctest --test-dir build
cd ..

```
It is recommended to us the `sims` directory for your simulations. Forks from github.com/gemini3d are use to allow GEMINI to output conductivity volumes and for other, minor adjustments, e.g. reading additional configuration namelists. Forks will be updated semi-regularly.

## Sample Simulation
Inside `data/sample` is an example set of data required for a simulation of a Swarm-over-Poker event. To setup the sample simulation, start within MATLAB in aurora_gemini and run
```matlab
init
aurogem.swop.setup('data/sample', 2, 'AC')
```
The initial condition can take a few hours to run.

## Related Publications:
van Irsel, J., Lynch, K., Mule, A., Zettergren, M., (2024), Generation of top boundary conditions for 3D ionospheric models constrained by auroral imagery and plasma flow data, _Journal of Geophysical Research: Space Physics_.\
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JA032722

van Irsel, J., Lynch, K., Mule, A., Zettergren, M., Burchill, J., (2025), Data-Driven 3D Simulations of Auroral Arc Systems, _Journal of Geophysical Research: Space Physics_.\
Manuscript in preparation.
