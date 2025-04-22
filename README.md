# Aurora Gemini
Boundary driver development and postprocessing tools for Gemini3D simulations of auroral arc systems.

## Quick BASH Install
Run the following:
```sh
bash install_aurogem.sh
```

## Manual Install
From the directory `aurora_gemini/..` run the following:
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

## Related Publications:
van Irsel, J., Lynch, K., Mule, A., Zettergren, M., (2024), Generation of top boundary conditions for 3D ionospheric models constrained by auroral imagery and plasma flow data, _Journal of Geophysical Research: Space Physics_.\
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JA032722

van Irsel, J., Lynch, K., Mule, A., Zettergren, M., Burchill, J., (2025), Data-Driven 3D Simulations of Auroral Arc Systems, _Journal of Geophysical Research: Space Physics_.\
Manuscript in preparation.
