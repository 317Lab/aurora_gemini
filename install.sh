#!/bin/bash

#read -e -p "Install location: " install_direc
#install_direc="${install_direc/#\~/$HOME}"
#install_direc="$(realpath "$install_direc" 2>/dev/null)"

#if [ -d "$install_direc" ]; then
#	cd "$install_direc"
#else
#	echo "Directory does not exist: $install_direc" >&2
#	exit 1
#fi

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

