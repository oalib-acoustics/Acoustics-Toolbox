# acoustic-toolbox
public Github port of the official Acoustics Toolbox provided at [oalib.hlsresearch.com/AcousticsToolbox](https://oalib.hlsresearch.com/AcousticsToolbox)

### changes to run
My system is Ubuntu 22.04 and Matlab 2020b; I had to fix where the MATLAB gfortran library was pointing towards. You can read more about this [here](https://stackoverflow.com/questions/9628273/libgfortran-version-gfortran-1-4-not-found). Check what package of MATLAB and gfortran you are using to modify the following bash commands as needed.

	cd /usr/local/MATLAB/R2020b/sys/os/glnxa64
	sudo ln -sf /usr/lib/x86_64-linux-gnu/libgfortran.so.5.0.0 libgfortran.so.5

### how to use this repo

1. Download this repo by cloning it or as a zipfile. This contains executables! If you want to re-compile for yourself (command line):

		rm -rf bin/
		mkdir bin/
		make all
		make install

2. Open the startup.m file for MATLAB. If you've never done that before, you can do `open startup` in MATLAB's command window.
3. Add `addpath(genpath('path/to/acoustic-toolbox'));`
4. Test in MATLAB!

		which kraken
		which kraken.exe
		kraken MunkK
		
### quick description

#### BELLHOP
ray tracing program

#### KRAKEN
_real_ normal modes

#### KRAKENC
_complex_ normal modes

#### RAMGEO
low frequency transmission loss in range-dependent envionments w/ _fluid_ seabeds

#### RAMSGEO
low frequency transmission loss in range-dependent envionments w/ _elastic_ seabeds

#### SCOOTER
depth dependent Greens function for range-independent environments

It should be noted that KRAKEN and KRAKENC access FIELD; SCOOTER access FIELDS.
