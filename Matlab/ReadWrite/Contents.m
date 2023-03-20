% READWRITE
%
% Files
%   read_env          - Read an environmental file
%   read_env_core     - Read the core of the environmental file
%   read_bell         - Read the rest of the environmental file

%   read_flp          - Read the field.flp file

%   readbdry3d        - Read a boundary file (3D case)

%   read_arrivals_asc - Read the ASCII format arrivals data file written by BELLHOP or BELLHOP3D
%   read_arrivals_bin - Read the BINARY format arrivals data file written by BELLHOP or BELLHOP3D
%   readrc            - Read the top and bottom reflection coefficients

%   read_modes        - Read the modes produced by KRAKEN
%   read_modes_asc    - Read the modes produced by KRAKEN
%   read_modes_bin    - Read the modes produced by KRAKEN

%   read_shd          - Read the shade file
%   read_shd_asc      - Read an ascii shade file
%   read_shd_bin      - Read a binary shade file

%   readRcvrBearings  - Read receiver bearings
%   read_flp3d        - Read the field3d.flp file
%   read_ram_tlgrid   - Read the binary tl.grid output file written by RAM
%   read_ts           - Read the time-series file
%   readpat           - Read in a source beam pattern
%   readsxsy          - Read source x and y coordinates
%   readszrz          - Read source depths and receiver depths
%   readtheta         - Read receiver angles
%   readvector        - Reads a row-vector

%   readssp2d         - Read the SSPFIL

%   readati           - Read an atifil into the workspace variables
%   readbty           - Read a  btyfil into the workspace variables
%   readr             - Read receiver ranges

%   write_env         - Write an environmental file
%   write_bell        - Write out the rest of the environmental file
%   write_env_RAM     - Write an environmental file

%   write_fieldflp    - Write a field-parameters file
%   write_fieldsflp   - Write a field-parameters file
%   write_field3dflp  - 

%   writebdry         - Write a boundary file (bathymetry or altimetry) from the workspace variables
%   writebdry3d       - Write a 3d bathymetry file for use by BELLHOP3D
%   writessp          - Write an SSP matrix

%   CmplxSSP          - Set the imaginary part of the SSP based on the frequency
%   equally_spaced    - Test whether the vector x is composed of roughly equally-spaced values
%   get_component     - extract a single component out of the stress-displacement vector

%   readssp3d         - Read the 3D SSPFIL used by Bellhop3D
