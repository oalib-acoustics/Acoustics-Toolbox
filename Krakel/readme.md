#KRAKEN Normal modes for ocean acoustics

   Copyright (C) 2009 Michael B. Porter

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


UPDATES:


August 2010
Lots of modifications to bring KRAKEL up to date with the rest of the
Acoustic Toolbox package, and use modern Fortran capabilities. Found 4
bugs:

	1) EXTRAP was being added to k in the main routine. This is done in KRAKEN and KRAKENC as a way to add scatter losses back in. However, k was never set in KRAKEL,

	2) A factor of rho was missing in the acoustic-elastic interfaces, causing Scholte modes to disappear,

	3) The elastic-elastic interface condition was shifted two diagonals to the left in the matrix set up. This caused big errors for cases with more than one elastic layer (not counting halfspaces).

	4) The default spectral limits were set to .9 times the slowest wave speed; however, some Scholte waves are .85 times the slowest wave speed and were missed.

Added some capabilities to plotmode.m to display elastic modes.
