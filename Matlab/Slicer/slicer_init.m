function [ ConfigParams, BathyData ] = slicer_init

% SET_SLICER_CONFIGPARAMS - set various slicer parameters

% Copies the supported configuration parameters to the ConfigParams struct.

calc_topdir = [ pwd(), filesep(), 'Work' ];     % top level directory for calculations


%%

% construct the full path to the BELLHOP executable and verify it really exists

%exe_bellhop = [ at_path, filesep(), 'bin', filesep(), 'bellhop.exe' ];
exe_bellhop = which( 'bellhop.exe' );

if ~exist( exe_bellhop, 'file' )
   error( [ mfilename, ': unable to locate the BELLHOP executable: ', exe_bellhop ] );
end

ConfigParams.exe_bellhop = exe_bellhop;

% construct the full path to the RAM executable and verify it really exists

%exe_ram = [ at_path, filesep(), 'bin', filesep(), 'ram.exe' ];
exe_ram = which( 'ram.exe' );

if ~exist( exe_ram, 'file' )
   error( [ mfilename, ': unable to locate the RAM executable: ', exe_ram ] );
end

ConfigParams.exe_ram = exe_ram;
%%

% top directory

if ~exist( calc_topdir, 'dir' )
   error( [ mfilename, ': the specified ''calc_topdir'' directory: ', ...
      calc_topdir, ' does not exist' ] );
end

ConfigParams.calc_topdir = calc_topdir;
%%

% configure how executables from the Acoustics Toolbox are launched

if ismac() || isunix()   % Mapple and Unix/Linux systems
   
   [ status, result ] = system( 'which ionice' );
   if status == 0
      % use ionice when available
      ionice   = strrep( result, sprintf( '\r' ), '' );
      ionice   = strrep( ionice, sprintf( '\n' ), '' );
      exec_cmd = [ 'env -i ', ionice, ' -c 3' ];
   else
      exec_cmd = 'env -i';
   end
   
else            % other operating systems
   exec_cmd = '';
end

ConfigParams.exec_cmd = exec_cmd;

%%
% Load in the Bathymetry for the area of interest

BathyData = LoadBathymetry( 'Stellwagen.xyz' );

