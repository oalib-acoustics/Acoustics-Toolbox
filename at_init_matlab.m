% -------------------------------------------------
%   at_init_matlab
%
%   Must be run FIRST to initialize paths for
%     running AT Matlab routines
%
%   Usually I set the path from the 'Set Path' tab in the Matlab window.
%
%   However, some people that use a command line interface to Matab find
%   this a useful alternative.
% ------------------------------------------------

% This is your Home where all AT routines are located
Home = pwd;

%addpath( fullfile(Home, 'bin' ) );             % AT binaries

Matdir = fullfile( Home, 'Matlab' );          % AT Matlab routines
addpath( Matdir );

% addpath for all extra folders/routines in Matlab
dir_list = dir( Matdir );

for j = 1 : length( dir_list )
   if ( dir_list( j ).isdir )
       % ignore CWD . and parent dir ..
       if ~strcmp( dir_list( j ).name, '.' ) && ~strcmp( dir_list( j ).name, '..' )
           addpath( fullfile( Matdir, dir_list( j ).name ) );
       end
   end
end

