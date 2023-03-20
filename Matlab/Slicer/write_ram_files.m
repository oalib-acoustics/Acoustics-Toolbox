function write_ram_files( ConfigParams, Slice )

% WRITE_RAM_FILES - write all input files needed to run RAM

cwd_entry = pwd( );   % save the current working directory

% extract the needed fields from the configuration and slice structures

calc_topdir = ConfigParams.calc_topdir;
calc_subdir = 'ram';

cd( calc_topdir );    % change to top of the directory tree where all model calculations are done

% check if the sub-directory for this specific case exists

path_subdir = sprintf( '%s/%s', calc_topdir, calc_subdir );

if ~exist( path_subdir, 'dir' )
  mkdir( calc_subdir )   % sub-directory for this case doesn't exist, create it
end

cd( calc_subdir ) % change to the sub-directory for this case

% write the RAM input file

write_ram_in( Slice.TitleEnv, Slice.freq, Slice.bathy_rd, ...
              Slice.SSP, Slice.Bdry, Slice.Pos );

cd( cwd_entry );  % change back to the working directory where we were on entry
