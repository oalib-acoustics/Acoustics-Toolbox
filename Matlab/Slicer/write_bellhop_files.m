function write_bellhop_files( ConfigParams, Slice )

% WRITE_BELLHOP_FILES - write all input files needed to run BELLHOP

cwd_entry = pwd( );   % save the current working directory

% extract the needed fields from the configuration and slice structures

calc_topdir = ConfigParams.calc_topdir;
calc_subdir = 'bellhop';
envfile     = 'temp';

% change to top of the directory tree where all model calculations are done

cd( calc_topdir );

% check if the sub-directory for this specific case exists

path_subdir = sprintf( '%s/%s', calc_topdir, calc_subdir );

if ~exist( path_subdir, 'dir' )
  mkdir( calc_subdir )  % sub-directory doesn't exist, create it
end

cd( calc_subdir )   % change to the sub-directory for this case

% if needed, write the range dependent SSP file

if ( ~isempty( Slice.SSP_range ) )
  ssp_filename = [ envfile, '.ssp' ];
  writessp( ssp_filename, Slice.SSP_range, Slice.SSP_mat );
end

% if needed, write the bathymetry file

if ( ~isempty( Slice.bathy_rd ) )
  bty_filename = [ envfile, '.bty' ];
  writebty( bty_filename, Slice.bathy_interp, Slice.bathy_rd );
end

% write the BELLHOP input file

env_filename = [ envfile, '.env' ];

write_env( env_filename, 'BELLHOP', Slice.TitleEnv, ...
           Slice.freq, Slice.SSP, Slice.Bdry, Slice.Pos, ...
           Slice.Beam, Slice.cInt, Slice.RMax );

cd( cwd_entry );   % change back to the working directory where we were on entry

