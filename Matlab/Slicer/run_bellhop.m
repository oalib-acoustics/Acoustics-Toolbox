function [ rv, rc ] = run_bellhop( ConfigParams, Slice )

% RUN_BELLHOP - run BELLHOP from Matlab and extract the desired results
%
% The complexity of this code is mainly because of the processing of a
% return code on different possible platforms


% set default values for the output arguments

rv = [];
rc = 0;

cwd_entry = pwd( );   % save the current working directory (so we can go back there on exit)

% extract the needed fields from the configuration and slice structures

calc_topdir  = ConfigParams.calc_topdir;
calc_subdir  = 'bellhop';
env_basename = 'temp';
exe_bellhop  = ConfigParams.exe_bellhop;
exec_cmd     = ConfigParams.exec_cmd;
run_type     = Slice.Beam.RunType( 1 );

% change to the directory tree where this calculation will be performed

path_subdir = sprintf( '%s/%s', calc_topdir, calc_subdir );

if ( exist( path_subdir, 'dir' ) )
   cd( path_subdir )
else
   error( [ mfilename, ': error: case specific calculation dir does not exist: ', path_subdir ] );
end

% remove the status file if it exists (it's used to check for certain errors)

rc_filename = [ env_basename, '.rc' ];

if ( exist( rc_filename, 'file' ) )
   delete( rc_filename )
end

% executing BELLHOP is host operating system dependent

if ismac() || isunix() % write the Bourne shell script to run BELLHOP
   
   sh_filename  = [ env_basename, '.sh'  ];
   log_filename = [ env_basename, '.log' ];
   
   fid = fopen( sh_filename, 'w' );
   
   fprintf( fid, '#!/bin/sh\n' );
   fprintf( fid, '#\n' );
   fprintf( fid, '# Be certain we are in the correct directory\n' );
   fprintf( fid, 'cd "%s"\n', path_subdir );
   fprintf( fid, '# Remove the BELLHOP status file to insure a clean start\n' );
   fprintf( fid, 'rm -f %s\n', rc_filename );
   fprintf( fid, '# Run BELLHOP and pipe any errors to the log file\n' );
   fprintf( fid, '%s "%s" %s > %s 2>&1\n', ...
      exec_cmd, exe_bellhop, env_basename, log_filename );
   fprintf( fid, '# Save the BELLHOP exit status to the status file\n' );
   fprintf( fid, 'echo $? > %s\n', rc_filename );
   fprintf( fid, 'exit 0\n');
   
   fclose( fid );
   
   % run the Bourne shell wrapper script
   
   [ status, ~ ] = system( [ '/bin/sh ', sh_filename ] );
   
elseif ispc()   % write the BATCH script to run BELLHOP
   
   bat_filename = [ env_basename, '.bat' ];
   log_filename = [ env_basename, '.log' ];
   
   fid = fopen( bat_filename, 'wt' );
   
   fprintf( fid, 'REM Be certain we are in the correct directory\n');
   fprintf( fid, 'CD "%s\\%s"\n', calc_topdir, calc_subdir );
   fprintf( fid, 'REM Remove the BELLHOP status file to insure a clean start\n' );
   fprintf( fid, 'DEL %s /q\n', rc_filename );
   fprintf( fid, 'REM Run BELLHOP and pipe any errors to the log file\n' );
   fprintf( fid, '"%s" %s > %s\n', exe_bellhop, env_basename, log_filename );
   fprintf( fid, 'REM Save the BELLHOP exit status to the status file\n' );
   fprintf( fid, 'ECHO %%ERRORLEVEL%% > %s\n', rc_filename );
   fprintf( fid, 'EXIT\n' );
   
   fclose( fid );
   
   % run the BATCH wrapper script
   
   [ status, ~ ] = system( bat_filename );
   
else
   arch_str = computer( 'arch' );
   error( [ mfilename, ': error: your platform is currently unsupported: ', arch_str ] );
end

% check for success or failure of the BELLHOP wrapper script

if status ~= 0
   rc = 1;   % bad result from call to system function
   return
else
   % check the contents of the status file for the BELLHOP exit status
   fid   = fopen( rc_filename, 'rt' );
   bh_rc = fscanf( fid, '%d');
   fclose( fid );
   
   if bh_rc ~= 0
      rc = 2;   % bad exit status from BELLHOP itself
      return
   end
end

% extract BELLHOP calculations

switch run_type
   case 'A'             % arrivals calculation, ascii format output file
      arr_filename = [ env_basename, '.arr' ];
      [ rv, Pos ] = read_arrivals_asc( arr_filename );
      
   case 'a'             % arrivals calculation, binary format output file
      arr_filename = [ env_basename, '.arr' ];
      [ rv, Pos ] = read_arrivals_bin( arr_filename );
      
   case { 'C', 'I' }    % coherent or incoherent transmission loss calculation
      shd_filename = [ env_basename, '.shd' ];
      [ PlotTitle, PlotType, freq, atten, shdPos, p ] = read_shd_bin( shd_filename );
      
      % convert from pressure to transmission loss (assumed SL = 0 dB)
      tl = abs( squeeze( p( 1, 1, :, : ) ) );
      tl( isnan( tl ) ) = 1.0d-6;
      tl( isinf( tl ) ) = 1.0d-6;
      tl( tl < 1.0d-6 ) = 1.0d-6;
      rv = -20.0 * log10( tl );		% transmission loss in dB
      
end

cd( cwd_entry )   % change back to the directory where we were upon entry
