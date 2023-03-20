function [ rv, rc ] = run_ram( ConfigParams )

% RUN_RAM - run RAM from Matlab and extract the desired results
%
% The complexity of this code is mainly because of the processing of a
% return code on different possible platforms

% set default values for the output arguments

rv = [];
rc = 0;

cwd_entry = pwd( );   % save the current working directory (so we can go back there on exit)

% extract the needed fields from the configuration and slice structures

calc_topdir = ConfigParams.calc_topdir;
calc_subdir = 'ram';
exe_ram     = ConfigParams.exe_ram;
exec_cmd    = ConfigParams.exec_cmd;

% change to the directory tree where this calculation will be performed

path_subdir = sprintf('%s/%s', calc_topdir, calc_subdir);

if ( exist( path_subdir, 'dir' ) )
  cd( path_subdir )
else
  error( [ mfilename, ': error: case specific calculation dir does not exist: ', path_subdir ] );
end

% remove the status file if it exists (it's used to check for certain errors)

rc_filename = 'ram.rc';

if ( exist( rc_filename, 'file' ) )
  delete( rc_filename )
end

% executing RAM is host operating system dependent

if ismac() || isunix()

  % write the Bourne shell script to run RAM

  sh_filename  = 'ram.sh';
  log_filename = 'ram.log';

  fid = fopen( sh_filename, 'w' );

  fprintf( fid, '#!/bin/sh\n' );
  fprintf( fid, '#\n' );
  fprintf( fid, '# Be certain we are in the correct directory\n');
  fprintf( fid, 'cd "%s"\n', path_subdir);
  fprintf( fid, '# Remove the RAM status file to insure a clean start\n');
  fprintf( fid, 'rm -f %s\n', rc_filename);
  fprintf( fid, '# Run RAM and pipe any errors to the log file\n');
  fprintf( fid, '%s "%s" > %s 2>&1\n', exec_cmd, exe_ram, log_filename);
  fprintf( fid, '# Save the RAM exit status to the status file\n');
  fprintf( fid, 'echo $? > %s\n', rc_filename);
  fprintf( fid, 'exit 0\n');

  fclose( fid );

  % run the Bourne shell wrapper script

  [ status, ~ ] = system( [ '/bin/sh ', sh_filename ] );

elseif ispc()

  % write the BATCH script to run RAM

  bat_filename = 'ram.bat';
  log_filename = 'ram.log';

  fid = fopen( bat_filename, 'wt' );

  fprintf( fid, 'REM Be certain we are in the correct directory\n');
  fprintf( fid, 'CD "%s\\%s"\n', calc_topdir, calc_subdir);
  fprintf( fid, 'REM Remove the RAM status file to insure a clean start\n');
  fprintf( fid, 'DEL %s /q\n', rc_filename);
  fprintf( fid, 'REM Run RAM and pipe any errors to the log file\n');
  fprintf( fid, '"%s" > %s\n', exe_ram, log_filename);
  fprintf( fid, 'REM Save the RAM exit status to the status file\n');
  fprintf( fid, 'ECHO %%ERRORLEVEL%% > %s\n', rc_filename);
  fprintf( fid, 'EXIT\n');

  fclose( fid );

  % run the BATCH wrapper script

  [ status, ~ ] = system( bat_filename );

else

  arch_str = computer( 'arch' );
  error( [ mfilename, ': error: your platform is currently unsupported: ', arch_str ] );

end

% check for success or failure of the BELLHOP wrapper script

if status ~= 0
  % bad result from call to system function
  rc = 1;
  return
else
  % check the contents of the status file for the RAM exit status
  fid   = fopen( rc_filename, 'rt' );
  bh_rc = fscanf( fid, '%d' );
  fclose( fid );
  if bh_rc ~= 0
    % bad exit status from RAM itself
    rc = 2;
    return
  end
end

% extract RAM calculations for this case, this is platform independent

% rv = read_ram_tlgrid( );
[ PlotTitle, PlotType, freq, atten, Pos, pressure ] = read_ram_tlgrid( );
rv = -20 * log10( squeeze( pressure ) );
size( pressure )
size( rv )

cd( cwd_entry )   % change back to the directory where we were upon entry
