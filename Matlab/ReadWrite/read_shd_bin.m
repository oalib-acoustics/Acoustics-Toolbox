function [ title, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_bin( varargin )

% Read a binary shade file
%
% Useage:
% ... = read_shd_bin( filename, xs, ys )
% where (xs, ys) is the source coordinate in km
% (xs, ys) are optional
%
% Output is a 4-D pressure field p( Ntheta, Nsd, Nrd, Nrr )
%
% Original version by Chris Tiemann, Feb. 2001
% Lots of mods ... mbp
% Laurel added extension to multiple source lat/longs 2011

%error( nargchk( 1, 3, nargin, 'struct' ) );
narginchk( 1, 3 )

filename = varargin{1};

% optional frequency
if nargin == 2
    freq = varargin{ 2 };
end

% optional source (x,y) coordinate
if nargin >= 3
    xs = varargin{ 2 };
    ys = varargin{ 3 };
else
    xs = NaN;
    ys = NaN;
end

%%
fid = fopen( filename, 'rb' );
if ( fid == -1 )
    error( 'read_shd_bin.m: No shade file with that name exists' );
end

recl     = fread( fid,  1, 'int32' );     %record length in bytes will be 4*recl
title    = fread( fid, 80, '*char' )';

fseek( fid, 4 * recl, -1 ); %reposition to end of first record
PlotType = fread( fid, 10, '*char'   );
PlotType = PlotType';

fseek( fid, 2 * 4 * recl, -1 ); %reposition to end of second record
Nfreq  = fread( fid, 1, 'int32'   );
Ntheta = fread( fid, 1, 'int32'   );

Nsx    = fread( fid, 1, 'int32'   );
Nsy    = fread( fid, 1, 'int32'   );
Nsz    = fread( fid, 1, 'int32'   );

Nrz    = fread( fid, 1, 'int32'   );
Nrr    = fread( fid, 1, 'int32'   );
freq0  = fread( fid, 1, 'float64' );
atten  = fread( fid, 1, 'float64' );

fseek( fid, 3 * 4 * recl, -1 ); %reposition to end of record 3
freqVec = fread( fid, Nfreq, 'float64' );

fseek( fid, 4 * 4 * recl, -1 ); %reposition to end of record 4
Pos.theta   = fread( fid, Ntheta, 'float64' );

if ( PlotType( 1 : 2 ) ~= 'TL' )
    fseek( fid, 5 * 4 * recl, -1 ); %reposition to end of record 5
    Pos.s.x     = fread( fid, Nsx, 'float64' );
    
    fseek( fid, 6 * 4 * recl, -1 ); %reposition to end of record 6
    Pos.s.y     = fread( fid, Nsy, 'float64' );
else   % compressed format for TL from FIELD3D
    fseek( fid, 5 * 4 * recl, -1 ); %reposition to end of record 5
    Pos.s.x     = fread( fid, 2,    'float64' );
    Pos.s.x     = linspace( Pos.s.x( 1 ), Pos.s.x( end ), Nsx );
    
    fseek( fid, 6 * 4 * recl, -1 ); %reposition to end of record 6
    Pos.s.y     = fread( fid, 2,    'float64' );
    Pos.s.y     = linspace( Pos.s.y( 1 ), Pos.s.y( end ), Nsy );
end

fseek( fid, 7 * 4 * recl, -1 ); %reposition to end of record 7
Pos.s.z = fread( fid, Nsz, 'float32' );

fseek( fid, 8 * 4 * recl, -1 ); %reposition to end of record 8
Pos.r.z = fread( fid, Nrz, 'float32' );

fseek( fid, 9 * 4 * recl, -1 ); %reposition to end of record 9
Pos.r.r = fread( fid, Nrr, 'float64' );
Pos.r.r = Pos.r.r';   % make it a row vector

%%
% Each record holds data from one source depth/receiver depth pair

switch PlotType
    case 'rectilin  '
        pressure = zeros( Ntheta, Nsz, Nrz, Nrr );
        Nrcvrs_per_range = Nrz;
    case 'irregular '
        pressure = zeros( Ntheta, Nsz,   1, Nrr );
        Nrcvrs_per_range = 1;
    otherwise
        pressure = zeros( Ntheta, Nsz, Nrz, Nrr );
        Nrcvrs_per_range = Nrz;
end

%%

if isnan( xs )    % Just read the first xs, ys, but all theta, sz, and rz
    % get the index of the frequency if one was selected
    ifreq = 1;
    if exist( 'freq', 'var' )
       freqdiff = abs( freqVec - freq );
       [ ~, ifreq ] = min( freqdiff );
    end

    for itheta = 1 : Ntheta
        %disp( [ 'Reading data for receiver at bearing ' num2str( itheta ) ' of ' num2str( Ntheta ) ] );
        for isz = 1 : Nsz
            %disp( [ 'Reading data for source at depth ' num2str( isz ) ' of ' num2str( Nsz ) ] );
            for irz = 1 : Nrcvrs_per_range
                %disp( [ 'Reading data for receiver at depth ' num2str( irz ) ' of ' num2str( Nrcvrs_per_range ) ] );
                recnum = 10 + ( ifreq  - 1 ) * Ntheta * Nsz * Nrcvrs_per_range + ...
                              ( itheta - 1 )          * Nsz * Nrcvrs_per_range + ...
                              ( isz    - 1 )                * Nrcvrs_per_range + ...
                                irz    - 1;

                status = fseek( fid, recnum * 4 * recl, -1 ); % Move to end of previous record
                if ( status == -1 )
                    error( 'Seek to specified record failed in read_shd_bin' )
                end
                
                temp = fread( fid, 2 * Nrr, 'float32' );    % Read complex data
                pressure( itheta, isz, irz, : ) = temp( 1 : 2 : 2 * Nrr ) + 1i * temp( 2 : 2 : 2 * Nrr );
                % Transmission loss matrix indexed by  theta x sd x rd x rr
                
            end
        end
    end
else              % read for a source at the desired x, y, z.
    
    xdiff = abs( Pos.s.x - xs * 1000. );
    [ ~, idxX ] = min( xdiff );
    ydiff = abs( Pos.s.y - ys * 1000. );
    [ ~, idxY ] = min( ydiff );
    
    % show the source x, y that was found to be closest
    % [ Pos.s.x( idxX ) Pos.s.y( idxY ) ]
    for itheta = 1 : Ntheta
        for isz = 1 : Nsz
            % disp( [ 'Reading data for source at depth ' num2str( isd ) ' of ' num2str( Nsd ) ] )
            for irz = 1 : Nrcvrs_per_range
                recnum = 10 + ( idxX   - 1 ) * Nsy * Ntheta * Nsz * Nrcvrs_per_range + ...
                              ( idxY   - 1 )       * Ntheta * Nsz * Nrcvrs_per_range + ...
                              ( itheta - 1 )                * Nsz * Nrcvrs_per_range + ...
                              ( isz    - 1 )                      * Nrcvrs_per_range + irz - 1;
                status = fseek( fid, recnum * 4 * recl, -1 ); % Move to end of previous record
                if ( status == -1 )
                    error( 'Seek to specified record failed in read_shd_bin' )
                end
                
                temp = fread( fid, 2 * Nrr, 'float32' );    %Read complex data
                pressure( itheta, isz, irz, : ) = temp( 1 : 2 : 2 * Nrr ) + 1i * temp( 2 : 2 : 2 * Nrr );
                % Transmission loss matrix indexed by  theta x sd x rd x rr
                
            end
        end
    end
end

fclose( fid );

