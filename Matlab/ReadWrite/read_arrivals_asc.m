function [ Arr, Pos ] = read_arrivals_asc( ARRFile )

% Read the ASCII format arrivals data file written by BELLHOP or BELLHOP3D
%
% Usage: [ Arr, Pos ] = read_arrivals_asc( ARRFile );
%
% Arr is a structure array containing all the arrivals information
% Pos is a structure containing the positions of source and receivers
%
% ARRFile is the name of the ASCII format Arrivals File
%
% mbp Sep 1996
% jcp Jul 2018

fid = fopen( ARRFile, 'r' );	% open the file

if ( fid == -1 )
   error( [ mfilename, ': Arrivals file cannot be opened' ] )
end

% read the 2D/3D flag to determine the file format

flag = fscanf( fid, '%s',  1 );

% check for the case of erroneously reading a BINARY format arrivals file

if ~strcmp( flag, '''2D''' ) && ~strcmp( flag, '''3D''' )
   error( [ mfilename, ': not an ASCII format Arrivals file?' ] )
end

% proceed accordingly for the Bellhop 2D vs 3D format

if strcmp( flag, '''2D''' )
   
   % read the 2D format % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
   
   Pos.freq = fscanf( fid, '%f',  1 );	% acoustic frequency
   
   Nsz     = fscanf( fid, '%i',  1 );	% number of source   depths
   Pos.s.z = fscanf( fid, '%f', Nsz );	% source   depths
   
   Nrz     = fscanf( fid, '%i',  1 );	% number of receiver depths
   Pos.r.z = fscanf( fid, '%f', Nrz );	% receiver depths
   
   Nrr     = fscanf( fid, '%i',  1 );	% number of receiver ranges
   Pos.r.r = fscanf( fid, '%f', Nrr );	% receiver ranges
   
   % pre-allocate memory for the Arr arrivals structure array
   
   Arr = repmat( struct( 'Narr', { int16(0) }, ...
      'A',             { single(1.0i) }, ...
      'delay',         { single(1.0i) }, ...
      'SrcDeclAngle',  { single(0.0) }, ...
      'RcvrDeclAngle', { single(0.0) }, ...
      'NumTopBnc',     { int16(0) }, ...
      'NumBotBnc',     { int16(0) } ), Nrr, Nrz, Nsz );
   % loop over sources, rcv depths and rcv ranges to read all arrival info
   
   for isd = 1 : Nsz
      
      % read the maximum number of arrivals for this source
      Narrmx2 = fscanf( fid, '%i', 1 );
      
      for irz = 1 : Nrz
         
         for irr = 1 : Nrr
            
            % read the number of arrivals at this receiver
            Narr = fscanf( fid, '%i', 1 );
            Arr( irr, irz, isd ).Narr = int16( Narr );
            
            % read and store all the arrivals, if there are any
            if Narr > 0
               
               da = single( fscanf( fid, '%f', [ 8, Narr ] ) );
               
               Arr( irr, irz, isd ).A = single( ...
                  da( 1, 1 : Narr ) .* exp( 1.0i * da( 2, 1 : Narr ) * pi / 180.0 ) );
               Arr( irr, irz, isd ).delay = single( ...
                  da( 3, 1 : Narr ) + 1.0i * da( 4, 1 : Narr ) );
               Arr( irr, irz, isd ).SrcDeclAngle  = da( 5, 1 : Narr );
               Arr( irr, irz, isd ).RcvrDeclAngle = da( 6, 1 : Narr );
               Arr( irr, irz, isd ).NumTopBnc     = int16( da( 7, 1 : Narr ) );
               Arr( irr, irz, isd ).NumBotBnc     = int16( da( 8, 1 : Narr ) );
               
            end	% if any arrivals
            
         end	% next receiver range
         
      end	% next receiver depth
      
   end	% next source
   
else	% end of read 2D file format data
   
   % read the 3D format % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
   
   Pos.freq    = fscanf( fid, '%f',  1 );		% acoustic frequency
   
   Nsx         = fscanf( fid, '%i',  1 );		% number of source x-coordinates
   Pos.s.x     = fscanf( fid, '%f', Nsx );   %           source x-coordinates
   
   Nsy         = fscanf( fid, '%i',  1 );		% number of source y-coordinates
   Pos.s.y     = fscanf( fid, '%f', Nsy );   %           source y-coordinates
   
   Nsz         = fscanf( fid, '%i',  1 );		% number of source z-coordinates
   Pos.s.z     = fscanf( fid, '%f', Nsz );   %           source z-coordinates
   
   Nrz         = fscanf( fid, '%i',  1 );		% number of receiver z-coordinates
   Pos.r.z     = fscanf( fid, '%f', Nrz );   %           receiver z-coordinates
   
   Nrr         = fscanf( fid, '%i',  1 );		% number of receiver ranges
   Pos.r.r     = fscanf( fid, '%f', Nrr );   % receiver ranges
   
   Nrtheta     = fscanf( fid, '%i',  1 );		% number of receiver bearings
   Pos.r.theta = fscanf( fid, '%f', Nrtheta );	%        receiver bearings
   
   % pre-allocate memory for the Arr arrivals structure array
   
   Arr = repmat( struct( 'Narr', { int16(0) }, ...
      'A',             { single(1.0i) }, ...
      'delay',         { single(1.0i) }, ...
      'SrcDeclAngle',  { single(0.0) }, ...
      'SrcAzimAngle',  { single(0.0) }, ...
      'RcvrDeclAngle', { single(0.0) }, ...
      'RcvrAzimAngle', { single(0.0) }, ...
      'NumTopBnc',     { int16(0) }, ...
      'NumBotBnc',     { int16(0) } ), Nrr, Nrz, Nrtheta, Nsz );
   
   % loop over sources, rcv bearings, depths and ranges to read all arrival info
   
   for isd = 1 : Nsz
      
      % read the maximum number of arrivals for this source
      Narrmx2 = fscanf( fid, '%i', 1 );
      
      for irtheta = 1 : Nrtheta
         
         for irz = 1 : Nrz
            
            for irr = 1 : Nrr
               
               % read the number of arrivals at this receiver
               Narr = fscanf( fid, '%i', 1 );
               Arr( irr, irz, irtheta, isd ).Narr = int16( Narr );
               
               % read and store all the arrivals, if there are any
               if Narr > 0
                  
                  da = single( fscanf( fid, '%f', [ 10, Narr ] ) );
                  
                  Arr( irr, irz, irtheta, isd ).A = single( ...
                     da( 1, 1:Narr ) .* exp( 1.0i * da( 2, 1:Narr ) * pi/180.0 ) );
                  Arr( irr, irz, irtheta, isd ).delay = single( ...
                     da( 3, 1:Narr ) + 1.0i * da( 4, 1:Narr ) );
                  Arr( irr, irz, irtheta, isd ).SrcDeclAngle  = da( 5, 1:Narr );
                  Arr( irr, irz, irtheta, isd ).SrcAzimAngle  = da( 6, 1:Narr );
                  Arr( irr, irz, irtheta, isd ).RcvrDeclAngle = da( 7, 1:Narr );
                  Arr( irr, irz, irtheta, isd ).RcvrAzimAngle = da( 8, 1:Narr );
                  Arr( irr, irz, irtheta, isd ).NumTopBnc = int16(da( 9, 1:Narr));
                  Arr( irr, irz, irtheta, isd ).NumBotBnc = int16(da(10, 1:Narr));
                  
               end	% if any arrivals
               
            end	% next receiver range
            
         end	% next receiver depth
         
      end	% next receiver bearing
      
   end	% next source depth
   
end	% end of read 3D file format data

% close the arrivals file

fclose( fid );
