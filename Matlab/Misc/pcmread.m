
function [ Y, bits ] = pcmread( pcmfile )

% reads a 'PCM' file (intel, little endian)
%
% usage: [ Y, bits ] = pcmread( pcmfile );
% mbp 6/00

if ~isstr( pcmfile )
  error( 'pcmfile must be a string.' );
end

% add default extension
if isempty( findstr( pcmfile, '.' ) )
  pcmfile = [ pcmfile '.pcm' ];
end

fid_in = fopen( pcmfile, 'r', 'ieee-le' );
if ( fid_in == -1 )
  error('Can''t open PCM file for input.');
end

[ Y, bits ] = fread( fid_in, inf, 'int16' );
fclose( fid_in );

% scaled to [-1, +1]

Y = Y / 2^15;

