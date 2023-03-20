function pcmwrite( Y, pcmfile )

% PCMWRITE write a 'PCM' file (intel, little endian)
%
%    pcmwrite( Y, PCMFILE ); writes data Y to a PCM file specified by the file name
%    PCMFILE. Amplitude values outside the range [-1, +1] are clipped.
%    mbp 6/00

if ~isstr( pcmfile )
  error( 'pcmfile must be a string.' );
end

% add default extension
if isempty( findstr( pcmfile, '.' ) )
  pcmfile = [ pcmfile '.pcm' ];
end

% open output file

fid_out = fopen( pcmfile, 'w', 'ieee-le' );
if ( fid_out == -1 )
  error('Can''t open PCM file for output.');
end

% clip

ii = find( Y < -1 );
if ~isempty( ii ),
   warning( 'data values < -1 clipped' );
   Y( ii ) = -1;
end

ii = find( Y > 1 );
if ~isempty( ii ),
   warning( 'data values > +1 clipped' );
   Y( ii ) = 1;
end

fwrite( fid_out, 2^15 * Y, 'int16' );
fclose( fid_out );
