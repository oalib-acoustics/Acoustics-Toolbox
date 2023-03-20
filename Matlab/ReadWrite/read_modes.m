function [ Modes ] = read_modes( modfil, freq, modes )

% Read the modes produced by KRAKEN
% usage:
%    [ Modes ] = read_modes( filename, freq, modes )
%
% filename should include the extension
% freq is the frequency
% modes is an optional vector of mode indices

% mbp, May 2001

% identify the file type

% !!! put the following line in your calling routine, for the first call
% clear read_modes_bin % to force rewind to beginning of mode file

%endchar = length( modfil );

% if ( endchar >= 6 )
%    switch modfil( 1 : 6 )
%       case 'MODFIL'
%          ext = '.mod';
%       case 'MOAFIL'
%          ext = '.moa';
%    end
% % else
% %    FileType = 'mod';
% %    if ( endchar >= 4 )
% %       FileType = lower( filename( endchar-2 : endchar ) )
% %    end
% end

[ ~, FileRoot, ext ] = fileparts( modfil );

if ( isempty( ext ) )
  ext = '.mod';   % use this as the default extension if none specified
else
   % handle cases where the modfil is FileRoot.mod.mat
   % above fileparts call only pulls off the .mat part as the extension
   if ( ext == '.mat' )
      [ ~, FileRoot, ~ ] = fileparts( FileRoot );
      ext = '.mod.mat';
   end
end

modfil = [ FileRoot ext ];

% read the modal data

switch ext
   case '.mod' % binary format
      if nargin == 2
         Modes = read_modes_bin( modfil, freq );
      else
         Modes = read_modes_bin( modfil, freq, modes );
      end
   case '.mod.mat' % Matlab format
      load( modfil );
   case '.moa' % ascii format
      if nargin == 1
         Modes = read_modes_asc( modfil );
      else
         Modes = read_modes_asc( modfil, modes );
      end
      
   otherwise
      %errordlg( 'Unrecognized file extension', 'Warning' )
      error( 'read_modes.m: Unrecognized file extension' )
end

% identify the index of the frequency closest to the user-specified value
freqdiff = abs( Modes.freqVec - freq );
[ ~, freq_index ] = min( freqdiff );

%%
% calculate wavenumbers in halfspaces (if there are any modes)

if ( Modes.M ~= 0 )
   
   if ( Modes.Top.BC == 'A' )   % top
      Modes.Top.k2     = ( 2 * pi * Modes.freqVec( 1 )  / Modes.Top.cp )^2;
      gamma2           = Modes.k .^ 2 - Modes.Top.k2;
      Modes.Top.gamma  = PekerisRoot( gamma2 );   % vertical wavenumber
      Modes.Top.phi    = Modes.phi( 1, : );       % mode value at halfspace
   else
      Modes.Top.rho   = 1.0;
      Modes.Top.gamma = zeros( size( Modes.k ) );
      Modes.Top.phi   = zeros( size( Modes.phi( 1, : ) ) );
   end
   
   if ( Modes.Bot.BC == 'A' )   % bottom
      Modes.Bot.k2    = ( 2 * pi * Modes.freqVec( freq_index ) / Modes.Bot.cp )^2;
      gamma2          = Modes.k .^ 2 - Modes.Bot.k2;
      Modes.Bot.gamma = PekerisRoot( gamma2 );    % vertical wavenumber
      Modes.Bot.phi   = Modes.phi( end, : );      % mode value at halfspace
   else
      Modes.Bot.rho   = 1.0;
      Modes.Bot.gamma = zeros( size( Modes.k ) );
      Modes.Bot.phi   = zeros( size( Modes.phi( end, : ) ) );
   end
end