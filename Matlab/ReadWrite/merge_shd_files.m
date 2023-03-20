function merge_shd_files(shd_files_in, shd_file_out)

% MERGE_SHD_FILES - function to merge broadband shd files into a single output
%
% The merge_shd_files( ) function accepts a list of broadband shd files
% as input, and merges them into a single output shd file. It presently
% only handles Matlab shd files. It checks that the receiver grids of the
% input shd files are consistent, and for any duplicate frequencies.
%
% The input shd files do not have to be in any particular order as the
% ordering of the frequencies (ascending) and the associated pressure field
% is checked and those fields will be automatically rearranged as needed to
% make them consistent. For large or numerous input shd files, it is best
% to present them in the correct order as dynamic (RAM) memory requirements
% can become significant if any reordering is required.
%
% Usage:
%
%  merge_shd_files(shd_files_in, shd_file_out)
%
%  Inputs:
%   shd_files_in  - cell array with the broadband input shd filenames
%   shd_files_out - character array containing the output shd filename
%
%  Outputs:
%   none
%
% $Id: merge_shd_files.m,v 1.4 2019/09/16 18:55:21 jcp Exp jcp $

% Check for missing input arguments

if nargin < 2
  error( [mfilename, ': one or more required input arguments are missing.'] );
end

% Sanity check on the first input argument, it should be a cell array

if ~iscell(shd_files_in)
  error( [mfilename, ': input argument shd_files_in is not a cell array.'] );
end

% Load all the input shd files to verify that they are mutually consistent

num_files = length(shd_files_in);

shd_structs = repmat( struct( 'PlotTitle', { ' ' }, ...
                                  'atten', {  0  }, ...
                                    'Pos', { struct( 'theta', { [ ] }, ...
                                                         's', { [ ] }, ...
                                                         'r', { [ ] }  ) }, ...
                                  'freq0', { [ ] }, ...
                                'freqVec', { [ ] }, ...
                               'pressure', { [ ] }, ...
                               'PlotType', { [ ] } ), num_files, 1 );

for j = 1:num_files
  shd_structs(j) = load(shd_files_in{j});
  % verify that the size of the pressure fields are consistent
  sizePressure_j = size(shd_structs(j).pressure);
  if j > 1
    % check the size of this pressure field (ignoring the number of frequencies)
    if ~isequal(sizePressure_j(2:end), sizePressure(2:end))
      error( [mfilename, ': size of pressure fields of inputs do not agree.'] );
    end
  else
    % save the size of the 1-st pressure field
    sizePressure = sizePressure_j;
  end
  % free memory that's not needed
  shd_structs(j).pressure = [ ];
end

PlotTitle = shd_structs(1).PlotTitle;
atten = shd_structs(1).atten;
Pos = shd_structs(1).Pos;
freq0 = shd_structs(1).freq0;
freqVec = shd_structs(1).freqVec;
PlotType = shd_structs(1).PlotType;

% Sanity check on the atten, Pos (src, rcv) structs, they should all agree!

num_freqs = length(freqVec);

for j = 2:num_files
  % verify the atten fields are consistent
  if ~isequal(atten, shd_structs(j).atten)
    error( [mfilename, ': atten fields of inputs do not agree.'] );
  end
  % verify the Pos fields are consistent
  if ~isequal(Pos, shd_structs(j).Pos)
    error( [mfilename, ': Pos (src, rcv) fields of inputs do not agree.'] );
  end
  % verify the PlotType fields are consistent
  if ~isequal(PlotType, shd_structs(j).PlotType)
    error( [mfilename, ': PlotType fields of inputs do not agree.'] );
  end
  % update the number of frequencies
  num_freqs = num_freqs + length(shd_structs(j).freqVec);
  % free memory that's no longer needed
  shd_structs(j).Pos = [ ];
end

% Collate all the frequencies into a single array, check for any duplicates

freqVec = zeros(num_freqs, 1);

j_pnt = 0;
for j = 1:num_files
  lenFreq = length(shd_structs(j).freqVec);
  % transfer these frequencies to get a complete list
  freqVec(j_pnt+1:j_pnt+lenFreq) = shd_structs(j).freqVec;
  % free memory that's no longer needed
  shd_structs(j).freqVec = [ ];	
  % update the pointer
  j_pnt = j_pnt + lenFreq;
end

if ~isequal(sort(freqVec), unique(freqVec))
  error( [mfilename, ': duplicate frequencies found in the inputs.'] );
end

% Load all the input shd files again, collate the pressures into a single matrix

sizePressure(1) = num_freqs;	% update the total number of frequencies
pressure = zeros(sizePressure);
j_pnt = 0;

for j = 1:num_files
  shd_structs(j) = load(shd_files_in{j});
  sizePressure = size(shd_structs(j).pressure);
  % store the pressures from this shd file in the collated pressure matrix
  pressure(j_pnt+1:j_pnt+sizePressure(1),:,:,:) = shd_structs(j).pressure;
  % free memory that's no longer needed
  shd_structs(j).pressure = [ ];
  % update the pointer
  j_pnt = j_pnt + sizePressure(1);
end

clear shd_structs;

% Sort freqVec, pressure by frequency if needed

[sortFreq, itags] = sort(freqVec);

if ~isequal(freqVec, sortFreq)
  % reorder freqVec and pressure accordingly
  freqVec = freqVec(itags);
  pressure = pressure(itags,:,:,:);
end

clear itags sortFreq

% Output the merged shd file (note: version v7.3 fmt is slower than Molasses)

save(shd_file_out, 'PlotTitle', 'atten', 'Pos', ...
                   'freq0', 'freqVec', 'pressure', 'PlotType', '-v7.3');

%
% End of merge_shd_files.m
