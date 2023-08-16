% deletes all temporary files in the tests directory

Noise_cases = [ ...
   'halfspace       '; ...
   'Pekeris         '; ...
   'KupermanIngenito'; ...
   'ATOC            ';
   ];

for icase = 1: size( Noise_cases, 1 )
   directory = deblank( Noise_cases( icase, : ) )
   eval( [ 'cd ' directory ] );
   clean
   cd ..
end

