function phi = get_component( Modes, comp )

% extract a single component out of the stress-displacement vector
% mbp 8/14/2010
%
% Note that KRAKEL tabulates the modes at the finite difference grid
% while the other codes subtabulate to the specified receiver depths

% step through each medium

jj = 1;
k  = 1;

for medium = 1 : Modes.Nmedia
   
   for ii = 1 : length( Modes.z )
   % choose the following line instead of the one above for a KRAKEL run
   %for ii = 1 : Modes.N( medium ) + 1
      % Jump out if the modes are not tabulated in the elastic media (KRAKEN or KRAKENC run)
      
      if ( k > size( Modes.phi, 1 ) )
         return
      end
      
      switch Modes.Mater( medium, : )
         case 'ACOUSTIC'
            phi( jj, : ) = Modes.phi( k, : );
            k = k + 1;
         case 'ELASTIC '
            switch comp
               case 'H'
                  phi( jj, : ) = Modes.phi( k    , : );
                  
               case 'V'
                  phi( jj, : ) = Modes.phi( k + 1, : );
                  
               case 'T'
                  phi( jj, : ) = Modes.phi( k + 2, : );
                  
               case 'N'
                  phi( jj, : ) = Modes.phi( k + 3, : );
               otherwise
                  disp( 'Fatal Error in get_component: Unknown component' )
            end   % switch comp
            k = k + 4;
            
         otherwise
            disp( 'Fatal Error in get_component: Unknown material type' )
            
      end   % switch Mater
      jj = jj + 1;
   end   % next ii
end   % next medium
