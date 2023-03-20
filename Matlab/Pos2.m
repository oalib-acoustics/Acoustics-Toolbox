classdef readsdrd
    
    % Variable 'Pos' is a structure:
    % Pos.r.z = vector of receiver depths
    % Pos.Nrd     = number of receiver depths
    % Pos.s.z = vector of source depths
    % Pos.Nsd     = number of source depths
    
    properties
        Ntheta
        Nsd
        Nrd
        Nrr
        
        theta
        s
        r
    end
    
    methods
        function Pos = readsdrd( fid )
            
            % Read in source depths and receiver depths
            
            % source info
            Pos.Nsd = fscanf( fid, '%i', 1 );
            fprintf( '\nNumber of sources   = %i \n', Pos.Nsd )
            fgetl( fid );
            
            Pos.s.z = fscanf( fid, '%f', Pos.Nsd );
            
            if ( Pos.Nsd < 10 )
                fprintf( '%8.2f  \n', Pos.s.z )   % print all the depths
            else
                fprintf( '%8.2f ... %8.2f \n', Pos.s.z( 1 ), Pos.s.z( end ) ) % print first, last depth
            end
            fgetl( fid );
            
            if Pos.Nsd > 2
                Pos.s.z = linspace( Pos.s.z( 1 ), Pos.s.z( 2 ), Pos.Nsd ); % generate vector of receiver depths
                warning( 'Producing source depths by interpolation between sd(1) and sd(2)' )
            end
            
            % receiver info
            Pos.Nrd = fscanf( fid, '%i', 1 );
            fprintf( '\nNumber of receivers = %i \n', Pos.Nrd )
            fgetl( fid );
            
            Pos.r.z = fscanf( fid, '%f', Pos.Nrd );
            
            if ( Pos.Nrd < 10 )
                fprintf( '%8.2f  \n', Pos.r.z )   % print all the depths
            else
                fprintf( '%8.2f ... %8.2f \n', Pos.r.z( 1 ), Pos.r.z( end ) ) % print first, last depth
            end
            
            disp( '  ' )
            
            if Pos.Nrd > 2
                Pos.r.z = linspace( Pos.r.z( 1 ), Pos.r.z( 2 ), Pos.Nrd ); % generate vector of receiver depths
                warning( 'Producing receiver depths by interpolation between rd(1) and rd(2)' )
            end
            fgetl( fid );
            fprintf( '\n' )
        end
    end
end
