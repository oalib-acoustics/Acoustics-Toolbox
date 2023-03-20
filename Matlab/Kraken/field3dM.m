function field3dM( FileRoot )

'This is not currently functioning'
'There seem to be some minor errors related to the new braodaband freqVec in the modefiles'
error( 'Not currently functioning' )

% calculates the field using modes produced by KRAKEN
% parameters read from field3d.flp
%
% usage: field3d( filename )
% mbp November 2012 based on earlier Fortran version

clear read_modes_bin % to force rewind to beginning of mode file
MaxM = 5001;
iCorner = [ 1, 2 ; 3, 2; 3, 1 ];

[ Title, Option, Mlimit, R, Rmin, Rmax, Pos, sx, sy, Nsx, Nsy, theta, Ntheta, NNodes, NElts, x, y, ModeFileName, Node ] = ...
    read_flp3d( FileRoot );

sd  = Pos.s.z;
rd  = Pos.r.z;   % receiver depths
Nsd = length( sd );
Nrd = length( rd );
Nr  = length( R );   % receiver ranges
Pos.r.r = R;
Pos.s.x     = sx;
Pos.s.y     = sy;
Pos.theta   = theta;
Pos.Ntheta  = Ntheta;
Pos.Nrr     = Nr;

% optionally read in a source beam pattern
%SBP = Opt( 7 : 7 );
%SrcBmPat = readpat( FileRoot, SBP )

AdjacentElement = zeros( 3, NElts );   % need to initialize here so that array is in the scope of this routine
BuildAdjacentElementTable

P = zeros( Ntheta, Nr );
pressure = zeros( Nsx, Nsy, Nsd, Nrd, Ntheta, Nr );

% MAIN Loop: over receiver depth

for ird1 = 1 : Nrd
    [ k, phiR, M, freq, TitleEnv, iSet ] = ReadAllModeSets( rd( ird1 ) );   % Get modes at the receiver depth
    
    for isx = 1 : Nsx   % source x coordinate
        
        for isy = 1 : Nsy
            IElementSource = IdentifySourceElement( sx( isx ), sy( isy ) );
            % write( *, * ) 'source position, element', sx( isx ), sy( isy ), IElementSource
            
            % Read modes at all the source depths, for the corner nodes of the source triangle
            
            MS   = zeros( 3, 1 );
            phiS = zeros( MaxM, Nsd, 3 );
            
            for ICorner1 = 1 : 3
                INode = Node( ICorner1, IElementSource );
                temp = char( ModeFileName{ INode } );
                if ( temp( 1 : 5 ) ~= 'DUMMY' )
                    clear read_modes_bin
                    [ Modes ] = read_modes( [ temp '.mod' ], 0.0 );
                    phiS( 1 : Modes.M, :, ICorner1 ) = interp1( Modes.z, Modes.phi, sd' );   % subsample phi onto the receiver grid
                    MS( ICorner1 ) = Modes.M;
                    freq           = Modes.freqVec( 1 );
                else
                    MS( ICorner1 ) = 0;   % Dummy mode file, set number of modes to zero
                end
            end
            
            Msource = min( MS );   % limited by the # at the corners of the surrounding triangle
            MProp   = Msource;
            
            for isd1 = 1 : Nsd   % source depths
                
                P = 0;
                if ( Msource > 0 )
                    
                    phiST( 1 : Msource, : ) = phiS( 1 : Msource, isd1, : ); % Get modes at the source depth for this particular index
                    
                    % the appropriate routine to evaluate the field
                    Option( 1 : 3 ) = 'PDQ';   % force us of PDQ as this is the only one implemented
                    switch ( Option( 1 : 3 ) )
                        case ( 'STD' )  % STandDard (ignores hor. refraction)
                            %Evaluate3D(  freq, k, phiR, phiST, M, MaxM, IElementSource, sx( isx ), sy( isy ), ...
                            %     theta, Ntheta, Rmin, Rmax, Nr, Mlimit, P );
                        case ( 'PDQ' )  % Pretty Darn Quick (ignores and ignores)
                            P = EvaluatePDQ( freq, k, phiR, phiST, M, MaxM, IElementSource, sx( isx ), sy( isy ), ...
                                theta, Ntheta, Rmin, Rmax, Nr, Mlimit, AdjacentElement );
                        case ( 'GBT' )  % Gaussian Beam Tracing (include hor. refraction)
                            %EvaluateGB(  k, phiR, phiST, M, MaxM, IElementSource, sx( isx ), sy( isy ), ...
                            %     theta, Ntheta, Rmin, Rmax, Nr, Mlimit, Option, P, freq );
                        otherwise
                            %ERROUT( PRTFile, 'F', 'FIELD3D', 'Unknown option' );
                    end
                end
                
                % save the field
                pressure( isx, isy, isd1, ird1, :, : ) = P;
                
            end   % source depth
        end % Source_y
    end % Source_x
end % RcvrDepth

PlotTitle   = Title;
PlotType    = 'rectilin  ';
freq0       = 0.0; % Modes.freq;
atten       = 0.0;
% Pos.r.r = R;

save( [ FileRoot '.shd.mat' ], 'PlotTitle', 'PlotType', 'freq0', 'atten', 'Pos', 'pressure' )

%%
    function BuildAdjacentElementTable
        
        % Constructs a table AdjacentElement( iElt, iside ) which gives the
        % element number that shares iside with element number iElt
        
        % The logic is slightly different to the Fortran version because
        % Matlab does not have the CYCLE command to specify a jump to a new step in
        % an outer loop
        
        AdjacentElement = zeros( 3, NElts );
        
        for iElt = 1 : NElts  % Loop over each trianglular element
            
            for iside = 1 : 3 % Loop over triangle sides
                
                if ( AdjacentElement( iside, iElt ) == 0 )
                    
                    Node1 = Node( iCorner( iside, 1 ), iElt );
                    Node2 = Node( iCorner( iside, 2 ), iElt );
                    
                    % Search other elements to find common side
                    for iEltT = 1 : NElts
                        if ( iEltT ~= iElt )
                            Node1T = Node( iCorner( :, 1 ), iEltT );
                            Node2T = Node( iCorner( :, 2 ), iEltT );
                            
                            for isideT = 1 : 3
                                
                                % Do iElt and iEltT share this side?
                                if ( ( Node1 == Node1T( isideT ) && Node2 == Node2T( isideT ) ) || ...
                                        ( Node1 == Node2T( isideT ) && Node2 == Node1T( isideT ) ) )
                                    AdjacentElement( iside,  iElt  ) = iEltT;
                                    AdjacentElement( isideT, iEltT ) = iElt;
                                    break
                                end
                                
                            end % Side2
                        end
                        
                        % did we find the adjacent element? If so, go to next side
                        if ( AdjacentElement( iside, iElt ) ~= 0 )
                            break
                        end
                    end % Element2
                end
                
            end % Side
            
        end % Element1
        
    end

%%
    function [ k, phiR, M, freq, TitleEnv, iSet ] = ReadAllModeSets( rd )
        
        % Reads in the values of the modes at the given receiver depth
        % for every node in the triangulation.
        % Counts up NSets in the process
        
        % The logic is different to the Fortran version because Matlab
        % doesn't have an appropriate CYCLE command
        
        NSets = 0;
        'should use NSets not NNodes; MaxM'
        iSet  = zeros( NNodes, 1 );
        M     = zeros( NNodes, 1 );
        phiR  = zeros( 500, NNodes );
        
        for INode = 1 : NNodes
            
            % Check if the modes have already been read
            if ( INode >= 2 )
                for JNode = 1 : INode - 1
                    if ( strcmp( char( ModeFileName{ INode } ), char( ModeFileName{ JNode } ) ) )
                        iSet( INode ) = iSet( JNode );  % Copy previously read modes
                        break   % we're done for that node
                    end
                end
            end
            
            if ( iSet( INode ) ~= 0 )
                continue   % go do the next node
            end
            
            NSets = NSets + 1;
            iSet( INode ) = NSets;
            
            % Check for 'DUMMY' elts (acoustic absorbers)
            temp = char( ModeFileName{ INode } );
            if ( temp( 1 : 5 ) == 'DUMMY' )
                M( NSets ) = 0;
            else  % Necessary to read in modes
                clear read_modes_bin

                [ Modes ] = read_modes( [ temp '.mod' ], 0.0 );   % 0.0 is a dummy frequency
                
                k(    1 : Modes.M, NSets ) = Modes.k;
                phiR( 1 : Modes.M, NSets ) = interp1( Modes.z, Modes.phi, rd' );   % subsample phi onto the receiver grid
                M( NSets ) = Modes.M;
                freq       = Modes.freqVec( 1 );
                TitleEnv   = Modes.title;
            end
            
        end % NodeLoop
        
    end

%%
    function IElementSource = IdentifySourceElement( xs, ys )
        
        % Identifies the element containing the source at ( xs, ys )
        
        % We define a function Enclosure( xs, ys ) which is 1 when ( xs, ys )
        % is inside a given triangle and decreases from 1 the further the
        % source moves away from the triangle.
        
        % The element with the highest value of Enclosure is identified as the source element.
        % If several elements enclose, the highest numbered element is given posession.
        
        IdentifySourceElement = 0;
        EnclosureMax = 0.0;
        
        for iElt = 1 : NElts
            Node1 = Node( 1, iElt );
            Node2 = Node( 2, iElt );
            Node3 = Node( 3, iElt );
            
            x1 = x( Node1 )   ;   y1 = y( Node1 );
            x2 = x( Node2 )   ;   y2 = y( Node2 );
            x3 = x( Node3 )   ;   y3 = y( Node3 );
            
            % Compute areas of triangles
            Delta = ( x2 * y3 - y2 * x3 ) - ( x1 * y3 - y1 * x3 ) + ( x1 * y2 - y1 * x2 );
            A1    = ( x2 * y3 - y2 * x3 ) - ( xs * y3 - ys * x3 ) + ( xs * y2 - ys * x2 );
            A2    = ( xs * y3 - ys * x3 ) - ( x1 * y3 - y1 * x3 ) + ( x1 * ys - y1 * xs );
            A3    = ( x2 * ys - y2 * xs ) - ( x1 * ys - y1 * xs ) + ( x1 * y2 - y1 * x2 );
            
            Enclosure = abs( Delta ) / ( abs( A1 ) + abs( A2 ) + abs( A3 ) );
            
            if ( Enclosure > EnclosureMax )
                IdentifySourceElement  = iElt;
                EnclosureMax = Enclosure;
            end
            
        end % Element
        
        IElementSource = IdentifySourceElement;
        
        %if ( IdentifySourceElement == 0 ) ERROUT( PRTFile, 'F', 'FIELD3D - IdentifySourceElement', 'Source not inside the grid' )
        
    end
%%
    function [ Outside, RIn, ROut, kIn, kOut, phiIn, phiOut, const, MProp ] = SourceElement( iElement, xs, ys, tsx, tsy, ...
            M, maxM,k, phiR, phiS )
        
        % Given an element number for the source
        % Computes modal information
        
        MProp = intmax;
        
        % Coordinates of the centroid of the source element
        
        xCenter = sum( x( Node( 1 : 3, iElement ) ) ) / 3.0;
        yCenter = sum( y( Node( 1 : 3, iElement ) ) ) / 3.0;
        
        for ISide = 1 : 3
            
            Node1 = Node( iCorner( ISide, 1 ), iElement );
            Node2 = Node( iCorner( ISide, 2 ), iElement );
            iSet1 = iSet( Node1 );
            iSet2 = iSet( Node2 );
            MProp = min( [ MProp, M( iSet1 ), M( iSet2 ) ] );
            
            x1S = x( Node1 ) - xs;
            y1S = y( Node1 ) - ys;
            txC = x( Node1 ) - xCenter;
            tyC = y( Node1 ) - yCenter;
            tx  = x( Node2 ) - x( Node1 );
            ty  = y( Node2 ) - y( Node1 );
            
            Delta = tsx * ty - tsy * tx;
            
            % *** Radial parallel to side? ***
            
            if ( Delta == 0.0 )
                SV(  ISide ) = realmax;
            else
                RVC( ISide ) = ( txC * ty  - tyC * tx  ) / Delta;
                RV(  ISide ) = ( x1S * ty  - y1S * tx  ) / Delta;
                SV(  ISide ) = ( x1S * tsy - y1S * tsx ) / Delta;
            end
        end
        
        % Identify two good sides and one bad side based on the intercept point
        
        IBad = 1;
        if ( abs( SV( 2 ) - 0.5 ) > abs( SV( IBad ) - 0.5 ) )
            IBad = 2;
        end
        
        if ( abs( SV( 3 ) - 0.5 ) > abs( SV( IBad ) - 0.5 ) )
            IBad = 3;
        end
        
        IGood1 = 1;
        IGood2 = 2;
        if ( IBad == 1 )
            IGood1 = 3;
        else
            if ( IBad == 2 )
                IGood2 = 3;
            end
        end
        
        % The side with the lesser RVC is the inside
        
        if ( RVC( IGood1 ) < RVC( IGood2 ) )
            Inside  = IGood1;
            Outside = IGood2;
        else
            Inside  = IGood2;
            Outside = IGood1;
        end
        
        SIn      = SV( Inside );
        SIn      = min( max( SIn,  0.0 ), 1.0 );
        RIn      = RV( Inside );
        
        SOut     = SV( Outside );
        SOut     = min( max( SOut, 0.0 ), 1.0 );
        ROut     = RV( Outside );
        
        % Get values of modes at Z = SD and (x, y) = intercept points
        
        iCorner1 = iCorner(  Inside, 1 );
        iCorner2 = iCorner(  Inside, 2 );
        iCorner3 = iCorner( Outside, 1 );
        iCorner4 = iCorner( Outside, 2 );
        
        % Interpolate to get modal values at source
        R = 0.0;
        
        if ( RIn ~= ROut )
            alpha = ( R - RIn ) / ( ROut - RIn );
            alpha = min( max( alpha,  0.0 ), 1.0 );
        else
            alpha = 0.0;
        end
        
        L = 1 : MProp;
        phiIn  = phiS( L, iCorner1 ) + SIn   * ( phiS( L, iCorner2 ) - phiS( L, iCorner1 ) );
        phiOut = phiS( L, iCorner3 ) + SOut  * ( phiS( L, iCorner4 ) - phiS( L, iCorner3 ) );
        const  = phiIn( L )          + alpha * ( phiOut( L )         - phiIn( L ) );
        
        % Obtain values of modes at z = Rd and (x, y)=intercept points
        
        [ kIn,  phiIn  ] = InterpolateModes( iElement, Inside,  SIn,  M, k, phiR );
        [ kOut, phiOut ] = InterpolateModes( iElement, Outside, SOut, M, k, phiR );
        
    end
%%

    function [ Outside, ROut, kOut, phiOut ] = OUT( iElement, NewElement, xs, ys, tsx, tsy, M, k, phiR )
        % Given an element number and a side through which prop path enters
        % computes mode info at exit point
        
        % If no adj elt then exit. Previous values of KOut, phiOut retained
        
        if ( NewElement == 0 )
            ROut = realmax;
            return
        end
        
        % *** Loop over sides to find outside ***
        
        sOut = realmax;
        
        for ISide = 1 : 3
            
            % Don't try an outside the same as the inside
            if ( AdjacentElement( ISide, NewElement ) ~= iElement )
                
                Node1 = Node( iCorner( ISide, 1 ), NewElement );
                Node2 = Node( iCorner( ISide, 2 ), NewElement );
                
                x1S = x( Node1 ) - xs;
                y1S = y( Node1 ) - ys;
                Tx  = x( Node2 ) - x( Node1 );
                Ty  = y( Node2 ) - y( Node1 );
                
                Delta = tsx * Ty - tsy * Tx;
                
                % *** Radial parallel to side? ***
                
                if ( Delta == 0.0 )
                    sT    = realmax;
                else
                    ROutT = ( x1S * Ty  - y1S * Tx  ) / Delta;
                    sT    = ( x1S * tsy - y1S * tsx ) / Delta;
                end
                
                % If intercept is in segment, store the side number
                % This logic takes the side crossing that happens closest to
                % the middle of the side
                if ( abs( sT - 0.5 ) < abs( sOut - 0.5) )
                    Outside = ISide;
                    sOut    = sT;
                    ROut    = ROutT;
                end
            end
        end
        
        % Obtain values of modes at intercept points
        
        [ kOut, phiOut ] = InterpolateModes( NewElement, Outside, sOut, M, k, phiR );
        
    end

%%

    function [ kInt, phiInt ] = InterpolateModes( iElement, ISide, s, M, k, phiR )
        
        % Given an element, a side, and a proportional distance along the side
        % Returns interpolated modal values
        
        Node1 = Node( iCorner( ISide, 1 ), iElement );
        Node2 = Node( iCorner( ISide, 2 ), iElement );
        iSet1 = iSet( Node1 );
        iSet2 = iSet( Node2 );
        MProp = min( [ MProp, M( iSet1 ), M( iSet2 ) ] );
        
        st = min( max( s, 0.0 ), 1.0 ); % Extrapolation blocked by making sure s in [0, 1]
        
        ii = 1 : MProp;
        kInt   =    k( ii, iSet1 ) + st * (    k( ii, iSet2 ) -    k( ii, iSet1 ) );
        phiInt = phiR( ii, iSet1 ) + st * ( phiR( ii, iSet2 ) - phiR( ii, iSet1 ) );
        phiInt = phiInt.';   % make it a row vector
    end

%%

    function P = EvaluatePDQ( freq, k, phiR, phiS, M, maxM, IelementSource, xs, ys, ...
            theta, Ntheta, RminM, RmaxM, Nr, MLimit, AdjacentElement )
        
        % Computes 3-D pressure field using adiabatic mode theory.
        % Normalized to pressure of point source at 1 meter.
        % Note RminM must be zero.
        
        P       = zeros( Ntheta, Nr );    % zero out the pressure field
        delta_r = ( RmaxM - RminM ) / ( Nr - 1 );
        
        %  *** Loop over angle ***
        
        for Itheta = 1 : Ntheta
            tsx = cosd( theta( Itheta ) );
            tsy = sind( theta( Itheta ) );
            
            % Get modal values
            iElement = IelementSource;
            
            [ Outside, RIn, ROut, kIn, kOut, phiIn, phiOut, const, MProp ] = SourceElement( iElement, xs, ys, tsx, tsy, M, maxM, ...
                k, phiR, phiS );
            
            MProp = min( MLimit, MProp );
            
            if ( MProp > 0 )   % any propagating modes at the source?
                
                % apply the source beam pattern
                % using kIn for wavenumber; would be more precise to use interpolated value at source
                %             if ( SBPFlag == '*' && Itheta == 1 )
                %                kz2    = zeros( MProp, 1 );
                %                thetaT = zeros( MProp, 1 );
                %                S      = zeros( MProp, 1 );
                %
                %                c0    = 1500;   % reference sound speed, should be speed at the source depth
                %                omega = 2 * pi * freq;
                %                kz2   = real( omega ^ 2 / c0 ^ 2 - kIn( 1 : MProp ) ^ 2 );      % vertical wavenumber squared
                %                kz2( kz2 < 0 ) = 0;                                         % remove negative values
                %
                %                thetaT = RadDeg * atan( sqrt( kz2 ) / real( kIn( 1 : MProp ) ) );  % calculate the angle in degrees
                %                S = interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT );
                %                const( 1 : MProp ) = const( 1 : MProp ) * S;                       % apply the shading
                %             end
                
                const    = 1i * sqrt( 2.0 * pi ) * exp( 1i * pi / 4.0 ) * const( 1 : MProp ).';
                kAverage = 0.5 * ( kIn( 1 : MProp ) + kOut( 1 : MProp ) );
                PhaseInc = exp( -1i * kAverage( 1 : MProp ) * delta_r ).';
                phiIn    = phiIn( 1 : MProp );   % reduce size of phiIn to the number allowed to propagate
                
                for ir = 1 : Nr % March forward in range
                    RM = RminM + ( ir - 1 ) * delta_r;
                    if ( RM == 0.0 )
                        RM = min( 1.0, delta_r );
                    end
                    
                    % Crossing into new element?
                    while ( RM > ROut )
                        
                        % Copy outside info to inside
                        NewElement = AdjacentElement( Outside, iElement );
                        if ( NewElement == 0 )
                            break   % jump out if there is no adjacent element
                        end
                        
                        % RIn                   = ROut;
                        kIn      = kOut(   1 : MProp );
                        phiIn    = phiOut( 1 : MProp );
                        kAverage = 0.5 * ( kIn + kOut( 1 : MProp ) );
                        PhaseInc = exp( -1i * kAverage * delta_r ).';
                        
                        % Get new outside info
                        [ Outside, ROut, kOut, phiOut ] = OUT( iElement, NewElement, xs, ys, tsx, tsy, M, k, phiR );
                        
                        if ( MProp <= 0 )
                            break   % jump out if there are no propagating modes
                        end
                        iElement = NewElement;
                        
                    end
                    
                    % Compute modal contribution at this range
                    
                    const              = const .* PhaseInc;   % advance the phase
                    P( Itheta, ir )    = ( const * phiIn.' ) / sqrt( RM * kIn( 1 ) );
                end   % RangeStepping
            end
        end   % Bearing
        
    end
end






