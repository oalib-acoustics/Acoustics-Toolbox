function h = reducestep( SSP, x0, Tray, Layer0, c, xTop, nTop, xBot, nBot, rTopSeg, rBotSeg, deltas, h )

% Reduces the ray step size to make sure we land on interfaces and boundaries

cTray = c * Tray;
x = x0 + h * cTray;  % make a trial step

% *** Detect interface or boundary crossing
% *** Adjust step to land on an interface or boundary

hLayer = inf;

%if( abs( cTray( 2 ) ) > eps( hLayer ) )
if ( x( 2 ) < SSP.z( Layer0 ) )
    hLayer = ( SSP.z( Layer0 ) - x0( 2 ) ) / cTray( 2 );
end
%end

Layer1 = Layer0 + 1;

if ( x( 2 ) > SSP.z( Layer1 ) )
    hLayer = ( SSP.z( Layer1 ) - x0( 2 ) ) / cTray( 2 );
end

% top crossing
hTop = inf;
d    = x - xTop;       % vector from top    to ray endpoint

if ( d * nTop' > 0 )   % top crossing
    e    =  x0 - xTop;     % vector from top    to ray origin
    hTop = -e * nTop' / ( cTray * nTop' );
end

% bottom crossing
hBot = inf;
d    = x - xBot;       % vector from bottom to ray endpoint

if ( d * nBot' > 0 )   % bottom crossing
    e    =  x0 - xBot;     % vector from bottom to ray origin
    hBot = -e * nBot' / ( cTray * nBot' );
end

% top segment crossing in range
h4 = inf;

if ( x( 1 ) < rTopSeg( 1 ) && cTray( 1 ) ~= 0 )
    h4 = -( x0( 1 ) - rTopSeg( 1 ) ) / cTray( 1 );
end
if ( x( 1 ) > rTopSeg( 2 ) && cTray( 1 ) ~= 0 )
    h4 = -( x0( 1 ) - rTopSeg( 2 ) ) / cTray( 1 );
end


% bottom segment crossing in range
h5 = inf;

if ( x( 1 ) < rBotSeg( 1 ) && cTray( 1 ) ~= 0 )
    h5 = -( x0( 1 ) - rBotSeg( 1 ) ) / cTray( 1 );
end
if ( x( 1 ) > rBotSeg( 2 ) && cTray( 1 ) ~= 0 )
    h5 = -( x0( 1 ) - rBotSeg( 2 ) ) / cTray( 1 );
end

h = min( [ h, hLayer, hTop, hBot, h4, h5 ] );   % step is limited by the smaller of an interface or bottom crossing
%h = max( h, 1e4 * eps( deltas ) )   % make sure we make some motion

if ( h < 1.0d-4 * deltas )   % is it taking an infinitesimal step?
     h = 1.0d-5 * deltas;    % make sure we make some motion
     %iSmallStepCtr = iSmallStepCtr + 1;   % keep a count of the number of sequential small steps
else
     %iSmallStepCtr = 0;   % didn't do a small step so reset the counter
end