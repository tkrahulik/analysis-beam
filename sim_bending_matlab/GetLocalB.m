% fuction returning magnetic field B at current position R
function [ BLocal ] = GetLocalB( BMap, R )

if ( strcmp( BMap , 'BSimple1' ) )
  % simple model: use z position and if-statement to get B field
    if ( R(3) < 0.5 )
        BLocal = [ 0 0 0 ];
    elseif ( R(3) > 0.75 ) && ( R(3) < 1.25 )
        BLocal = [ -0.001 0 0 ];
    else
        BLocal = [ 0 0 0 ];
    end

elseif ( strcmp( BMap , 'BGauss5' ) )
    gmean = 1    % in meter
    gsigma = 0.1 % in meter
    gpeak = 5    % in mT
    
    Bx = gpeak / 4000. * Gauss( R(3) , gmean , gsigma )
    
    BLocal = [ Bx 0 0 ]

else
    BLocal = [ 0 0 0 ];
end
    
end