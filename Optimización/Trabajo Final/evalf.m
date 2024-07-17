function [f] = evalf(x,n)

global a b c d problem A w

% ---------------  Problem 1 ---------------

if ( problem == 1 )
    f = 0;
    for i = 1:n
        f = f - a(i) * exp( - b(i) * x(i) ) + c(i) * log( x(i) )^2 + d(i) * log( x(i) );
    end
    return
end

% ---------------  Problem 2 ---------------

if ( problem == 2 )
    f = 0;
    for i = 1:n
        f = f + a(i) * log( x(i) ^ d(i) + b(i) ) - c(i) * log( x(i) );
    end
    return
end

% ---------------  Problem 3 ---------------

if ( problem == 3 )
    detx = det( x );
   
    f = a * log( detx )^2 - b * log( detx );
    return
end

% ---------------  Problem 4 ---------------

if ( problem == 4 )
    detx = det( x );
    f = a * log( detx^d + b ) - c * log( detx );
    return
end

% ---------------  Problem 5 ---------------

% Center of mass on the positive orthant
if ( problem == 5 )
    f = 0;
    for j = 1:size(w,2)
        f = f + sum( log( w(:,j)' ./ x ).^2 );
    end
    f = f / 2;
    return  
end

% ---------------  Problem 6 ---------------

% Center of mass on the SDP matrices cone
if ( problem == 6 )
       
    invsqrtx = inv( sqrtm(x) );
    
    f = 0;
    for i = 1:size(A,3)
        f = f + 1/2 * ( norm( logm( invsqrtx * A(:,:,i) * invsqrtx ) ,'fro') )^2;
    end
    return  
end