function [g] = evalg(x,n)

global a b c d problem detx A sumlog sqrtx w

% ---------------  Problem 1 ---------------

if ( problem == 1 )
    for i = 1:n
        g(i) = a(i) * b(i) * exp( - b(i) * x(i) ) + 2 * c(i) * log( x(i) ) / x(i) + d(i) / x(i);
    end
    return
end

% ---------------  Problem 2 ---------------

if ( problem == 2 )
    for i = 1:n
        g(i) = a(i) * d(i) * x(i)^( d(i) - 1 )/( x(i)^d(i) + b(i) ) - c(i) / x(i);
    end
    return
end

% ---------------  Problem 3 ---------------

if ( problem == 3 )
    detx = det(x);
    
    g = ( 2 * a * log( detx ) - b ) * inv(x);
    
    if ( ~isreal(g) )
        g = real(g);
    end
    
    if ( ~issymmetric(g) )
        g = ( g + g' )/2;
    end
    return
end

% ---------------  Problem 4 ---------------

if ( problem == 4 )
    detx = det(x);
    g = ( a * d * detx^d / ( detx^d + b ) - c ) * inv(x);
    
    if ( ~isreal(g) )
        g = real(g);
    end
    
    if ( ~issymmetric(g) )
        g = ( g + g' )/2;
    end
    return
end

% ---------------  Problem 5 ---------------

% Center of mass on the positive orthant
if ( problem == 5 )
    for i = 1:n
        g(i) = sum( log( x(i) ./ w(i,:) ) ) / x(i);
    end 
    return
end

% ---------------  Problem 6 ---------------

% Center of mass on the SDP matrices cone
if ( problem == 6 )
    
    sqrtx    = sqrtm( x );
    invsqrtx = inv( sqrtx );
    
    sumlog = 0;
    for i = 1:size(A,3)
        sumlog = sumlog + logm( sqrtx * inv( A(:,:,i) ) * sqrtx );
    end
    
    g = invsqrtx * sumlog * invsqrtx;
    
    if ( ~isreal(g) )
        g = real(g);
    end
    
    if ( ~issymmetric(g) )
        g = ( g + g' )/2;
    end   
    return
end