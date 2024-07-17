clear all
clc

global a b c d problem A w

% Set the Strategy for the Riemannian gradient method
% Strategy = 1: Lipschitz stepsize
% Strategy = 2: Adaptive stepsize
% Strategy = 3: Armijo's stepsize

strategy = 2;

% Set the problem to be solved:
% problem = 1: Problem 1 of the paper
% problem = 2: Problem 2 of the paper
% problem = 3: Problem 3 of the paper
% problem = 4: Problem 4 of the paper
% problem = 5: Center of mass on the positive orthant
% problem = 6: Center of mass on the SDP matrices cone

problem = 6;

if ( problem == 1 ) 
    
    % Set the dimension of the problem
    n = 100;
    
    % Set the problem data
    alpha = 0;
    beta  = 10;
    
    a(1,1:n) = alpha + rand * ( beta - alpha );
    b(1,1:n) = alpha + rand * ( beta - alpha );
    c        = 1.1 * a + rand * 3.9 * a;
    d(1,1:n) = alpha + rand * ( beta - alpha );
    
    % Set the stating point
    l = 0;
    u = 20;
    
    xini = l + rand(1,n) * ( u - l );
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);
    
elseif ( problem == 2 ) 
    
    % Set the dimension of the problem
    n = 100;
    
    % Set the problem data
    alpha = 0;
    beta  = 10;
    
    a(1,1:n) = alpha + rand * ( beta - alpha );
    b(1,1:n) = alpha + rand * ( beta - alpha );

    alpha = 2;
    beta  = 10;

    d(1,1:n) = alpha + rand * ( beta - alpha );
    c = rand * a .* d;
    
    % Set the stating point
    l = 0;
    u = 20;
    
    xini = l + rand(1,n) * ( u - l );
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);
    
elseif ( problem == 3 )
    
    % Set the dimension of the problem
    n = 5;  
    
    % Set the problem data
    alpha = 0;
    beta  = 10;

    a = alpha + rand * ( beta - alpha );
    b = alpha + rand * ( beta - alpha );
    c = [];
    d = [];
    
    % Set the stating point
    l = 0;
    u = 20;
    
    lambda = l + rand(1,n) * ( u - l );

    xini = orth(randn(n));
    xini = xini * diag(lambda) * xini';

    if ( ~isreal(xini) )
        xini = real(xini);
    end

    if ( ~issymmetric(xini) )
        xini = ( xini + xini' )/2;
    end   
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);
    
elseif ( problem == 4 )
    
    % Set the dimension of the problem
    n = 5;  
    
    % Set the problem data
    alpha = 0;
    beta  = 10;
    
    a = alpha + rand * ( beta - alpha );
    b = alpha + rand * ( beta - alpha );
    d = alpha + rand * ( beta - alpha );
    c = rand * a * d;

    % Set the stating point
    l = 0;
    u = 20;
    
    lambda = l + rand(1,n) * ( u - l );

    xini = orth(randn(n));
    xini = xini * diag(lambda) * xini';

    if ( ~isreal(xini) )
        xini = real(xini);
    end

    if ( ~issymmetric(xini) )
        xini = ( xini + xini' )/2;
    end   
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);
    
elseif ( problem == 5 )
    
    % Set the dimension of the problem
    n = 20;
    
    % Set the number of points
    m = 5;
    
    % Set the problem data
    l = 0;
    u = 100;

    for j = 1:m
        w(:,j) = l + rand(1,n) * ( u - l );
    end

    % Set the stating point
    xini = l + rand(1,n) * ( u - l );
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);

elseif ( problem == 6 )
    
    % Set the dimension of the matrices
    n = 100;
    
    % Set the number of matrices
    m = 20;
    
    % Set the problem data
    l = 0;
    u = 100;
    
    lambda = l + rand(1,n) * ( u - l );
    
    for i = 1:m

        A(:,:,i) = orth(randn(n));
        A(:,:,i) = A(:,:,i) * diag(lambda) * A(:,:,i)';

        if ( ~isreal(A(:,:,i)) )
            A(:,:,i) = real(A(:,:,i));
        end

        if ( ~issymmetric(A(:,:,i)) )
            A(:,:,i) = ( A(:,:,i) + A(:,:,i)' )/2;
        end      
        
    end
    
    % Set the stating point
    xini = zeros(n);
    for i = 1:m
        xini = xini + logm( A(:,:,i) );
    end
    xini = expm( xini / m );
    
    % Call the methods
    [x,solinfo,it,time,nfev,ngev] = GradRiemannian(xini,n,strategy);
    %[x,solinfo,it,time,nfev,ngev] = Gradient(xini,n);
    
end  