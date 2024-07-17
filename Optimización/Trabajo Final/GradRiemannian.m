function [x,solinfo,it,time,nfev,ngev] = GradRiemannian(x,n,strategy) 

global a b c d problem detx A sumlog sqrtx w

% -------------------------------------------------------------------------
%   Initialization
% -------------------------------------------------------------------------

% Start timing

tic;

% Set default parameters

opttol = 10^(-5);
maxit  = 1000;
ftol   = 1/2;
smallnum = 10^(-20);
iprint = true;

% Define the constant stepsize for strategy 1

if ( strategy == 1 )
    if ( problem == 1 )
        t = 1 / sum( ( a + 2*c ).^2 );
    elseif ( problem == 2 )
        t = 1 / sum( a.^2 .* d.^4 );
    elseif ( problem == 3 )
        t = 1 / ( 2 * a * sqrt(n) );
    elseif ( problem == 4 )
        t = 1 / ( a * d^2 * n );    
    elseif ( problem == 5 )
        radius = 0;
        for i = 1:size(w,2)
            radius = max( radius, norm( log( w(:,i)' ./ x ) ) );
        end
        t = 1 / ( 4 * radius * coth( 4 * radius ) );
        t = 1.99 * t;
    elseif ( problem == 6 )
        invsqrtx = inv( sqrtm(x) );    

        radius = 0;
        for i = 1:size(A,3)
            radius = max( radius, norm( logm( invsqrtx * A(:,:,i) * invsqrtx ) ,'fro') );
        end
        t = 1 / ( 4 * radius * coth( 4 * radius ) );
        t = 1.99 * t;
    end
end

% Set the initial trial value for the Lipschitz constant for strategy 3

if ( strategy == 3 )
    L0 = 1;
    t  = 1 / L0;
end
    
% Counters

it   = 0;
nfev = 0;
ngev = 0;

% Evaluate objective function and gradient

[f]  = evalf(x,n);
nfev = nfev + 1;

[g]  = evalg(x,n);
ngev = ngev + 1;

% Compute the sup-norm of the gradient

gsupn = norm(g,inf);

% -------------------------------------------------------------------------
%   Main loop
% -------------------------------------------------------------------------
    
while(1)
    
    % Print information
    
    if ( iprint ) 
        if ( it == 0 )
            fprintf('----------------------------------------------------------------\n')
            fprintf('                  Riemannian gradient method                    \n')
            fprintf('----------------------------------------------------------------\n')
        end
        
        if ( mod(it,10) == 0 )
            fprintf('\n')
            fprintf('  It    objective     gsupn          t           nfev \n')
        end
        
        if ( it == 0 )
            fprintf(' %4d   %+6.2e    %6.2e        -         %5d\n',it,f,gsupn,nfev)
        else
            fprintf(' %4d   %+6.2e    %6.2e     %6.2e     %5d\n',it,f,gsupn,t,nfev)
        end
    end
    
    % ---------------------------------------------------------------------
    %   Test stopping criteria
    % ---------------------------------------------------------------------
    
    % Test optimality		
    
    if ( gsupn <= opttol )
        solinfo = 0;
        
        % Stop timing
        
        time = toc;
        
        if ( iprint )
            fprintf('\n')
            fprintf('Solution was found \n\n')
            fprintf('Number of functions evaluations  : %6d \n',nfev)
			fprintf('Total CPU time in seconds: %10.4f \n\n',time)		
        end
                            
        return
    end
    
    % Test whether the number of iterations is exhausted
    
    if ( it >= maxit )
        solinfo = 1;
        
        % Stop timing
        
        time = toc;
        
        if ( iprint )
            fprintf('\n')
            fprintf('Maximum of iterations reached \n\n')
            fprintf('Number of functions evaluations  : %6d \n',nfev)
			fprintf('Total CPU time in seconds: %10.4f \n\n',time)		
        end
        
        return
    end
    
    % Increment the iteration counter
    
    it = it + 1;
    
    % ---------------------------------------------------------------------
    %   Use the constant stepsize for strategy 1 or perform the line search
    %   for strategys 2 and 3
    % ---------------------------------------------------------------------
    
    % Initialize the trial stepsize for strategy 2. For strategy 3 simply
    
    if ( strategy == 2 )
        t = 1;
    end
    
    % Define the trial point
    
    if ( problem == 1 || problem == 2 || problem == 5 )
        xgrad = x .* g;
        gradfeucn2 = dot(xgrad,xgrad);
        xtrial   = x .* exp( - t * xgrad );
        xtrial = max( xtrial , smallnum );
    elseif ( problem == 3 )
        gradfeucn2 = trace( ( x * g )^2 );
        %detx = det(x);
        xtrial = ( detx^( 2 * a ) * exp( - b ) )^( -t ) * x;   
    elseif ( problem == 4 )
        gradfeucn2 = trace( ( x * g )^2 );
        %detx = det(x);
        xtrial = exp( - t * ( ( a * d * detx^d / ( detx^d + b ) - c ) ) ) * x;       
    elseif ( problem == 6 )    
        gradfeucn2 = trace( ( x * g )^2 );        
        xtrial = sqrtx * expm( - t * sumlog ) * sqrtx;
        % See evafg for sqrtx and sumlog
    end  
    
    % Evaluate the objective function at the trial point
        
    [ftrial] = evalf(xtrial,n);
    nfev = nfev + 1;
    
    % ---------------------------------------------------------------------
    %   Line search
    % ---------------------------------------------------------------------
    
    % Perform the line search for strategys 2 and 3

    if ( strategy == 2 || strategy == 3 )
        while ( ftrial > f - ftol * t * gradfeucn2 )
            
            % Define the new trial stepsize
            
            t = t/2;
            
            % Define the new trial point

            if ( problem == 1 || problem == 2 || problem == 5 )
                xtrial   = x .* exp( - t * xgrad );
                xtrial = max( xtrial , smallnum );
            elseif ( problem == 3 )
                xtrial = ( detx^( 2 * a ) * exp( - b ) )^( -t ) * x;
            elseif ( problem == 4 )
                xtrial = exp( - t * ( ( a * d * detx^d / ( detx^d + b ) - c ) ) ) * x;
            elseif ( problem == 6 )
                xtrial = sqrtx * expm( - t * sumlog ) * sqrtx;
            end
            
            % Evaluate the objective function at the new trial point

            [ftrial] = evalf(xtrial,n);
            nfev = nfev + 1;
        end
    end
    
    % ---------------------------------------------------------------------
    %   Update x
    % ---------------------------------------------------------------------
    
    x = xtrial;
    f = ftrial;
        
    % Compute the gradient of f

    [g]  = evalg(x,n);
    ngev = ngev + 1;
    
    % Compute the sup-norm of the gradient

    gsupn = norm(g,inf);
    
    % ---------------------------------------------------------------------
    %   Iterate
    % ---------------------------------------------------------------------
end