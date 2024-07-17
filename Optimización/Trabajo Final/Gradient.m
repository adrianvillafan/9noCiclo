function [x,solinfo,it,time,nfev,ngev] = Gradient(x,n) 

global problem

% -------------------------------------------------------------------------
%   Initialization
% -------------------------------------------------------------------------

% Start timing

tic;

% Set default parameters

opttol = 10^(-5);
maxit  = 1000;
ftol   = 1/2;
eps    = 10^(-10);
iprint = true;

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
            fprintf('                  Euclidian gradient method                     \n')
            fprintf('----------------------------------------------------------------\n')
        end
        
        if ( mod(it,10) == 0 )
            fprintf('\n')
            fprintf('  It    objective     gsupn          t           nfev \n')
        end
        
        if ( it == 0 )        
            fprintf(' %4d   %+6.2e    %6.2e        -         %5d \n',it,f,gsupn,nfev)
        else
            fprintf(' %4d   %+6.2e    %6.2e     %6.2e     %5d \n',it,f,gsupn,t,nfev)
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
    %   Line search
    % ---------------------------------------------------------------------
    
    if ( problem == 1 || problem == 2 || problem == 5 )
        
        % Compute ||g||^2
        
        geucn2 = dot(g,g);
        
        % Compute the maximum value of t such that t <= 1 and x - t * g > 0
        
        xg = x ./ g;
        xg = xg( xg > 0 );
        
        if ( isempty(xg) )
            t = 1;
        else
            t = min( xg ) - eps;
            t = min( 1 , t ); 
        end
        
    elseif ( problem == 3 || problem == 4 || problem == 6 )
        
        % Compute ||g||^2
        
        geucn2 = trace( g ^ 2 );
        
        % Compute the maximum value of t such that t <= 1 and x - t * g > 0
        
        % Compute the largest eigvalue of g
        
        eigg = eigs(g,1,'LA');
        
        if ( eigg <= 0 )
            t = 1;
        else
            
            % Compute the smallest eigvalue of x
            
            eigx = eigs(x,1,'SA');
            
            if ( eigx <= 0 )
                
                solinfo = 2;
        
                % Stop timing
        
                time = toc;
        
                if ( iprint )
                    fprintf('\n')
                    fprintf('Error: the currente iterate is not positive definite \n\n')
                    fprintf('Number of functions evaluations  : %6d \n',nfev)
                    fprintf('Total CPU time in seconds: %10.4f \n\n',time)		
                end

                return
            end

            t = eigx / eigg;
            
            if ( t <= eps )
                t = 0.1 * t;
                disp('Warning: Allowable stepsize too small')
            else
                t = t - eps;  
            end
            
            t = min( 1 , t );
        end
         
    end   
    
    % Define the trial point
       
    xtrial   = x - t * g;
    
    % Evaluate the objective function at the trial point

    [ftrial] = evalf(xtrial,n);
    nfev = nfev + 1;

    while ( ftrial > f - ftol * t * geucn2 )
        
        % Define the new trial stepsize
        
        t = t/2;
        
         % Define the new trial point
         
        xtrial   = x - t * g;
        
        % Evaluate the objective function at the new trial point
        
        [ftrial] = evalf(xtrial,n);
        nfev = nfev + 1;
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