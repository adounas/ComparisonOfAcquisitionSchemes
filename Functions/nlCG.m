function [ x, f, relResVec ] = nlCG( y , A, Psi, varargin )
%//////////////////////////////////////////////////////////////////////////
% Non-linear conjugate gradient (CG) algorithm minimizing the following
% objective function:
%
% f(x) = ||A*x - y||^2_2 + lambda*||Psi*x||_1.
%
% References:  
% - Lustig M et al., "Sparse MRI: The Application of Compressed Sensing for
%   Rapid MR Imaging", MRM 58(6):1182-95, 2007.
% - SparseMRI V0.2 Matlab toolbox, 
%   http://www.eecs.berkeley.edu/~mlustig/Software.html
%
% Regularization weights of multiple operators composing A and Psi must be 
% passed within matrix/class object/function handle A and Psi. 
%
% Inputs:             
% -------
% y:      Measurement vector.
% A:      L2-term system matrix.
% Psi:    Sparsifying transform matrix.
% varargin{1-2}: - x0: Initial estimate.
%                - params: Reconstruction parameter struct.
%                  - lambda: L1-term penalty weight.
%                  - nIter: Maximum no. of main loop iterations.
%                  - nIterLS: Maximum no. of line search iterations.  
%                  - gradTol: Gradient tolerance.
%                  - l1Smooth: Smoothing parameter for absolute value  
%                              derivative approximation.   
%                  - alpha, beta, t0: Line search parameters.
% 
% Outputs:    
% --------
% x:            Reconstructed vector.
% f:            Objective function values.
% relResVec:    Relative residual vector (from L2-term).
%
% Function calls:    processInput, preObjective, objective, gradObjective.
% ---------------
%
% Claudio Santelli, 2015
%//////////////////////////////////////////////////////////////////////////

%==========================================================================
% Process input.
%==========================================================================
narginchk(3,5);
[x0, params] = processInput( y, A, Psi, varargin );
%--------------------------------------------------------------------------
% Get parameters.
%--------------------------------------------------------------------------
lambda   = params.lambda;
nIter    = params.nIter;
nIterLS  = params.nIterLS;
gradTol  = params.gradTol;
l1Smooth = params.l1Smooth;
alpha    = params.alpha;
beta     = params.beta;
t0       = params.t0;

%==========================================================================
% CG iterations.
%==========================================================================
x = x0; clear x0
if nIter>0
    %----------------------------------------------------------------------
    % Init.
    %----------------------------------------------------------------------
    g0 = gradObjective( x, A, Psi, y, lambda, l1Smooth );
    dx = -g0;
    for k=1:nIter
        if norm(g0(:))<gradTol, break; end % vs norm(dx(:))<gradToll
        %------------------------------------------------------------------
        % Backtracking line-search.
        %------------------------------------------------------------------
        t = t0;
        [ A_x, A_dx, Psi_x, Psi_dx ] = preObjective( x, dx, A, Psi, lambda );
        f0                   = objective( A_x, A_dx, Psi_x, Psi_dx, y, lambda, l1Smooth, 0 );
        [f1, ~, relRes] = objective( A_x, A_dx, Psi_x, Psi_dx, y, lambda, l1Smooth, t );

        i = 0;
        while ( f1 > f0 - alpha*t*abs(g0(:)'*dx(:)) ) && i<nIterLS  % vs. "... +alpha*t*real(..."
            i = i+1;
            t = beta*t; 
            [f1, ~, relRes] = objective( A_x, A_dx, Psi_x, Psi_dx, y, lambda, l1Smooth, t );
        end
        if i==nIterLS
            disp('Maximum number of line search iterations reached!');
            return;
        end

        %------------------------------------------------------------------
        % If necessary, adapt the initial step size, t0, to control the 
        % number line search iterations.
        %------------------------------------------------------------------
        if i>2, t0 = beta*t0; end
        if i<0, t0 = t0/beta; end

        f(k) = f1;
        relResVec(k) = relRes;
%         fprintf( 'NLCG iteration no. ... %d   , f: %f, rel. Res.: %f, L-S: %d \n', k, f1, relRes, i );

        %------------------------------------------------------------------
        % Update vectors
        %------------------------------------------------------------------
        x = x + t*dx;
        if k<nIter
            g     = gradObjective( x, A, Psi, y, lambda, l1Smooth ); % Gradient.
            gamma = (g(:)'*g(:)) / (g0(:)'*g0(:)+eps);
            dx    = -g + gamma*dx;                                   % Search direction.
            g0    = g;
        end
    end
end

end

function [x0, params] = processInput( y, A, Psi, argInArray ) 
%//////////////////////////////////////////////////////////////////////////
% Checks and, if neccessary, intitializes input and reconstruction
% parameters.
%
% argInArray{1,2} = {x0, params}.
%//////////////////////////////////////////////////////////////////////////

%==========================================================================
% y.
%==========================================================================
if isempty(y) || ~isnumeric(y) 
    error('Measurement vector y must be a numeric array!');
end
%==========================================================================
% x0, params.
%==========================================================================
x0 = []; params = [];
if numel(argInArray)==1
    if isstruct(argInArray{1})
        params = argInArray{1};
    elseif isnumeric(argInArray{1}) 
        x0 = argInArray{1};
    elseif ~isempty(argInArray{1})
        error('Fourth input argument must be a numeric array (x0) or a struct (params)!');
    end
elseif numel(argInArray)==2
    if isstruct(argInArray{1})
        params = argInArray{1};
        if isnumeric(argInArray{2})
            x0 = argInArray{2};
        elseif ~isempty(argInArray{2})
            error('Fifth input argument must be a numeric array (x0)!');
        end
    elseif isnumeric(argInArray{1})
        x0 = argInArray{1};
        if isstruct(argInArray{2})
            params = argInArray{2};
        elseif ~isempty(argInArray{2})
            error('Fifth input argument must be a struct (params)!');
        end
    elseif ~isempty(argInArray{1})
        error('Fourth input argument must be a numeric array (x0) or a struct (params)!');
    end
end
if isempty(x0)
    if isa(A,'function_handle')
        x0 = zeros( size( A(y,'transp') ), class(y) );
    else
        x0 = zeros( size( A'*y ), class(y) );
    end
end
if isempty(params), params = struct(); end 
%==========================================================================
% A.
%==========================================================================
if isempty(A)
    params.A = 1;
    warning('');
end
%==========================================================================
% Psi.
%==========================================================================
if isempty(Psi)
    params.Psi = 0;
    params.lambda = 0;
    warning('');
end
%==========================================================================
% params.
%==========================================================================
%--------------------------------------------------------------------------
% lambda.
%--------------------------------------------------------------------------
if ~isfield(params,'lambda') || isempty(params.lambda)
    params.lambda = 0;
%     warning('Regularization parameter lambda equals to 0. Reconstruction is not L1-regularized!');
elseif ~isnumeric(params.lambda) || ~isscalar(params.lambda) || ...
       ~isreal(params.lambda) || params.lambda<0
    error('Regularization parameter lambda must be a non-negativ real-valued scalar!');
elseif params.lambda==0
%     warning('Regularization parameter lambda equals to 0. Reconstruction is not L1-regularized!');
end  
%--------------------------------------------------------------------------
% nIter.
%--------------------------------------------------------------------------
if ~isfield(params,'nIter') || isempty(params.nIter)
    params.nIter = 30;
elseif ~isnumeric(params.nIter) || ~isscalar(params.nIter) || ...
       ~isreal(params.nIter) || params.nIter<0
    error('No. of iterations, nIter, must be a scalar larger than or equal to 0!');
end
%--------------------------------------------------------------------------
% nIterLS.
%--------------------------------------------------------------------------
if ~isfield(params,'nIterLS') || isempty(params.nIterLS)
    params.nIterLS = 150;
elseif ~isnumeric(params.nIterLS) || ~isscalar(params.nIterLS) || ...
       ~isreal(params.nIterLS) || params.nIterLS<=0
    error('No. of line-search iterations, nIterLS, must be a scalar larger than 0!');
end
%--------------------------------------------------------------------------
% gradTol.
%--------------------------------------------------------------------------
if ~isfield(params,'gradTol') || isempty(params.gradTol)
    params.gradTol = 1.0000e-030; 
elseif ~isnumeric(params.gradTol) || ~isscalar(params.gradTol) || ...
       ~isreal(params.gradTol) || params.gradTol<0
    error('');
end
%--------------------------------------------------------------------------
% l1Smooth.
%--------------------------------------------------------------------------
if ~isfield(params,'l1Smooth') || isempty(params.l1Smooth)
    params.l1Smooth = 1e-15; 
elseif ~isnumeric(params.l1Smooth) || ~isscalar(params.l1Smooth) || ...
       ~isreal(params.l1Smooth) || params.l1Smooth<0
    error('');
end
%--------------------------------------------------------------------------
% alpha.
%--------------------------------------------------------------------------
if ~isfield(params,'alpha') || isempty(params.alpha)
    params.alpha = 0.01; 
elseif ~isnumeric(params.alpha) || ~isscalar(params.alpha) || ...
       ~isreal(params.alpha) || params.alpha<0
    error('');
end
%--------------------------------------------------------------------------
% beta.
%--------------------------------------------------------------------------
if ~isfield(params,'beta') || isempty(params.beta)
    params.beta = 0.6; 
elseif ~isnumeric(params.beta) || ~isscalar(params.beta) || ...
       ~isreal(params.beta) || params.beta<0
    error('');
end
%--------------------------------------------------------------------------
% t0.
%--------------------------------------------------------------------------
if ~isfield(params,'t0') || isempty(params.t0)
    params.t0 = 1; 
elseif ~isnumeric(params.t0) || ~isscalar(params.t0) || ...
       ~isreal(params.t0) || params.t0<0
    error('');
end

end

function [ A_x, A_dx, Psi_x, Psi_dx ] = preObjective( x, dx, A, Psi, lambda )
%//////////////////////////////////////////////////////////////////////////
% Precalculates operators applied on x and dx in objective function.
%//////////////////////////////////////////////////////////////////////////

%==========================================================================
% L2-term.
%==========================================================================
if isa(A,'function_handle')
    A_x  = A(x,'notransp');
    A_dx = A(dx,'notransp');
else
    A_x  = A*x; 
    A_dx = A*dx;
end

%==========================================================================
% L1-term.
%==========================================================================
if lambda
    if isa(Psi,'function_handle')
        Psi_x  = Psi(x,'notransp');
        Psi_dx = Psi(dx,'notransp');
    else
        Psi_x  = Psi*x;
%         Psi_x = wavelet_transform(Psi,x);
        Psi_dx = Psi*dx;
%         Psi_dx = wavelet_transform(Psi,dx);
    end
else
    Psi_x = 0; Psi_dx = 0;
end

end

function [ f, l2Term, relRes ] = objective( A_x, A_dx, Psi_x, Psi_dx, y, lambda, l1Smooth, t )
%//////////////////////////////////////////////////////////////////////////
% Evaluates objective function at x + t*dx.
%//////////////////////////////////////////////////////////////////////////

f = 0; % Init.

%==========================================================================
% L2-term.
%==========================================================================
l2Term = norm( A_x(:) + t*A_dx(:) - y(:), 2 )^2;
f   = f + l2Term;

%==========================================================================
% L1-term.
%==========================================================================
if lambda
    tmp = Psi_x(:) + t*Psi_dx(:);
    tmp = sqrt( conj(tmp).*tmp + l1Smooth );
    f   = f + lambda*sum(tmp);
end

%==========================================================================
% Relative residual (of L2-term).
%==========================================================================
if nargout==3
    relRes = sqrt(l2Term) / norm(y(:),2);
end

end

function grad = gradObjective( x, A, Psi, y, lambda, l1Smooth )
%//////////////////////////////////////////////////////////////////////////
% Calculates gradient of objective function f at x.
%//////////////////////////////////////////////////////////////////////////

%==========================================================================
% L2-term.
%==========================================================================
if isa(A,'function_handle')
    grad = A( A(x,'notransp')-y, 'transp' );
else
    grad = A' * (A*x - y);
end

%==========================================================================
% L1-term.
%==========================================================================
if lambda
    if isa(Psi,'function_handle')
        x    = Psi(x,'notransp');
        x    = ( sqrt( conj(x).*x + l1Smooth ).^-1 ) .* x;
        grad = grad + lambda * Psi(x,'transp');
    else
        x    = Psi*x;
%         sz = size(x);
%         x = wavelet_transform(Psi,x);
        x    = ( sqrt( conj(x).*x + l1Smooth ).^-1 ) .* x;
        grad = grad + lambda * (Psi'*x);
%         grad = grad + lambda * inverse_wavelet_transform(Psi,x,sz);
    end
end

grad = 2*grad;

end


 




