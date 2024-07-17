%% QMGRIT Solution for Periodic Heat Equation in 1D
function u = Qmgrit(coarsing, grids, paral, maxiter, gamma, nu, TOL, Phi, app, plot_res)
% Solves periodic equations using QMGRIT method.
%
% Detailed algorithm explanation, performance optimization tips, and robust error handling are included
% to assist users in adapting the script to their specific needs.
%
% Inputs:
%   rhs - Function handle for the RHS of the heat equation.
%   coarsing - Time grid coarsening factor.
%   grids - Number of multigrid levels.
%   paral - Flag indicating whether to use parallel computation.
%   maxiter - Maximum number of QMGRIT iterations allowed.
%   gamma - Maximum number of eternal wandern iterations.
%   nu - Correction factor, adjust based on convergence needs.
%   TOL - Stopping criterion based on residual norm tolerance.
%   Xend - Final spatial point for the simulation.
%   Phi - Function handle for the time integration step.
%   app - Structure containing problem-specific parameters, has to have Nt,
%         T, Nx, tc, g.
%
% Validate input parameters
if app.Nx <= 0 || app.Nt <= 0 || coarsing <= 0
    error('Nx, Nt, and coarsing must be positive integers');
end

%% Initialization of QMGRIT
par = paral;
m = coarsing;                              % Standard coarsening factor for time grid refinement
L = grids;                                 % Number of grid levels in the multigrid hierarchy
maxit = maxiter;                           % Maximum number of QMGRIT iterations allowed
qum = gamma;                               % Maximum number of eternal wandern iterations
cf = nu;                                   % Correction factor, adjust based on convergence needs
tol = TOL;                                 % Stopping criterion based on residual norm tolerance
conv = [];                                 % Convergence 
nt = app.Nt;                               % Number of time steps
t = app.T; %linspace(0, T, nt+1)';         % Time domain discretization
% Adjust the time vector to exclude the last time point because of periodicity 
t = t(1:end-1);
% Initial solution and RHS setup for each grid level
u = cell(L, 1); % Preallocate cell array for solution vectors at different levels
for l = 1:L
    % Set the seed for the random number generator
    seed = 42; % Use any integer value
    rng(seed);
    % Generate random numbers
    u{l} = 1 * (rand(app.Nx-1, numel(app.tc{l})) - rand(app.Nx-1, numel(app.tc{l})));
end
%% Residual computation for initial guess
res = zeros(1, maxit+1); % Preallocate array for storing residuals
r = zeros(size(u{1})); % Preallocate array for the initial residual computation
    % Compute residual for the first time step
    r(:, 1) = app.g{1}(:, 1) + Phi(u{1}(:, end), app, 1) - u{1}(:, 1);
    % Compute residuals for subsequent time steps
    for j = 2:length(t)
        r(:, j) = app.g{1}(:, j) + Phi(u{1}(:, j-1), app, 1) - u{1}(:, j);
    end
% Store the norm of the initial residual
g = app.g;
res(1) = norm(r);
fprintf(1, 'iter  1: norm(r) = %e\n', res(1));
%% Main QMGRIT iteration loop
for nus = 2:maxit+1
    % Perform a QMGRIT iteration step
    [u, g] = qmgrit(1, L, u, g, Phi, app, m, app.tc, nus, qum, cf);
    r = zeros(size(u{1})); % Reset the residual vector for the current iteration

    r(:, 1) = g{1}(:, 1) + Phi(u{1}(:, end), app, 1) - u{1}(:, 1);
    % Compute residuals for subsequent time steps
    for j = 2:length(t)
        r(:, j) = g{1}(:, j) + Phi(u{1}(:, j-1), app, 1) - u{1}(:, j);
    end

    % Update and print the residual norm for the current iteration
    res(nus) = norm(r);
    fprintf(1, 'iter %2d: norm(r) = %e\n', nus, res(nus));
    fprintf(1, "%e\n", res(nus) / res(nus-1));
    conv = [conv, res(nus) / res(nus-1)];
    % Break the loop if the solution converges based on the tolerance
    if res(nus) < tol
        miter = nus;
        break;
    end
end
u = u{1};
%% Post-processing and visualization
% Convergence behavior plot
    if plot_res
    % Convergence behavior plot
    figure;
    semilogy(0:miter, res(1:miter+1), '*-');
    title([int2str(L), '-level Lin. QMGRIT, m=', int2str(m)]);    
    end
% Display the geometric mean of the convergence factor
fprintf(1, 'Conv fac is %e\n', geomean(conv(2:end-1)));
%% Recursive QMGRIT Application
function [u,g] = qmgrit(l, L, u, g, Phi, app, m, t, iter, qum, cf)          
if l == L
   %% coarse-grid solve
   for ew = 1:qum
    u{l}(:,1)=g{l}(:,1)...
                  + feval(Phi,u{l}(:,end),app,l);
    for k=2:length(t{l})
        u{l}(:,k) = g{l}(:,k)...
                  + feval(Phi,u{l}(:,k-1),app,l);       
    end
   end   
else
    %% QMGRIT algorithm
    %% F- and C-points on current grid level
    nt           = length(t{l});
    Cpts         = 1:m:nt;
    Fpts         = 1:nt;
    Fpts(Cpts)   = [];
    %% F-relaxation
   if ((l > 1) || ( (iter == 2) && (l == 1) ))
      if par == 0
       for k=Fpts
            u{l}(:,k) = g{l}(:,k)...
                + feval(Phi,u{l}(:,k-1),app,l);
        end   
      end
    if par == 1
        tmpa = circshift(u{l}, [0,0]);
        tmpa = tmpa(:,Fpts-1);
        tmpg = circshift(g{l}, [0,0]);
        tmpg = tmpg(:,Fpts);
        Fpts_range = 1:numel(Fpts);
        tmp = zeros(size(u{l},1),numel(Fpts));
           
        parfor k = Fpts_range
               tmp(:,k) = tmpg(:,k) + feval(Phi, tmpa(:,k), app, l);
        end

        u{l}(:,Fpts) = tmp;
    end 
   end
%% CF-relaxations
for p=1:cf
   %% C-relaxation
%% C-relaxation with parallel option
if par == 0 % Sequential
    % Handle the first Cpt
    u{l}(:,Cpts(1)) = g{l}(:,Cpts(1)) + feval(Phi, u{l}(:,Fpts(end)), app, l);
    % Sequentially update the rest of the Cpts
    for k = Cpts(2:end)
        u{l}(:,k) = g{l}(:,k) + feval(Phi, u{l}(:,k-1), app, l);
    end
elseif par == 1 % Parallel
    % Sequentially handle the first Cpt due to its unique dependency
    u{l}(:,Cpts(1)) = g{l}(:,Cpts(1)) + feval(Phi, u{l}(:,Fpts(end)), app, l);
    % Prepare for parallel execution
    tmpa = circshift(u{l}, [0,1]); % Shift u{l} right to align previous u{l}(:,k-1) for k in Cpts
    tmpg = g{l}(:,Cpts(2:end)); % Corresponding g{l} values for the parallel section
    Cpts_range = 2:numel(Cpts); % Indices for the parallel section, excluding the first handled Cpt
    tmp = zeros(size(u{l},1), numel(Cpts_range)); % Storage for parallel results
    parfor k = 1:numel(Cpts_range)
        idx = Cpts_range(k); % Actual index in Cpts
        % Adjust indexing within tmpa and tmpg to match parfor's 1-based indexing
        tmp(:,k) = tmpg(:,k) + feval(Phi, tmpa(:,Cpts(idx)), app, l);
    end
    % Assign the computed values back to u{l}, excluding the first Cpt
    u{l}(:,Cpts(2:end)) = tmp;
end
    %% F-relaxation
    if par == 0
       for k=Fpts
            u{l}(:,k) = g{l}(:,k)...
                + feval(Phi,u{l}(:,k-1),app,l);
        end   
    end
    if par == 1
        tmpa = circshift(u{l}, [0,0]);
        tmpa = tmpa(:,Fpts-1);
        tmpg = circshift(g{l}, [0,0]);
        tmpg = tmpg(:,Fpts);
        Fpts_range = 1:numel(Fpts);
        tmp = zeros(size(u{l},1),numel(Fpts));

        parfor k = Fpts_range
               tmp(:,k) = tmpg(:,k) + feval(Phi, tmpa(:,k), app, l);
        end

        u{l}(:,Fpts) = tmp;
    end  
end
   %% Compute Residual at Coarse-Points
% Handle the first Cpt
g{l+1}(:,1) = g{l}(:,Cpts(1)) + feval(Phi,u{l}(:,Fpts(end)),app,l) - u{l}(:,Cpts(1));
% Sequentially update the rest of the Cpts
for k = 2:length(Cpts)
    g{l+1}(:,k) = g{l}(:,Cpts(k)) + feval(Phi,u{l}(:,Cpts(k)-1),app,l) - u{l}(:,Cpts(k));
end
    %% next level
    u{l+1}  = zeros(size(u{l+1}));
    [u,g] = qmgrit(l+1, L, u, g, Phi, app, m, t, iter, qum, cf);    
    %% correct the approximation u at C-points
    u{l}(:,Cpts) = u{l}(:,Cpts) + u{l+1};
    %% carry out F-relax to correct the approximation u at F-points
     if par == 0
       for k=Fpts
            u{l}(:,k) = g{l}(:,k)...
                + feval(Phi,u{l}(:,k-1),app,l);
       end   
     end
     if par == 1
        tmpa = circshift(u{l}, [0,0]);
        tmpa = tmpa(:,Fpts-1);
        tmpg = circshift(g{l}, [0,0]);
        tmpg = tmpg(:,Fpts);
        Fpts_range = 1:numel(Fpts);
        tmp = zeros(size(u{l},1),numel(Fpts));

           parfor k = Fpts_range
                tmp(:,k) = tmpg(:,k) + feval(Phi, tmpa(:,k), app, l);
           end

        u{l}(:,Fpts) = tmp;
     end 
end
end
end
