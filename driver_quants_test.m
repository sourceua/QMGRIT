% QMGRIT Tests for the Periodic Complex Valued Equation towards 
% Schroedinger Equation in 1D
% This script provides a framework for testing various RHS functions of the Schroedinger equation
% using the QMGRIT method, with documentation and visualization features.
%
% Usage:
%   - Run this script in MATLAB to visualize solutions for predefined RHS functions.
%   - Adjust 'coarsing', 'grids', 'Nx', and 'Nt' parameters to explore their impact on the solution.
%
% Inputs:
%   ch - Base resolution for spatial and temporal discretization.
%   Nx - Number of spatial discretization points.
%   Nt - Number of temporal discretization points.
%   coarsing - Coarsening factor for the QMGRIT method, affecting time grid refinement.
%   grids - Number of grid levels used in the multigrid hierarchy of QMGRIT.
%   maxiter - Maximum number of QMGRIT iterations allowed.
%   gamma - Maximum number of eternal wandern iterations.
%   nu - Correction factor, adjust based on convergence needs.
%   TOL - Stopping criterion based on residual norm tolerance.
%   diff - Diffusion coefficient in the partial differential equation (PDE).
%   Tend - Final time for the simulation.
%   dt - Time step size.
%   Xend - Final spatial point for the simulation.
%   dx - Spatial step size.
%
% Outputs:
%   Graphical visualizations of analytical and numerical solutions corresponding to the
%   tested RHS functions, and computed solution with convergence behavior;
%   parallel efficency testing
%
% Example Adjustments:
%   - To increase resolution, set 'ch' to a higher value (e.g., 2^9).
%   - For a more aggressive coarsening, increase 'coarsing' (e.g., 4).
%
% Performance Tips:
%   - 'Nx' and 'Nt' should be chosen based on the desired accuracy and available computational resources.

%close all
%clear all

ch = 2^17; app.Nx = 5; app.Nt = ch;       % Number of spatial and time intervals
coarsing = 2^1;                             % Coarsening level for time grid refinement
grids = 3;                                  % Number of grids in the multigrid hierarchy
plot_ana = 0;                               % 0: no plots
plot_num = 0;                               % 1: plot numerical solutions; 0: no plots
plot_err = 0;                               % 0: no plots
plot_res = 0;                               % 0: no plots
paralleltest = 1;                           % 1: performing a test with 2 workers; 0: no tests
maxiter = 5;                               % Maximum number of QMGRIT iterations allowed
gamma = 1;                                  % Maximum number of eternal wandern iterations
nu = 1;                                     % Maximum number of CF iterations
TOL = 9e-9;                                % Stopping criterion based on residual norm tolerance
app.diff = 1;                               % Coefficient in the PDE
app.Tend = 1;                            % Final time for the simulation
app.dt = app.Tend / app.Nt;                 % Time step size
app.T = linspace(0, app.Tend,app.Nt+1);     % Time domain discretization
app.Xend = 1;                               % Final spatial point for the simulation
app.dx = app.Xend / app.Nx;                 % Spatial step size
app.X = linspace(0, app.Xend, app.Nx+1)';   % Spatial domain discretization
 
%% Spatial operator setup for the time integrator
% Preallocate cell array for matrices at different levels
app.M = cell(grids, 1);
e = ones(app.Nx-1, 1);
% Assemble the matrix for each level considering the coarsening in time
for l = 1:grids
    r = (app.dt * coarsing^(l-1)) / (app.dx^2);
    app.M{l} = (1/1i)*spdiags( ...
        [app.diff * r * e, (1i - 2 * app.diff * r) * e, app.diff * r * e], ...
        -1:1, app.Nx-1, app.Nx-1);
end
%% Time integrator function definition
% Defines the time integration step using a simple backward Euler scheme
Phi = @(ustart, app, l) (app.M{l} \ ustart);

% RHS functions and corresponding analytical solutions
u_exact_functions = {
    @(x, t) 'idk'
};

rhs_functions = {
  @(x, t) 1 * (x) .* (x-1) * (1 * cos(3 * 2 * t) + 1 * 1i * sin(1 * 2 * t))
};

% Iterating over functions to solve and plot
for fun = 1:length(u_exact_functions)
    u_exact = u_exact_functions{fun};
    rhs = rhs_functions{fun};
    
    %% Time grid setup for multilevel integration
    % Adjust the time vector to exclude the last time point because of periodicity 
    Tp = app.T(1:end-1);
    % Preallocate cell array for time grids at different levels
    app.tc = cell(grids, 1);
    % Populate the cell array with coarsened time grids
    for l = 1:grids
        app.tc{l} = Tp(1:coarsing^(l-1):end);
    end
    %% Initial solution and RHS setup for each grid level
    app.g = cell(grids, 1); % Preallocate cell array for RHS vectors at different levels
    % Zero initialization for RHS
    for l = 1:grids
        app.g{l} = 0 * (rand(app.Nx-1, numel(app.tc{l})) ...
                        - rand(app.Nx-1, numel(app.tc{l})));
    end
    % Compute the initial RHS for the finest level based on the source term
    for j = 1:app.Nt
        app.g{1}(:, j) = app.M{1} \ (app.dt *(1i^-1)* rhs(app.X(2:end-1), app.T(j)));
    end

    times=[];
    delete(gcp('nocreate')); % Delete any existing parallel pool
    % Call output_driver_qmgrit for the current RHS function
    tic
    u = Qmgrit(coarsing, grids, 0, maxiter, gamma, nu, TOL, Phi, app, plot_res);
    times = [times, toc];
    
    if paralleltest
    numWorkersArray = 2:2; % Array to hold the number of workers to test
    for i = 1:length(numWorkersArray)
         numWorkers = numWorkersArray(i);
         delete(gcp('nocreate')); % Delete any existing parallel pool
         parpool(numWorkers); % Create a new parallel pool with the desired number of workers

    
         tic; % Start timing
         u = Qmgrit(coarsing, grids, 1, maxiter, gamma, nu, TOL, Phi, app, plot_res);
         times = [times, toc]; % Store the execution time
    end
    
    figure;
    plot(1:numWorkers, log10(times), '-o');
    xlabel('Number of Workers');
    ylabel('Logarithmic Execution Time');
    title('quants Algorithm Performance with Varying Number of Workers');
    grid on;
    xticks(1:numWorkers); % Set x-axis ticks to integers
    end

    % if plot_ana 
    %     % Plotting the analytical solution
    %     figure;
    %     surf(app.X(2:end-1), app.T, u_exact(app.X(2:end-1), app.T)');
    %     xlabel('x');
    %     ylabel('t');
    %     zlabel('u(x, t)');
    %     title(['Analytical Solution for Function ', num2str(fun)]);
    % end
    if plot_num
        % Computed solution plot
        if 0
        figure;
        surf(app.X(2:end-1),app.T(1:end-1),real(u'));
        title('computed solution real')
        xlabel('x'); ylabel('t');
        set(gca,options_labels{:})
        figure;
        surf(app.X(2:end-1),app.T(1:end-1),imag(u'));
        title('computed solution imag')
        xlabel('x'); ylabel('t');
        set(gca,options_labels{:})
        figure;
        surf(app.X(2:end-1),app.T(1:end-1),abs(u'));
        title('computed solution abs')
        xlabel('x'); ylabel('t');
        set(gca,options_labels{:})
        end
        %%plot solution on complex plane
        figure;
        % Plot the initial values (t = 0)
        plot(real(u(:,1)), imag(u(:,1)), 'ro', 'MarkerSize', 4, 'LineWidth', 1.5);
        hold on;
        % Plot the solution at each time step
        for i = 2:app.Nt
            plot(real(u(:,i)), imag(u(:,i)), 'b.-', 'MarkerSize', 4, 'LineWidth', 0.5);
        end
        % Add labels and title
        xlabel('Real');
        ylabel('Imaginary');
        title('Solution on Complex Plane');
        grid on;
        axis equal;
        % Add legend
        legend('Initial Values', 'Solution');
    end
    % if plot_err 
    %     % Computed solution plot
    %     figure;
    %     surf(app.X(2:end-1), app.T(1:end-1), ...
    %         u'-u_exact(app.X(2:end-1), app.T(1:end-1))');
    %     title('quants Error');
    %     xlabel('x'); ylabel('t');
    % end 

end
