% Define pott values
pott = [4, 5, 6, 7, 8, 9, 10];
%pott = [4,5,6];

% Initialize matrices to store execution times for each pott
% Each row corresponds to a different action (I, II, III)
% Each column corresponds to a different pott value
executionTimes = zeros(3, length(pott)); % 3 rows for I, II, III; columns for each pott

ggmmrreess = 1;
qqmmggrriitt = 1;
pprreeccoonndd = 1;

% Loop over pott values
for col = 1:length(pott)
    currentPott = pott(col);
    

if ggmmrreess
    try
tic;
%% Parameters
c = 1;
Lf = 1;
Ts = 0; Te = 1; T = Te - Ts;
Nx = 2^currentPott;%input('Enter N_x: '); % Example: 50
Nt = 2^currentPott;%input('Enter N_t: '); % Example: 50

dx = Lf / Nx;
dt = T / Nt;
x = linspace(0, Lf, Nx + 1).';
t = linspace(Ts, Te, Nt + 1);

%% Analytical Solution
u_exact = zeros(Nx + 1, Nt + 1);
for j = 1:Nx + 1
    for i = 1:Nt + 1
        u_exact(j, i) = sin(2 * pi * x(j)) * sin(4 * pi * t(i));
         %u_exact(j, i) = cos(2 * pi * x(j)) * cos(4 * pi * t(i));
    end
end

% Erstellen einer anonymen Funktion, die die zusätzlichen Parameter bindet
preconFunc = @(xx) dqW_s(xx, pott, mxt, cfc, ews);

% %% Plotting Exact Solution
% figure;
% surf(t, x, u_exact);
% xlabel('Time');
% ylabel('Spatial position');
% zlabel('Amplitude');
% title('Analytic solution');
% shading interp;
% colorbar;
t = t(1:end-1);
%% Matrix Construction
mainDiag = (1 + c^2 * dt^2 / dx^2) * ones(Nx + 1, 1);
mainDiag(1) = 1; mainDiag(end) = 1;
offDiag = -c^2 * dt^2 / (2 * dx^2) * ones(Nx , 1);
offDiag(1) = 0; offDiag(end) = 0;
A = diag(mainDiag) + diag(offDiag, 1) + diag(offDiag, -1);
A= sparse(A);
invA = sparse(inv(A));

D = (circshift(eye(Nx + 1), [1, 0]) - 2 * eye(Nx + 1) + circshift(eye(Nx + 1), [-1, 0])) / dx^2;
D(1,end) = 0; D(end,1) = 0;
%D(1,:) = 0; D(end,:) = 0;
D=sparse(D);

F_u = invA;
F_v = c^2 * dt * D * F_u;
F = [F_u, dt * F_u; F_v, eye(Nx + 1) + dt * F_v];

uV = ones(2 * (Nx + 1), 1);
uV(1)=0; uV(end)=0; uV(Nx + 1)=0; uV(Nx + 2)=0; 

F  =  diag(uV)*(F) *diag(uV) ;
%F  = (F) ;

F=sparse(F);

I_t = diag(ones(Nt - 1, 1), -1);
I_t(1,end) = 1;
I_t = sparse(I_t);
I_x = sparse(eye(2 * (Nx + 1)));
F_t = sparse(diag(ones(Nt, 1)));

% Constructing G matrix
G = sparse(kron(I_t, -F)) + sparse(kron(F_t, I_x));
G = sparse(G);

%% Right Hand Side Vector
vec = zeros(2 * (Nx + 1) * Nt, 1);
for n = 1:Nt
    rhs = [zeros(Nx + 1, 1); 1 *dt* f(x(1:end), t(n))];
    vec((2 * (Nx + 1) * (n - 1) + 1):2 * (Nx + 1) * n) = rhs;
end

%% Solving the System
%uv = G \ vec;
%uv = gmres(G,vec,[],10^-10, 300);

tol = 1e-5;      % Tolerance for convergence
maxit = Nx*Nt/2;     % Maximum number of iterations
%x0 = zeros(size(vec));
%x0 = reshape(uv(:,2:end), 1, [])';
x0 = rand(size(vec))-rand(size(vec));
%x0 = 1*ones(size(vec));

[uuvv, flag1, relres1, iter1, resvec1] = gmres(G, vec, [], tol, maxit, [], [], x0);
executionTimes(1, col) = toc;
catch ME
% Handle error
    warning('Error occurred: %s', ME.message);
    % Initialize outputs to some default or error values if necessary
    xxx = []; % Or some other default/error value
    flag = -1; % Custom flag to indicate error
    relres = []; % Or some other default/error value
    iter = -1; % Custom value to indicate error
    resvec = []; % Or some other default/error value
%[uuvv, flag] = gmres(G, vec, [], tol, maxit, [], [], x0);
executionTimes(1, col) = -1; % Store time for action I
    end
end

if qqmmggrriitt
%clear all
tic;
driver_qmgrit_Wellen(2^currentPott, 2^currentPott, 2, 1, 1, 3);
executionTimes(2, col) = toc;
end

if pprreeccoonndd
%clear all
try
tic;
BEwgf_s(currentPott, 1, 1, 1);
executionTimes(3, col) = toc;
catch ME
% Handle error
    warning('Error occurred: %s', ME.message);
    % Initialize outputs to some default or error values if necessary
    xxx = []; % Or some other default/error value
    flag = -1; % Custom flag to indicate error
    relres = []; % Or some other default/error value
    iter = -1; % Custom value to indicate error
    resvec = []; % Or some other default/error value
%[uuvv, flag] = gmres(G, vec, [], tol, maxit, [], [], x0);
executionTimes(3, col) = -1; % Store time for action I
end


end

end

% Create table headers reflecting Nt=Nx=2^pott
pottHeaders = "Nt=Nx=2^" + string(pott);

% Create a table
resultsTable_n = array2table(executionTimes, 'VariableNames', pottHeaders, 'RowNames', {'GMRES', 'QMGRIT', 'QMGRIT+GMRES'});

% Display the table
disp(resultsTable_n);

% % Display the execution times
% if ggmmrreess
% fprintf('\nExecution time for GMRES: %f seconds\n', gmrs);
% end
% if qqmmggrriitt
% fprintf('Execution time for QMGRIT: %f seconds\n', qmgrit);
% end
% if pprreeccoonndd
% fprintf('Execution time for QMGRIT + GMRES: %f seconds\n', precond);
% end


%% Force Function
function val = f(x, t)
    val = -12 * pi^2 * sin(4 * pi * t) .* sin(2 * pi * x);
    %val = -12 * pi^2 * cos(4 * pi * t) .* cos(2 * pi * x);
end