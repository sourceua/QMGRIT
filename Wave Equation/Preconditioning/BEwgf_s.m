function [vidpovid, flag, iter, MxDst] = BEwgf_s(pott, mxt, cfc, ews)
%% Parameters
c = 1;
Lf = 1;
Ts = 0; Te = 1; T = Te - Ts;
Nx = 2^pott;%input('Enter N_x: '); % Example: 50
Nt = 2^pott;%input('Enter N_t: '); % Example: 50

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
F=sparse(F);

uV = ones(2 * (Nx + 1), 1);
uV(1)=0; uV(end)=0; uV(Nx + 1)=0; uV(Nx + 2)=0; 

F  =  diag(uV)*(F) *diag(uV) ;
%F  = (F) ;

F=sparse(F);

I_t = (diag(ones(Nt - 1, 1), -1));
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
maxit = 5;%Nt*Nx;     % Maximum number of iterations
%x0 = zeros(size(vec));
%x0 = reshape(uv(:,2:end), 1, [])';
%x0 = rand(size(vec))-rand(size(vec));
%x0 = 1*ones(size(vec));
%tic;
[xxx, flag, relres, iter, ~] = gmres(G, vec, [], tol, maxit, preconFunc, []);
%toc;
%Analyzing the results
if flag == 0
    disp('GMRES converged to the desired tolerance');
else
    disp('GMRES did not converge');
end
disp(['Number of iterations: ', num2str(iter)]);
disp(['Final relative residual: ', num2str(relres)]);

U = xxx;
%% Extract Solution and Reshape
% U = zeros(Nx + 1, Nt + 1);
% for n = 1:Nt
%     index_start = (2 * (n - 1) * (Nx + 1) + 1);
%     index_end = 2 * (n - 1) * (Nx + 1) + (Nx + 1);
% 
%     u = [uuvv(index_start:index_end)];
%     U(:, n) = u;
% end

U = reshape(U, 2 * (Nx + 1), []);
U = U(1:(Nx + 1),:);
U = [U, U(:,1)];
%U(:,1:end-1) = U(:,2:end);
% %U(:,2:end) = U(:,1:end-1);
% %U(:, Nt+1)=U(:, 1);
% %U(:, 1)=U(:, Nt+1);
% %U = circshift(U, [0,2]);
% U(:,end) = U(:,1);
% %U = circshift(U, [0,1]);
% vidpovid = U;
%% Plotting Numerical Solution
if 0
t = linspace(Ts, Te, Nt + 1);
figure;
surf(t, x, U);
xlabel('Time');
ylabel('Spatial position');
zlabel('Amplitude');
title('Numerical solution - 2 variablen global');
shading interp;
colorbar;

figure;
surf(t,x,(u_exact(:,:)-U));
xlabel('Time'); ylabel('Spatial position'); zlabel('Amplitude');
title('difference'); shading interp; colorbar;

fprintf(1, "\n");
% fprintf(1, "%e\n", max((u_exact-U),[],"all"));
% fprintf(1, "\n");

MxDst = max((u_exact(:,1:end)-U),[],"all");
end
end

%% Force Function
function val = f(x, t)
    val = -12 * pi^2 * sin(4 * pi * t) .* sin(2 * pi * x);
    %val = -12 * pi^2 * cos(4 * pi * t) .* cos(2 * pi * x);
end


