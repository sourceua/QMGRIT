function [convergenceFactor, iterations] = driver_qmgrit_Wellen(Tn, Xn, m_cf, mu_cfcf, nu_ew, grids)
%close all;
%clear all;
%% plot options
options_plot   = {'LineWidth',3,'MarkerSize',14};
options_labels = {'FontSize',20};
%% problem parameters
plott=0;
c = 1;                   % diffusion coefficient
T  = 1;                    % final time
ch = 5;
Nt =Tn;        
t  = linspace(0,T,Nt+1);  % time domain
dt = T/Nt;                 % time-step size
X  = 1;                    % interval bound of spatial domain
Nx =Xn;            % number of dofs in space
dx = X/Nx;                 % spatial step size
x  = linspace(0,X,Nx+1)';  % spatial domain 
for potenz = 1:1
maxit = Nt*Nx/2;                % maximum number of iterations
test_len = 1;
convf = zeros(maxit,test_len);
convp = zeros(1,test_len);
for test=1:test_len
%% QMGRIT parameters
m     = m_cf; %2^potenz;                 % coarsening factor
L     = grids;                % number of grid levels
qum = nu_ew;%(2^(ch-1))/; %  number of eternal wandern iterations
cf = mu_cfcf;%test*(2^(ch-potenz)-1); %  number of cf smoothing iterations
tol   = 1e-5;% norm([dx,dt^2]);              % stopping tolerance
d = 0;
%% Force
app.b = @(x,t) ( -12 * pi^2 * sin(4*pi*t) .* sin(2*pi*x));
%% analytic solution
app.u_exact = zeros(Nx+1,Nt+1);
for j = 1:Nx+1
    for i = 1:Nt+1
        app.u_exact(j,i) = sin(4*pi*t(i)) .* sin(2*pi*x(j));
    end
end
%% algo initialisierung
% time values at each grid level
t = t(1:end-1);
tc = cell(L,1);
for l=1:L
    tc{l} = t(1:m^(l-1):end);
end
% setup matrix that acts in space for time integrator Phi
app.M = cell(L,1);
%e     = ones(Nx-1,1);
for l=1:L
 %% Matrix Construction
mainDiag = (1 + 1*c^2 * (((dt*m^(l-1))^2)) / dx^2) * ones(Nx + 1, 1);
mainDiag(1) = 1; mainDiag(end) = 1;
offDiag = - c^2 * (((dt*m^(l-1))^2)) / (2 * dx^2) * ones(Nx , 1);
offDiag(1) = 0; offDiag(end) = 0;
A = diag(mainDiag) + diag(offDiag, 1) + diag(offDiag, -1);
invA = inv(A);
%invA(invA<1e-3) = 0;
D = (circshift(eye(Nx + 1), [1, 0]) - 2 * eye(Nx + 1) + circshift(eye(Nx + 1), [-1, 0])) / dx^2;
%D(1,:) = 0; D(end,:) = 0;
D(1,end) = 0; D(end,1) = 0;

F_u = invA;
F_v = c^2 * (dt*m^(l-1)) * D * F_u;
F = [F_u, (dt*m^(l-1)) * F_u; F_v, eye(Nx + 1) + (dt*m^(l-1)) * F_v];

 uV = ones(2 * (Nx + 1), 1);
uV(1)=0; uV(end)=0; uV(Nx + 1)=0; uV(Nx + 2)=0; 

 app.M{l}  =diag(uV) * (F) ;
end
%% solution and right-hand side on each grid level
u = cell(L,1);
v = cell(L,1);
g = cell(L,1);
for l=1:L
   u{l} = 1*1e+1*rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l}));
   v{l} = 1*1e+1*rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l}));
   %u{l}(1)=0; u{l}(end)=0; 
 % v{l}(1)=0; v{l}(end)=0; 
   uv{l} = [u{l};v{l}];
  % g{l} = [zeros(size(u{l}, 1), 1); zeros(size(u{l}, 1), 1)];
   g{l} = zeros(size(uv{l}));
end
% modiefeied RHS
for i=1:Nt
       g{1}(:,i) = [zeros(size(u{1}, 1), 1); dt * app.b(x, t(i))];
end
%Phi = @(uv, tstom, tstart, app, l) PhiExecutor(uv, tstom, tstart, app, l, @(x) x);  
Phi = @(uv, tstom, tstart, app, l) PhiExecutor(uv, tstom, tstart, app, l);  
%% residual
res  = zeros(1,maxit+1);
r = zeros(size(uv{1}));
 % r(:,i) = g{1}(:,i)...
 %         + Phi(uv{1}(:,i-1),t(i),t(i-1),app,l)...
 %         - uv{1}(:,i);
for i=2:length(t)
  r(:,i) = g{1}(:,i)...
         + Phi(uv{1}(:,i-1),t(i),t(i-1),app,l)...
         - uv{1}(:,i);
end
% r(:,1) = g{1}(:,1)...
%        + Phi(uv{1}(:,length(t)-1),t(length(t)),t(length(t)-1),app,1)...
%        - uv{1}(:,length(t));

res(1) = norm(r);
tt=res(1);
fprintf(1,'iter  1: norm(r) = %e\n', res(1))
iteration = 1;
%% QMGRIT iterations
for nu = 2:maxit+1

     [uv,g] = qmgrit_wellen_direct(1,L,uv,g,Phi,app,m,tc, qum, cf, iteration,d);   

 %r = [zeros(size(u{l}, 1), 1); zeros(size(u{l}, 1), 1)];
 %r = [zeros(size(u{1})); zeros(size(u{1}))];
 r = zeros(size(uv{1}));
 
 % r(:,i) = g{1}(:,i)...
 %         + Phi(uv{1}(:,i-1),t(i),t(i-1),app,l)...
 %         - uv{1}(:,i);
    for i=2:length(t)
         r(:,i) = g{1}(:,i)...
                + Phi(uv{1}(:,i-1),t(i),t(i-1),app,1)...
                - uv{1}(:,i);
    end
    % r(:,1) = g{1}(:,1) ...
    %    + Phi(uv{1}(:,length(t)-1),t(length(t)),t(length(t)-1),app,1)...
    %    - uv{1}(:,length(t));

    res(nu) = norm(r);
    fprintf(1,'iter %2d: norm(r) = %e\n', nu, res(nu));
    fprintf(1, "%e\n", res(nu)/res(nu-1));
    convf(nu-1,test) = res(nu)/res(nu-1);

    if res(nu) < tol
        miter = nu;
        break;
    end
 iteration=iteration+1;   
end
mu=norm(uv{1}(:,1)-uv{1}(:,end));
fprintf(1,'Periodisch ? %e\n', mu);
%% Plots
if plott
% Convergency
figure;
semilogy(0:nu-1,res(1:nu),'*-',options_plot{:});
title([int2str(L),'-level Lin. QMGRIT, m=',int2str(m)])
set(gca,options_labels{:})
% plot solution
figure;
surf(t,x,uv{1}(1:size(u{l}, 1),:));
%title('computed solution - 2 variablen iterativ')
xlabel('t'); ylabel('x');
%set(gca,options_labels{:})
zlabel('Amplitude');
title('Numerical solution - 2 variablen iterativ');
shading interp;
colorbar;

fprintf(1, "\n");
fprintf(1, "%e\n", max((app.u_exact-uv{1}(1:size(u{l}, 1),:)),[],"all"));
fprintf(1, "\n");
end
convp(1,test)=norm(eigs(app.M{1}^((cf+1)*m),1,'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', 2*(Nx+1)));
%eigs((app.M{1}^m-app.M{2})*app.M{1}^(cf*m),1,'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', 66)
end
convnum = zeros(1, test_len);
% for i=1:test_len
%   convf(convf(:,i)==0) = [];
% end
convf(convf==0) = [];
for i=1:test_len
    %convnum(1,i) = median(convf(end-50:end,i));
    convnum(1,i) = geomean(convf(:,i)); 
end
convergenceFactor=convnum;
iterations = length(convf);

if plott
figure;
stem(linspace(1,test_len,test_len), convnum, 'r');
hold on
stem(linspace(1,test_len,test_len), convp, 'b');
legend('nume','pred')
title([int2str(L),'-level Lin. QMGRIT, m=',int2str(m)])
grid on; zoom on;
hold off
end
end
%% Plotting Difference
% figure;
% surf(t,x,(app.u_exact-uv{1}(1:size(u{l}, 1),:)));
% xlabel('Time'); ylabel('Spatial position'); zlabel('Amplitude');
% title('difference'); shading interp; colorbar
%% Phi - Time Integrator
function uv_new = PhiExecutor(uv, tstom, tstart, app, l)%, postProcessFunc)
    u = uv(1:length(uv)/2);
    v = uv(length(uv)/2+1:end);
  %  u(1) = 0; u(end) = 0;
   % v(1) = 0; v(end) = 0;

    uv_old = [u; v];
    uv_new = app.M{l} * uv_old;

    u_new = uv_new(1:length(u));
    v_new = uv_new(length(u)+1:end);
 %   u_new(1) = 0; u_new(end) = 0;
  %  v_new(1) = 0; v_new(end) = 0;
    uv_new = [u_new; v_new];
end
end
