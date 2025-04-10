QQ = [];

%% problem parameters
tic
plott = 0;
plottt = 0;
gifff = 0;
c = 1;                   % diffusion coefficient
T  = 1;                    % final time
ch = 4;
Nt = 2^ch;        
dt = T/Nt;                 % time-step size
X  = 1;                    % interval bound of spatial domain
Nx = 2^ch;            % numberof dofs in space
dx = X/Nx;                 % spatial step size
x  = linspace(0,X,Nx+1)';  % spatial domain 
tatat = [];

dimension = (Nt) * 2 * (Nx + 1) * 1;

for basvec =  1 : dimension

    basisVector = zeros(dimension, 1);
    basisVector(basvec) = 1;

for potenz = 1:1
maxit = 1;                % maximum number of iterations
test_len = 1;
convf = zeros(maxit,test_len);
convp = zeros(1,test_len);
for test=1:test_len
%% QMGRIT parameters
m     = 2^1;%potenz;                 % coarsening factor
L   = 2;                % number of grid levels
EWeff = [];
EWiter = [];
for qum =1:1
convf = zeros(maxit,test_len);
cf = 1;%(2^(ch-potenz)-1);%2;%test*(2^(ch-potenz)-1); %  number of cf smoothing iterations
d = 0; %directe solve on coarse grid
tol   =  3.258000e-12; %norm([dx,dt])/1;              % stopping tolerance
%% Force
app.b = @(x,t) ( -12 * pi^2 * sin(4*pi*t) .* sin(2*pi*x));
%% analytic solution
t  = linspace(0,T,Nt+1);
app.u_exact = zeros(Nx+1,Nt+1);
for j = 1:Nx+1
    for i = 1:Nt+1
        app.u_exact(j,i) = sin(4*pi*t(i)) .* sin(2*pi*x(j));
    end
end
%% algo initialisierung
t = t(1:end-1);
% time values at each grid level
tc = cell(L,1);
for l=1:L
    tc{l} = t(1:m^(l-1):end);
end
% setup matrix that acts in space for time integrator Phi
app.M = cell(L,1);
for l=1:L
 %% Matrix Construction
mainDiag = (1 + 1*c^2 * (((dt*m^(l-1))^2)) / dx^2) * ones(Nx + 1, 1);
mainDiag(1) = 1; mainDiag(end) = 1;
offDiag = - c^2 * (((dt*m^(l-1))^2)) / (2 * dx^2) * ones(Nx , 1);
offDiag(1) = 0; offDiag(end) = 0;
A = diag(mainDiag) + diag(offDiag, 1) + diag(offDiag, -1);
invA = inv(A);
D = (circshift(eye(Nx + 1), [1, 0]) - 2 * eye(Nx + 1) + circshift(eye(Nx + 1), [-1, 0])) / dx^2;
D(1,end) = 0; D(end,1) = 0;

F_u = invA;
F_v = c^2 * (dt*m^(l-1)) * D * F_u;
F = [F_u, (dt*m^(l-1)) * F_u; F_v, eye(Nx + 1) + (dt*m^(l-1)) * F_v];

uV = ones(2 * (Nx + 1), 1);
uV(1)=0; uV(end)=0; uV(Nx + 1)=0; uV(Nx + 2)=0; 

 app.M{l}  = diag(uV) * (F) * diag(uV) ;
end
%% solution and right-hand side on each grid level
u = cell(L,1);
v = cell(L,1);
g = cell(L,1);
for l=1:L
   u{l} = 0*1e1*(rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l})));
   v{l} = 0*1e1*(rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l})));
  uv{l} = [u{l};v{l}];
   g{l} = zeros(size(uv{l}));
end
uv{1}=reshape(basisVector, 2 * (Nx + 1), []);
% modiefeied RHS
for i=1:Nt
     g{1}(:,i) = [zeros(size(u{1}, 1), 1); 0*dt * app.b(x, t(i))];
end
Phi = @(uv, tstom, tstart, app, l) PhiExecutor(uv, tstom, tstart, app, l);  
%% residual
res  = zeros(1,maxit+1);
r = zeros(size(uv{1}));
 r(:,1) = g{1}(:,1) ...
       + Phi(uv{1}(:,length(t)),t(length(t)),t(length(t)-1),app,1)...
       - uv{1}(:,1);
for i = 2:length(t)
  r(:,i) = g{1}(:,i)...
         + Phi(uv{1}(:,i-1),t(i),t(i-1),app,l)...
         - uv{1}(:,i);
end
res(1) = norm(r);
tt = res(1);
fprintf(1,'iter  1: norm(r) = %e\n', res(1))
iteration = 1;
%% QMGRIT iterations
for nu = 2:maxit+1

     [uv,g] = qmgrit_wellen_direct(1,L,uv,g,Phi,app,m,tc, qum, cf, iteration, d);   

 r = zeros(size(uv{1}));

 r(:,1) = g{1}(:,1) ...
       + Phi(uv{1}(:,length(t)),t(length(t)),t(length(t)-1),app,1)...
       - uv{1}(:,1);
    for i = 2:length(t)
         r(:,i) = g{1}(:,i)...
                + Phi(uv{1}(:,i-1),t(i),t(i-1),app,1)...
                - uv{1}(:,i);
    end
   

    res(nu) = norm(r);
    fprintf(1,'iter %2d: norm(r) = %e\n', nu, res(nu));
    fprintf(1, "%e\n", res(nu)/res(nu-1));
    convf(nu-1,test) = res(nu)/res(nu-1);

    if res(nu) < tol
        miter = nu;
        break;
    end
 iteration = iteration+1;   
end

QQ = [QQ, reshape(uv{1}(:,1:end),[],1)];
end
end
end
end

%% Phi - Time Integrator
function uv_new = PhiExecutor(uv, tstom, tstart, app, l)%, postProcessFunc)
    u = uv(1:length(uv)/2);
    v = uv(length(uv)/2+1:end);
    u(1) = 0; u(end) = 0;
    v(1) = 0; v(end) = 0;

    uv_old = [u; v];
    uv_new = app.M{l} * uv_old;

    u_new = uv_new(1:length(u));
    v_new = uv_new(length(u)+1:end);
    u_new(1) = 0; u_new(end) = 0;
   v_new(1) = 0; v_new(end) = 0;
    uv_new = [u_new; v_new];
end