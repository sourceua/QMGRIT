%% problem parameters
tic
plott_ana=1;
plott_sol =1;
c = 1;                   % diffusion coefficient
T  = 1;                    % final time
ch = 6;
Nt = 2^ch;        
t  = linspace(0,T,Nt+1);  % time domain
dt = T/Nt;                 % time-step size
X  = 1;                    % interval bound of spatial domain
Nx = 2^ch;            % number of dofs in space
dx = X/Nx;                 % spatial step size
x  = linspace(0,X,Nx+1)';  % spatial domain 
for potenz = 1:1
maxit = 100;                % maximum number of iterations
test_len =1;
convf = zeros(maxit,test_len);
convp = zeros(1,test_len);
for test=1:test_len
%% QMGRIT parameters
m     = 16^potenz;                 % coarsening factor
L   = 2;                % number of grid levels
%qum = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30];%(2^(ch-1))/; %  number of eternal wandern iterations
qum = 1;
EWeff = [];
EWiter = [];
%for qum = 1:1
%convf = zeros(maxit,test_len);
cf = test*(2^(ch-potenz)-1); %  number of cf smoothing iterations
d = 0; %directe solve on coarse grid
tol   =  1e-7; %norm([dx,dt])/1;              % stopping tolerance
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
%if d ==1 || l==2 || l==3 
%D(1,:) = 0; D(end,:) = 0;
%end
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
   u{l} = 1*1e+1*(rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l})));
   v{l} = 1*1e+1*(rand(Nx+1,numel(tc{l}))-rand(Nx+1,numel(tc{l})));
   %u{l}(1)=0; u{l}(end)=0; 
 % v{l}(1)=0; v{l}(end)=0; 
  uv{l} = [u{l};v{l}];
   g{l} = zeros(size(uv{l}));
end
if 0
uv{1}=randi([0,1], 2 * (Nx + 1) , numel(tc{1}) );
end
if 0
A = zeros(2 * (Nx + 1), numel(tc{1}));  % Nullmatrix
iii = randi(size(A,1));  % zufällige Zeile
jjj = randi(size(A,2));  % zufällige Spalte
%iii=5;
%jjj=1;
A(iii,1) = 1;  % genau eine 1 setzen
A(jjj,1) = 1;  % genau eine 1 setzen
uv{1}=A;
end
% modiefeied RHS
for i=1:Nt+1
       g{1}(:,i) = [zeros(size(u{1}, 1), 1); dt * app.b(x, t(i))];
end
%Phi = @(uv, tstom, tstart, app, l) PhiExecutor(uv, tstom, tstart, app, l, @(x) x);  
Phi = @(uv, tstom, tstart, app, l) PhiExecutor(uv, tstom, tstart, app, l);  
%% residual
res  = zeros(1,maxit+1);
r = zeros(size(uv{1}));
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

     [uv,g] = qmgrit_wellen_direct(1,L,uv,g,Phi,app,m,tc, qum, cf, iteration, d);   

 r = zeros(size(uv{1}));

    for i=2:length(t)
         r(:,i) = g{1}(:,i)...
                + Phi(uv{1}(:,i-1),t(i),t(i-1),app,1)...
                - uv{1}(:,i);
    end
    r(:,1) = g{1}(:,1) ...
       + Phi(uv{1}(:,length(t)-1),t(length(t)),t(length(t)-1),app,1)...
       - uv{1}(:,length(t));

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
u{1} = [u{1}(:,end-1),u{1}(:,1:end-1)];
mu=norm(uv{1}(:,1)-uv{1}(:,end));
fprintf(1,'Periodisch ? %e\n', mu);
%% Plots
if plott_sol
% Convergency
figure;
semilogy(0:nu-1,res(1:nu),'*-');
title([int2str(L),'-level Lin. QMGRIT, m=',int2str(m)])
% plot solution
figure;
surf(t,x,uv{1}(1:size(u{l}, 1),:));
title('computed solution - 2 variablen iterativ')
xlabel('t'); ylabel('x');
%set(gca,options_labels{:})
zlabel('Amplitude');
%title('Numerical solution - 2 variablen iterativ');
shading interp;
colorbar;
fprintf(1, "\n");
fprintf(1, "%e\n", max((app.u_exact-uv{1}(1:size(u{l}, 1),:)),[],"all"));
fprintf(1, "\n");
end
if L==1
 nuka = norm(eigs(F,1,'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', 2*(Nx+1)));
convp(1,test) = nuka^Nt;
else
  convp(1,test)=norm(eigs(app.M{1}^((cf+1)*m),1,'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', 2*(Nx+1)));
end
end

convnum = zeros(1, test_len);
for i = 1:test_len
    % Select non-zero elements of the i-th column, excluding the first 4 rows if needed
    nonZeroElements = convf(5:end, i);
    nonZeroElements = nonZeroElements(nonZeroElements ~= 0);
    % Calculate the geometric mean of the non-zero elements
    if isempty(nonZeroElements)
        convnum(1, i) = NaN; % Use NaN or another placeholder if there are no non-zero elements
    else
        convnum(1, i) = geomean(nonZeroElements);
    end
end
convergenceFactor=convnum;
iterations = length(convf);

EWeff = [EWeff; convergenceFactor];
EWiter = [EWiter, iterations];
%end
if plott_ana
figure;
% Create stem plot for convp with blue markers
stem(1:test_len, convp, 'Color', 'b', 'Marker', '>', 'MarkerSize', 14);
hold on;
% Create stem plot for convnum with orange markers
stem(1:test_len, convnum, 'Color', [0.8500, 0.3250, 0.0980], 'Marker', '<', 'MarkerSize', 14);

% Adjust x-axis to have integer ticks
xticks(1:test_len);
%xlim([1 test_len]);

% Concatenate L, m, and mu into the title, converting numbers to strings
title([int2str(L), ' - level QMGRIT, m = ', int2str(m), ', N_t = N_x = 2^', int2str(ch)], 'FontSize', 16);

% Increase font size for axes numbering
ax = gca; % Current axes
ax.FontSize = 16; % Change axis font size

xlabel('q, whereas \nu = qN_c-1');
ylabel('convergence factor');

% Increase font size for legend
legend('predicted', 'numerical',  'FontSize', 16, 'Location', 'southeast');

% Turn on grid and zoom
grid on; 
zoom on;

% Release the hold on the current figure
hold off;
if 0
 % Save figures to files
     figName = sprintf(['Study/L%d/ConvergenceAnalysis_CF' ...
         '_m%d_ch%d.fig'], L, m, ch);
     figName1 = sprintf(['Study/L%d/ConvergenceAnalysis_CF' ...
         '_m%d_ch%d.jpeg'], L, m, ch);
     saveas(gcf, figName);
     saveas(gcf, figName1);
end
end

%% Plotting Difference
if 0
figure;
surf(t,x,(app.u_exact-uv{1}(1:size(u{l}, 1),:)));
xlabel('Time'); ylabel('Spatial position'); zlabel('Amplitude');
title('difference'); shading interp; colorbar
end
toc
end
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