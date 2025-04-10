c = 1;                   % diffusion coefficient
T  = 1;                    % final time
ch = 4;
Nt =2^ch;        
t  = linspace(0,T,Nt+1);  % time domain
dt = T/Nt;                 % time-step size
X  = 1;                    % interval bound of spatial domain
Nx =2^ch;            % number of dofs in space
dx = X/Nx;                 % spatial step size
x  = linspace(0,X,Nx+1)';  % spatial domain 
for potenz = 1:1
maxit = 1;                % maximum number of iterations
test_len = 1;
for test=1:test_len
%% QMGRIT parameters
m     = 2^potenz;                 % coarsening factor
L     = 2;                % number of grid levels
qum = 1;%(2^(ch-1))/; %  number of eternal wandern iterations
cf = 1;%1*test*(2^(ch-potenz)-1); %  number of cf smoothing iterations
tol   = 1e-12;              % stopping tolerance
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

 app.M{l}  = diag(uV) * (F) *  diag(uV);
end
I = diag(ones(2*(Nx+1),1));
% Given matrix 'a' and identity matrix 'I'
n = size(app.M{1}, 1); % Assuming 'a' is n by n
Nc = Nt/m; % Number of blocks
% Creating the identity matrix of size n*n
I_n = eye(n*Nt);
I_nN = eye(n*Nc);
% Creating the block circulant matrix
A = zeros(n*Nt);
Ap = zeros(n*Nt); % Initialize the block circulant matrix with zeros
Ac = zeros(n*Nc); % Initialize the coarse block circulant matrix with zeros
As = zeros(n*Nc); % Initialize the coarse block circulant matrix with zero
% Create the larger matrix 
P1 = zeros(m*n*(Nt/m), n*(Nt/m));
R1 = zeros( n*(Nt/m) , m*n*(Nt/m));

%% A
% Filling the diagonal with I blocks
for i = 1:Nt
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    A(idx, idx) = I;
end
% Placing 'a' on the subdiagonal
for i = 2:Nt
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    A(idx_row, idx_col) = -app.M{1};
end
%% Ap
% Filling the diagonal with I blocks
for i = 1:Nt
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    Ap(idx, idx) = I;
end
% Placing 'a' on the subdiagonal
for i = 2:Nt
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    Ap(idx_row, idx_col) = -app.M{1};
end
% Placing 'a' in the top right block
Ap(1:n, end-n+1:end) = -app.M{1}; % Top right corner

%% Ac
% Filling the diagonal with I blocks
for i = 1:Nc
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    Ac(idx, idx) = I;
end
% Placing 'a' on the subdiagonal
for i = 2:Nc
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    Ac(idx_row, idx_col) = -app.M{2};
end
Ac=sparse(Ac);

%% As
% Filling the diagonal with I blocks
for i = 1:Nc
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    As(idx, idx) = I;
end
% Placing 'a' on the subdiagonal
for i = 2:Nc
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    As(idx_row, idx_col) = -app.M{1}^m;
end
% Placing 'a' in the top right block
As(1:n, end-n+1:end) = -app.M{1}^m; % Top right corner
As=sparse(As);

%% P1
% Fill the matrix block by block
for colBlock = 0:(Nt/m)-1 % Three column blocks
    for rowBlock = 0:(m-1) % Each column block has 'm' row blocks
        rowStart = (colBlock*m + rowBlock)*n + 1;
        colStart = colBlock*n + 1;
        
        if rowBlock == 0
            % Identity matrix for the top of each column block
            P1(rowStart:rowStart+n-1, colStart:colStart+n-1) = eye(n);
        else
            % app.M{1} raised to the power of rowBlock for subsequent blocks
            P1(rowStart:rowStart+n-1, colStart:colStart+n-1) = app.M{1}^rowBlock;
        end
    end
end
P1=sparse(P1);

%% R1
% Fill MM row by row, with each row block [I, a, a^2] in sequence for each Nc
for colBlock = 1:Nc  % Iterate over the number of column blocks
    % Calculate the starting column index for the current block
    colStart = (colBlock-1)*m*n + 1; 
    % Place the identity matrix and the powers of 'a' in the correct positions
    R1((colBlock-1)*n+1:colBlock*n, colStart:colStart+n-1) = eye(n); % Place I
end
R1=sparse(R1);

EW = (I_n - A\Ap);
EW = sparse(EW);
CGC = (I_nN - Ac\As);
CGC=sparse(CGC);
S = (I_nN - As);
S=sparse(S);
Q2=P1 * CGC * S * R1;
matrices = {EW, Q2, P1 * CGC * S^2 * R1, P1 * CGC * S^3 * R1, P1 * CGC * S^4 * R1, P1 * CGC * S^5 * R1};
matrix_names = {'EW', 'P1 * CGC * S * R1', 'P1 * CGC * S^2 * R1', 'P1 * CGC * S^3 * R1','P1 * CGC * S^4 * R1','P1 * CGC * S^5 * R1'};
% Preallocate an array to store the fifth eigenvalues
first_eigenvalues = zeros(1, length(matrices));
second_eigenvalues = zeros(1, length(matrices));
third_eigenvalues = zeros(1, length(matrices));
fourth_eigenvalues = zeros(1, length(matrices));
fifth_eigenvalues = zeros(1, length(matrices));
% Iterate over each matrix and compute its fifth eigenvalue
for i = 1:length(matrices)
    eigenvalues = eigs(matrices{i}, 5, 'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', Nt*2*(Nx+1));
    first_eigenvalues(i) = norm(eigenvalues(1));
     second_eigenvalues(i) = norm(eigenvalues(2));
     third_eigenvalues(i) = norm(eigenvalues(3));
    fourth_eigenvalues(i) = norm(eigenvalues(4));
    fifth_eigenvalues(i) = norm(eigenvalues(5));
    % Display the matrix name and its fifth eigenvalue
    fprintf('1st. %s: %f\n', matrix_names{i}, first_eigenvalues(i));
    fprintf('2st. %s: %f\n', matrix_names{i}, second_eigenvalues(i));
  fprintf('3st. %s: %f\n', matrix_names{i}, third_eigenvalues(i));
  fprintf('4st. %s: %f\n', matrix_names{i}, fourth_eigenvalues(i));
     fprintf('5th. %s: %f\n', matrix_names{i}, fifth_eigenvalues(i));
     fprintf(1,'\n')
end
end
end