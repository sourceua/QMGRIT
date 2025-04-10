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
L     = 3;                % number of grid levels
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
%e     = ones(Nx-1,1);
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

 app.M{l}  =  diag(uV) * (F) * diag(uV);
end
I = sparse(diag(ones(2*(Nx+1),1)));
n = size(app.M{1}, 1); 
Nc = Nt/m; % Number of blocks (coarse)
Ncc = Nc/m; % Number of blocks (coarse-coarse)
% Creating the identity matrix of size n*n
I_n = sparse(eye(n*Nt));
I_nN = sparse(eye(n*Nc));
I_nNc = sparse(eye(n*Ncc));
% Creating the block circulant matrix
Ap = sparse(zeros(n*Nt)); % Initialize the block matrix with zeros
A = sparse(zeros(n*Nt)); % Initialize the block matrix with zeros
Apc = sparse(zeros(n*Nc)); % Initialize the matrix with zeros
Ac = sparse(zeros(n*Nc)); % Initialize the matrix with zeros
Acc = sparse(zeros(n*Ncc)); % Initialize the coarse coarse block  matrix with zeros
As = sparse(zeros(n*Nc)); % Initialize the coarse block circulant matrix with zero
Ass = sparse(zeros(n*Ncc)); % Initialize the coarse coarse block circulant matrix with zero
% Create the larger matrix 
P1 = sparse(zeros(m*n*(Nt/m), n*(Nt/m)));
P2 = sparse(zeros(m*n*(Ncc), n*(Ncc)));
R1 = sparse(zeros( n*(Nt/m) , m*n*(Nt/m)));
R2 = sparse(zeros( n*(Ncc) , m*n*(Ncc)));
% 
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
Ap=sparse(Ap);
% 
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

%% Apc
% Filling the diagonal with I blocks
for i = 1:Nc
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    Apc(idx, idx) = I;
end
% Placing 'a' on the subdiagonal
for i = 2:Nc
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    Apc(idx_row, idx_col) = -app.M{2};
end
% Placing 'a' in the top right block
Apc(1:n, end-n+1:end) = -app.M{2}; % Top right corner
Apc=sparse(Apc);

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

%% 3rd level
%% Acc
% Filling the diagonal with I blocks for Acc
for i = 1:Ncc
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    Acc(idx, idx) = I;
end
% Placing 'a' on the subdiagonal for Acc
for i = 2:Ncc
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    Acc(idx_row, idx_col) = -app.M{3};
end
Acc = sparse(Acc);

%% Ass
% Filling the diagonal with I blocks for Ass
for i = 1:Ncc
    idx = (i-1)*n + (1:n); % Indexing the correct block position
    Ass(idx, idx) = I;
end
% Placing 'a' on the subdiagonal for Ass
for i = 2:Ncc
    idx_row = (i-1)*n + (1:n);
    idx_col = (i-2)*n + (1:n);
    Ass(idx_row, idx_col) = -app.M{2}^m;
end
% Placing 'a' in the top right block for Ass
Ass(1:n, end-n+1:end) = -app.M{2}^m; % Top right corner
Ass = sparse(Ass);

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

%% 3rtd Level
%% P2
% Fill the matrix block by block for P2
for colBlock = 0:(Ncc-1)
    for rowBlock = 0:(m-1)
        rowStart = (colBlock*m + rowBlock)*n + 1;
        colStart = colBlock*n + 1;
        
        if rowBlock == 0
            % Identity matrix for the top of each column block
            P2(rowStart:rowStart+n-1, colStart:colStart+n-1) = eye(n);
        else
            % app.M{2} raised to the power of rowBlock for subsequent blocks
            P2(rowStart:rowStart+n-1, colStart:colStart+n-1) = app.M{2}^rowBlock;
        end
    end
end
P2 = sparse(P2);

%% R2
% Fill R2 row by row, with each row block [I, a, a^2] in sequence for each Ncc
for colBlock = 1:Ncc
    colStart = (colBlock-1)*m*n + 1;
    R2((colBlock-1)*n+1:colBlock*n, colStart:colStart+n-1) = eye(n); % Place I
end
R2 = sparse(R2);

 EW = (I_n - A\Ap);
 EW = sparse(EW);

%fprintf(1,'\n')
fprintf(1,'Three Level, disk=%i\n', Nt)

%%
% Calculate the length of the vector
vecLength = 2 * (Nx + 1) * Nt;
vecLengthG = 2 * (Nx + 1) * Nt/2;
vecLengthGG = 2 * (Nx + 1) * Nt/4;
% Initialize the vector with ones
unityVector = ones(vecLength, 1);
unityVectorG = ones(vecLengthG, 1);
unityVectorGG = ones(vecLengthGG, 1);
% Set the specified indices to zero
for k = 0:(2 * Nt - 1)
    baseIdx = k * (1 * (Nx + 1));
    indices = [baseIdx + 1, baseIdx + (Nx + 1), baseIdx + (Nx + 1) + 1];
    unityVector(indices) = 0;
end
for k = 0:(2 * Nt/2 - 1)
    baseIdx = k * (1 * (Nx + 1));
    indices = [baseIdx + 1, baseIdx + (Nx + 1), baseIdx + (Nx + 1) + 1];
    unityVectorG(indices) = 0;
end
for k = 0:(2 * Nt/4 - 1)
    baseIdx = k * (1 * (Nx + 1));
    indices = [baseIdx + 1, baseIdx + (Nx + 1), baseIdx + (Nx + 1) + 1];
    unityVectorGG(indices) = 0;
end
unityVector = unityVector(1:end-1);
unityVectorG = unityVectorG(1:end-1);
unityVectorGG = unityVectorGG(1:end-1);
ONE = diag(unityVector);
ONEG = diag(unityVectorG);
ONEGG = diag(unityVectorGG);
%%

S = (I_nN - As);

S2 = (I_nNc - Ass);

CGC2 =  sparse(I_nNc - Acc \ Ass) ; % coarse grid correction for the second and third levels

Q23_1 =  P2 * CGC2 * S2 * R2; % Q23 matrix

CGC3_1 =   (I_nN - ((eye(size(Q23_1)) - Q23_1) * (Apc^(-1) * As))) ; 

Q13_1 = P1 * CGC3_1 * S^1 * R1 ; % Q13 matrix


til = 1;
[Vec, Val] = eigs(Q13_1, til, 'largestabs', 'Tolerance', 1e-19, 'SubspaceDimension', Nt*2*(Nx+1));
eigenvalues = diag(Val);

% Print each eigenvalue with its norm
for k = 1:length(eigenvalues)
    fprintf('Q1 Eigenvalue %d: %f + %fi, Norm: %f\n', k, real(eigenvalues(k)), imag(eigenvalues(k)), abs(eigenvalues(k)));
end
end
end
