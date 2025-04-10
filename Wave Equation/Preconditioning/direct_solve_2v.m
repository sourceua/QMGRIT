function u_direkt_new = direct_solve_2v(l, g, app)
Nt = size(g{end}, 2);
dt = 1/Nt;  
% if l > 1 
% g{end} = (dt^1*g{end});
% end
f = g{end}(:,1:end);
%f = g{end}(:,2:end);
Nx = size(f, 1)/2 - 1;

F = sparse(app.M{end});

I_t = sparse(diag(ones(Nt-1, 1), -1));
I_t(1,end) = 1;
I_x = sparse(eye(2 * (Nx + 1)));
F_t = sparse(diag(ones(Nt, 1)));

% Constructing G matrix
G = sparse(kron(I_t, -F)) + sparse(kron(F_t, I_x));
G = sparse(G);

%% Right Hand Side Vector
%vec = zeros(2 * (Nx + 1) * Nt, 1);
vec = reshape(f, [], 1);
%for n = 1:Nt
   % rhs = [zeros(Nx + 1, 1); dt * f(:,n)];
  %  vec(( (Nx + 1) * (n - 1) + 1): (Nx + 1) * n) = f(:,n);% rhs;
%end

%% Solving the System
%uv = G \ vec;
uv = gmres(G,vec,[],1e-4, 300);

%% Extract Solution and Reshape
    U = zeros(Nx + 1, Nt);
    V = zeros(Nx + 1, Nt);
    for n = 1:Nt
        index_start_u = (2 * (n - 1) * (Nx + 1) + 1);
        index_end_u = (2 * (n - 1) * (Nx + 1)) + Nx + 1;
        index_start_v = index_end_u + 1;
        index_end_v = index_start_v + Nx;

        U(:, n ) = uv(index_start_u:index_end_u);
        V(:, n ) = uv(index_start_v:index_end_v);
    end

  %  U(:,1:end-1) = U(:,2:end);
    %U(:, Nt+1) = U(:, 1);
    %V(:,1:end-1) = V(:,2:end);
    %V(:, Nt+1) = V(:, 1);
% %
   u_direkt_new = [U;V];
end


