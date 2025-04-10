function [uv,g] = qmgrit_wellen_direct(l, L, uv, g, Phi, app, m, t, qum, cf, iteration, d)
direct = d;
omega=1;
if l == L
   %% coarse-grid solve
   for nu = 1:qum
       if direct
      uv{l} = direct_solve_2v(l, g, app);
       elseif direct==0
    uv{l}(:,1)=g{l}(:,1)...
                  + feval(Phi,uv{l}(:,end),t{l}(end),t{l}(end-1),app,l);
    for i=2:length(t{l})
        uv{l}(:,i) = g{l}(:,i)...
                  + feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);              
    end
       end
   end  
else
    
    % QMGRIT algorithm
    %% F- and C-points on current grid level
    nt           = length(t{l});
    Cpts         = 1:m:nt;
    Fpts         = 1:nt;
    Fpts(Cpts)   = [];
    %% F-relaxation
  if ((l > 1) || ( (iteration == 1) && (l == 1) ))
      for i=Fpts
        new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
      end
  end
for p=1:cf
%% C-relaxation 
 new_value = g{l}(:,(1)) + feval(Phi, uv{l}(:,(end)), t{l}(1), t{l}(length(t{l})-1), app, l);
    uv{l}(:,Cpts(1)) = uv{l}(:,(1)) + omega * (new_value - uv{l}(:,(1)));
    for i=Cpts(2:end)
  new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
    end
    %% F-relaxation
    for i=Fpts
  new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
    end
end
    %% Compute Residual at Coarse-Points

    g{l+1}(:,1) = g{l}(:,(1)) +...
                  feval(Phi,uv{l}(:,end),t{l}(Cpts(end)),t{l}(Cpts(end)-1),app,l)-...
                   uv{l}(:,(1));
    for i=2:length(Cpts)
   g{l+1}(:,i) = g{l}(:,Cpts(i)) + ...
   feval(Phi,uv{l}(:,Cpts(i)-1),t{l}(Cpts(i)),t{l}(Cpts(i)-1),app,l)...
               - uv{l}(:,Cpts(i));         
end

    %% next level
    uv{l+1}  = zeros(size(uv{l+1}));
    [uv,g] = qmgrit_wellen_direct(l+1, L, uv, g, Phi, app, m, t, qum, cf, iteration, d);
    %% correct the approximation u at C-points
    uv{l}(:,Cpts) = uv{l}(:,Cpts) + uv{l+1};
    %% carry out F-relax to correct the approximation u at F-points
  for i=Fpts
new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
  end
end