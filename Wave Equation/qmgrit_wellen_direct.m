function [uv,g] = qmgrit_wellen_direct(l, L, uv, g, Phi, app, m, t, qum, cf, iteration, d)
           %[u,g] = qmgrit(l, L, u, g, Phi, app, m, t, iter, qum,x,tt)
plot = 0;
uebrtrg = 0;
direct = d;
ress=0;
omega= 1;%0.61803398875;
%tmp = zeros(1,qum); 
if l == L
   %% coarse-grid solve
   for nu = 1:qum
       if direct
      uv{l} = direct_solve_2v(l, g, app);
       elseif direct==0
    for i=2:length(t{l})
        uv{l}(:,i) = g{l}(:,i)...
                  + feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);
    end    
     uv{l}(:,1)=g{l}(:,end)...
                  + feval(Phi,uv{l}(:,end-1),t{l}(end),t{l}(end-1),app,l);
  if uebrtrg
     uv{l}(:,1)=uv{l}(:,end);
  end
   
if ress
     r = zeros(size(uv{l}));
    for i=2:length(t{l})
         r(:,i) = g{l}(:,i)...
                + Phi(uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l)...
                - uv{l}(:,i);
    end
    r(:,1) = g{l}(:,1) ...
       + Phi(uv{l}(:,length(t{l})-1),t{l}(length(t{l})),t{l}(length(t{l})-1),app,l)...
       - uv{l}(:,length(t{l}));
    fprintf(1,'EWiter %2d: EWnorm(r) = %e\n', nu, norm(r));
end    
       end
 
  % if mod(nu,qum/1)==0 && 0 && mod(iteration, 10)==0
  if plot
     figure;
surf(t{l},linspace(0,1,2^5+1),uv{l}(1:length(uv{l})/2,:));
title('computed coarse displacement')
xlabel('t'); ylabel('x')
 figure;
surf(t{l},linspace(0,1,2^5+1),uv{l}(length(uv{l})/2+1:end,:));
title('computed coarse velocity')
xlabel('t'); ylabel('x')
   end
   end  
else
    
    % MGRIT algorithm
    %% F- and C-points on current grid level
    nt           = length(t{l});
    Cpts         = 1:m:nt;
    Fpts         = 1:nt;
    Fpts(Cpts)   = [];
    % fprintf(1,'********* start MGRIT *********\n');
    % fprintf(1,'  l                  = %d\n', l);
    % fprintf(1,'  nt                 = %d\n', nt);
    % fprintf(1,'  m                  = %d\n', m);
    % fprintf(1,'  number of F-points = %d\n', length(Fpts));
    % fprintf(1,'  number of C-points = %d\n', length(Cpts));

    %% F-relaxation
%   if ((l > 1) || ( (iter == 1) && (l == 1) ))
 %      for i=Fpts
   %    uv{l}(:,i) = g{l}(:,i)...
 %                + feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);
%       end
%   end
for p=1:cf
%% C-relaxation  
    for i=Cpts(2:end)
    %   uv{l}(:,i) = g{l}(:,i)...
      %           + feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);

% uv{l}(:,i) = (1 - omega) * uv{l}(:,i) + omega * (g{l}(:,i)...
  %           + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l));
  new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
    end
  %uv{l}(:,1) = g{l}(:,1) + ...
   %feval(Phi,uv{l}(:,length(t{l})-1),t{l}(length(t{l})),t{l}(length(t{l})-1),app,l);
  new_value = g{l}(:,1) + feval(Phi, uv{l}(:,length(t{l})-1), t{l}(1), t{l}(length(t{l})-1), app, l);
    uv{l}(:,1) = uv{l}(:,1) + omega * (new_value - uv{l}(:,1));
   
%uv{l}(:,1) =  (1 - omega) * uv{l}(:,1) + omega * (g{l}(:,1) + ...
  % feval(Phi,uv{l}(:,length(t{l})-1),t{l}(length(t{l})),t{l}(length(t{l})-1),app,l));
  if uebrtrg
   uv{l}(:,1) = uv{l}(:,end);
  end  
    %% F-relaxation
    for i=Fpts
         %   uv{l}(:,i) = g{l}(:,i)...
           %     + feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);
%uv{l}(:,i) = (1 - omega) * uv{l}(:,i) + omega * (g{l}(:,i)...
  %          + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l));
  new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
    end

%%
 if 0
% plot solution
figure;
surf(t{l},linspace(0,1,2^5+1),uv{1}(1:length(uv{1})/2,:));
title('computed solution presmoothing')
xlabel('t'); ylabel('x');
 end
end
    %% Compute Residual at Coarse-Points
for i=2:length(Cpts)
   g{l+1}(:,i) = g{l}(:,Cpts(i)) + ...
   feval(Phi,uv{l}(:,Cpts(i)-1),t{l}(Cpts(i)),t{l}(Cpts(i)-1),app,l)...
               - uv{l}(:,Cpts(i));         
end
g{l+1}(:,Cpts(1)) = g{l}(:,Cpts(1)) +...
                  feval(Phi,uv{l}(:,Cpts(end)-1),t{l}(Cpts(end)),t{l}(Cpts(end)-1),app,l)-...
                   uv{l}(:,Cpts(1));
if uebrtrg
g{l+1}(:,Cpts(1)) = g{l+1}(:,(end));
end
    %% next level
    uv{l+1}  = zeros(size(uv{l+1}));
    [uv,g] = qmgrit_wellen_direct(l+1, L, uv, g, Phi, app, m, t, qum, cf, iteration, d);
  %   uv{l+1}(:,1)=uv{l+1}(:,end);                                  
    %% correct the approximation u at C-points
    uv{l}(:,Cpts) = uv{l}(:,Cpts) + uv{l+1};
 if plot
    % plot solution
figure;
surf(t{l},linspace(0,1,2^5+1),uv{1}(1:length(uv{1})/2,:));
title('computed solution after CGC')
xlabel('t'); ylabel('x');
 end
   %u{l}(:,end)=u{l}(:,1);
%   uv{l}(:,1)=uv{l}(:,end);
    %% carry out F-relax to correct the approximation u at F-points
  for i=Fpts
%         if i==2
%             uv{l}(:,i)  = g{l}(:,i)...
%                 + omega*feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);
%         else 
%             uv{l}(:,i) = g{l}(:,i)...
%                 + omega*feval(Phi,uv{l}(:,i-1),t{l}(i),t{l}(i-1),app,l);
%         end
new_value = g{l}(:,i) + feval(Phi, uv{l}(:,i-1), t{l}(i), t{l}(i-1), app, l);
    uv{l}(:,i) = uv{l}(:,i) + omega * (new_value - uv{l}(:,i));
  end
   if plot
 figure;
surf(t{l},linspace(0,1,2^5+1),uv{1}(1:length(uv{1})/2,:));
title('computed solution after postsmothing')
xlabel('t'); ylabel('x');
   end
end