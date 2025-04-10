nn = size(QQ, 2);
QQQ_n = QQ;
if L==2
QQQ_a = Q2;
elseif L==3
QQQ_a = Q13_1;
elseif L==1
QQQ_a = EW;
end
figure; imagesc((QQQ_n)); title("Emerge: Numerical"); colorbar
figure; imagesc((QQQ_a)); title("Emerge: Analytical"); colorbar
figure; imagesc(abs(QQQ_n - QQQ_a)); title("Difference"); colormap('jet'); colorbar
eigenvaluesA = eig(full(QQQ_a));
eigenvaluesB = eig(QQQ_n);
s_eigenvaluesA = sort(eigenvaluesA);
s_eigenvaluesB = sort(eigenvaluesB);
max(s_eigenvaluesA)
max(s_eigenvaluesB)
if 1
maxDiff = abs(max(s_eigenvaluesA) - max(s_eigenvaluesB));
if maxDiff <= 1e-14
    disp('Die maximalen EV sind gleich bis zu einer Genauigkeit von 1e-14.');
else
    disp('Die maximalen EV sind nicht gleich bis zu dieser Genauigkeit.');
end
end
