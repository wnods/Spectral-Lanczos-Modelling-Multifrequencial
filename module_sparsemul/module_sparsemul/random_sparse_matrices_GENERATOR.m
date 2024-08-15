% create sparse random matrices
clc;

na = 4000; ma = 3000;
nb = 3000; mb = 5000;
a = full(sprand(na,ma,.3)); dimA = [na,ma,nnz(a)];
b = full(sprand(nb,mb,.4)); dimB = [nb,mb,nnz(b)];
%b = eye(3); dimB = [3,3,nnz(b)];

% na = 3; ma = 3;
% nb = 3; mb = 3;
% a = [0 10 0; 0 0 20; 30 0 0]; dimA = [na,ma,nnz(a)];
% b = [40 0 70; 50 0 0; 0 60 0];dimB = [nb,mb,nnz(b)];

save dimA.txt dimA  -ascii
save matA.txt a -double -ascii
save dimB.txt dimB  -ascii
save matB.txt b -double -ascii

disp('Program Finished!')