%Corr adjancency matrix

clear all;

listA = load('network2.dat');

[m lixo] = size(listA);

n = max(max(listA));

A = zeros(n,n);

for i=1:m
    A(listA(i,1),listA(i,2)) = 1;
    A(listA(i,2),listA(i,1)) = 1;
end;

c = corr(A);

[L D] = eig(c);
