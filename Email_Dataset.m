clear all;

% data = load('email-Eu-core-Dept3.txt');
data = load('email-Eu-core.txt');

[L C] = size(data);

m1 = max(max(data(:,1:2)));

A = zeros(m1,m1);

for e=1:L
    i = data(e,1) + 1;
    j = data(e,2) + 1;
    t = data(e,3);
    A(i,j) = 1;
    A(j,i) = 1;
end;

imagesc(A);