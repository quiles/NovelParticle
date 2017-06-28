clear all;

N = 10;
mat = rand(N,N);
mat = (mat + mat')/2;
for i=1:N
    mat(i,i) = 0;
end;

st = fopen('dados.txt', 'w');

fprintf(st,'%d\n',N);
for i=1:N-1
    for j=i+1:N
        fprintf(st,'%f ',mat(i,j));
    end;
end;
fclose(st);
        
mat