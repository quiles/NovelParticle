clear all;

data = load('UCI_Msg.dat');

A = zeros(1899,1899);
A(1,1) = 100;

for t=1:59
    for i=1000*(t-1)+1:1000*t
        A(data(i,1),data(i,2)) = 75-t;
        A(data(i,2),data(i,1)) = 75-t;
    end;
    imagesc(A);
    t
    pause(0.5);
end;
