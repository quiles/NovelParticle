clear all;

N = 1000;
inicio = 0;
step = 100;
time = 500;
s = 0; % box size
fname = 'data1000S/net_m0.60_r1';

si = ceil(time/step)

Vel = zeros(1,si);
Acc = zeros(1,si);
AccR = zeros(1,si);
AccA = zeros(1,si);

idX = [0];

nColors = 5;
cores = jet(nColors);
    
for i=inicio:step:time
    nomef = sprintf('%s.time_%d.par',fname,i);
    a = load(nomef);
    [N C] = size(a);
    [i N]
    nColors = max(a(:,2));
    cores = jet(nColors);

    for j=1:nColors
        index = find(a(:,2)==j);
        p = plot3(a(index,4),a(index,5),a(index,6),'.');
        hold on;
        set(p,'Color',cores(j,:), 'MarkerSize',30);
    end;
    box on;
    hold off;
    if (s>=1)
        axis([-s s -s s -s s]);     
    end;
    pause(0.1);
end;


xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
zlabel('x_3','FontSize',16);

