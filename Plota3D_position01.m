clear all;

inicio = 0;
step = 1;
time = 50;
s = 0; % box size
fname = 'network';

si = ceil(time/step)

Vel = zeros(1,si);
Acc = zeros(1,si);
AccR = zeros(1,si);
AccA = zeros(1,si);

idX = [0];

nColors = 5;
cores = jet(nColors);
 
Part = figure;
Forces = figure;

for i=inicio:step:time
    nomef = sprintf('%s.time_%d.par',fname,i);
    a = load(nomef);
    [N C] = size(a);
    [i N]
    nColors = max(a(:,2));
    cores = jet(nColors);
    figure(Part);
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

    FA = a(:,C-1);
    FR = a(:,C);
    
    DF = (FA - FR).^2;

    FA = sort(a(:,C-1));
    FR = sort(a(:,C));
    
    figure(Forces);
%     plot(DF);
    [norm(FA) norm(FR) norm(DF)]
    plot(FA,'b');
    hold on;
    plot(FR,'r');
    box on;
%     axis([1 N 0 1.1]);
    hold off;
    
    pause(1.5);
end;

figure(Part);
xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
zlabel('x_3','FontSize',16);

