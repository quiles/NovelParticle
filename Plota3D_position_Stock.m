clear all;

inicio = 0;
step = 1;
time = 100;
s = 0; % box size
fname = '/Users/quiles/Networks/Stock-Exchange/19871222.time_';

nColors = 5;
cores = jet(nColors);
 
Part = figure;

teste = zeros(step,3);

for i=inicio:step:time
    nomef = sprintf('%s%d.par',fname,i);
    a = load(nomef);
    [N C] = size(a);
    nColors = max(a(:,3));
    [i N nColors]

    cores = jet(nColors);
    figure(Part);
    hold off;
    if (nColors > 0) 
        for j=1:nColors
            index = find(a(:,3)==j);
            p = plot3(a(index,4),a(index,5),a(index,6),'.');
            hold on;
            set(p,'Color',cores(j,:), 'MarkerSize',30);
        end;
    else
        p = plot3(a(:,4),a(:,5),a(:,6),'.');
        set(p,'Color',[0 0 0], 'MarkerSize',30);
    end;

    box on;
    hold off;
    if (s>=1)
        axis([-s s -s s -s s]);     
    end;

    
    pause(0.1);
end;

