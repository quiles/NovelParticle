clear all;

inicio = 1;
step = 1;
time = 66;
s = 0; % box size
fname = '../../Drosofila/networks/drosophila_subset_';

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

teste = zeros(step,3);

for i=inicio:step:time
    nomef = sprintf('%st%d.par',fname,i);
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

    FA = a(:,C-1);
    FR = a(:,C);
    
    DF = (FA - FR).^2;

    FA = sort(a(:,C-1));
    FR = sort(a(:,C));
    
    figure(Forces);
%     plot(DF);
    teste(i,:) = [norm(FA) norm(FR) norm(DF)];
    plot(FA,'b');
    hold on;
    plot(FR,'r');
    box on;
%     axis([1 N 0 1.1]);
    hold off;
    
    pause(0.1);
end;

figure(Part);
xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
zlabel('x_3','FontSize',16);

figure;
time = 1:66;
plot(time(1:30),teste(1:30,1),'r');
hold on;
plot(time(31:40),teste(31:40,1),'c');
plot(time(41:58),teste(41:58,1),'b');
plot(time(59:66),teste(59:66,1),'g');
ylabel('FA','FontSize',16);

figure;
plot(time(1:30),teste(1:30,2),'r');
hold on;
plot(time(31:40),teste(31:40,2),'c');
plot(time(41:58),teste(41:58,2),'b');
plot(time(59:66),teste(59:66,2),'g');
ylabel('FR','FontSize',16);

figure;
plot(time(1:30),teste(1:30,3),'r');
hold on;
plot(time(31:40),teste(31:40,3),'c');
plot(time(41:58),teste(41:58,3),'b');
plot(time(59:66),teste(59:66,3),'g');
ylabel('DF','FontSize',16);

figure;
dados = load('tempo.txt');
dados(:,3) = dados(:,4) ./ dados(:,3);
plot(time(1:30),dados(1:30,3),'r');
hold on;
plot(time(31:40),dados(31:40,3),'c');
plot(time(41:58),dados(41:58,3),'b');
plot(time(59:66),dados(59:66,3),'g');
ylabel('Centroids','FontSize',16);

figure;
plot(time(1:30),dados(1:30,2),'r');
hold on;
plot(time(31:40),dados(31:40,2),'c');
plot(time(41:58),dados(41:58,2),'b');
plot(time(59:66),dados(59:66,2),'g');
ylabel('Transient','FontSize',16);
