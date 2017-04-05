clear all;

N = 128;
inicio = 0;
step = 10;
time = 300;
s = 0; % box size
fname = './net_m0.10_r1';

si = ceil(time/step)

Vel = zeros(1,si);
Acc = zeros(1,si);
AccR = zeros(1,si);
AccA = zeros(1,si);

idX = [0];

nColors = 5;

cores = jet(nColors);
% cores = [166 206 227; % 1
%     31 120 180; % 2
%     178 223 138; % 3
%     51 160 144; % 4
%     251 154 153; % 5
%     227 26 28; % 6
%     254 191 111; % 7
%     255 127 0; % 8
%     202 178 214]% 9
% cores = cores / 255;
    
% correc = [1; 2; 3; 4; 1];

count = 10;
for i=inicio:step:time
%     subplot(1,2,1);
    nomef = sprintf('%s.time_%d.par',fname,i);
    a = load(nomef);
    [N C] = size(a);
    [i N]

    FA = zeros(1,N);
    FR = zeros(1,N);
    FRes = zeros(1,N);
    FTotal = 0;
    subplot(1,2,1);
    idxcom = find(a(:,4) == mod(count,16)+1);
    for j=1:N
%         if (max(ismember(idX,j))==1)
%             p = plot3(a(j,1),a(j,2),a(j,3),'*');
%         else
            p = plot3(a(j,1),a(j,2),a(j,3),'.');
%         end;
            
        hold on;
%         set(p, 'MarkerSize',50);
%        if (a(j,4) == mod(count,16)+1) 
            set(p,'Color',cores(a(j,4)+1,:), 'MarkerSize',40);
%        else set(p,'Color',[0 0 0], 'MarkerSize',20);
 %       end;
        FA(j) = (a(j,6)^2 + a(j,7)^2 + a(j,8)^2)^(1/2);
        FR(j) = (a(j,9)^2 + a(j,10)^2 + a(j,11)^2)^(1/2);
        FRes(j) = ((a(j,6)-a(j,9))^2 + (a(j,7)-a(j,10))^2 + (a(j,8)-a(j,11))^2)^(1/2);
        FTotal = FTotal + FRes(j);
%         accA = accA + (a(j,12)^2 + a(j,13)^2 + a(j,14)^2)^(1/2);
%         accR = accR + (a(j,15)^2 + a(j,16)^2 + a(j,17)^2)^(1/2);
%         vel = a(j,6); %vel + (a(j,6)^2 + a(j,7)^2 + a(j,8)^2)^(1/2);
%         acc = a(j,9); %acc = (a(j,9)^2 + a(j,10)^2 + a(j,11)^2)^(1/2);
    end;
    FT(count) = FTotal;
%     Acc(count) = acc;
%     AccA(count) = accA;
%     AccR(count) = accR;
    count = count + 1;
    box on;
    hold off;

    subplot(1,2,2);
    plot(1:size(FRes(idxcom),2),FA(idxcom),'.');
%     plot(1:size(FRes(idxcom),2),[FRes(idxcom); FA(idxcom); FR(idxcom)],'.');

    
%     axis([-10 10 -10 10 -10 10]);
    if (s>=1)
        axis([-s s -s s -s s]);     
    end;
%     subplot(1,2,2);
%     plot(i,[vel, acc],'.');
%     plot(i,[a(1,6), a(1,9)],'.');
%     plot(i,a(1,6));
%     hold on;

    pause(0.1);
end;

take = find((FR(idxcom)>0.9) & (FA(idxcom)>0.9));
subplot(1,2,1);
hold on;
p = plot3(a(idxcom(take),1),a(idxcom(take),2),a(idxcom(take),3),'*r');
set(p,'MarkerSize',50);


xlabel('x_1','FontSize',16);

% Create ylabel
ylabel('x_2','FontSize',16);

% Create zlabel
zlabel('x_3','FontSize',16);



% figure;
% subplot(2,2,1);
%     plot(Vel(2:si))
%     xlabel('Sum norm(Vel)','FontSize',16);
% subplot(2,2,2);
%     plot(Acc(2:si))
%     xlabel('Sum norm(Acc)','FontSize',16);
% subplot(2,2,3);
%     plot(AccA(2:si))
%     xlabel('Sum norm(Acc_A)','FontSize',16);
% subplot(2,2,4);
%     plot(AccR(2:si))
%     xlabel('Sum norm(Acc_R)','FontSize',16);

% nomef = sprintf('./Particle02/netGN_sample.time_final.par');
% a = load(nomef);
% [N C] = size(a);
% 
% for j=1:N
%         p = plot3(a(j,1),a(j,2),a(j,3),'.');
%         hold on;
%         set(p,'Color',cores(a(j,5),:), 'MarkerSize',50);
% end;
% box on;



