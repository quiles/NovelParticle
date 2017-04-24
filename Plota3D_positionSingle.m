clear all;


nomef = sprintf('networkF2.time_final.par');
% nomef = sprintf('./Particle02/posicao1.dat');
% nomef = sprintf('./Hie/network160.time_100.par');
nomec = sprintf('networkF2.cen');
a = load(nomef);
[N C] = size(a);

c = load(nomec);
[NC Dim] = size(c);
% b = load(nomec);
cores = jet(max(a(:,2))+1);

%         cent = plot3(b(:,1),b(:,2),b(:,3),'o');
%         hold on;
% 
%         set(cent, 'MarkerSize',60);

p = plot3(c(:,1),c(:,2),c(:,3),'r*');
set(p, 'MarkerSize',70);
hold on;

    for j=1:N
        p = plot3(a(j,4),a(j,5),a(j,6),'.');
        set(p, 'MarkerSize',50);
        set(p,'Color',cores(a(j,2),:), 'MarkerSize',50);
%         set(pb,'Color',cores(a(j,4),:), 'MarkerSize',10);
%         if (a(j,7) ~= 0)
%             minatt = min(a(:,7))-5;
% %             set(p,'Color',cores(a(j,4),:), 'MarkerSize',((a(j,7)-minatt)*2)^1.3);
%             set(p,'Color',cores(a(j,4),:), 'MarkerSize',50);
%         else
%             set(p,'Color',cores(a(j,4),:), 'MarkerSize',50);
%         end;
    end;

%     nomef = sprintf('./Particle01/exp/centroid%d.txt',i);
%     b = load(nomef);    
%     q = plot3(b(:,2),b(:,3),b(:,4),'o');
% 
%     set(q, 'MarkerSize',63);

    box on;
    hold off;
%     axis([-2 2 -2 2 -2 2]);

xlabel('x_1','FontSize',16);

% Create ylabel
ylabel('x_2','FontSize',16);

% Create zlabel
zlabel('x_3','FontSize',16);

