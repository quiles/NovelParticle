clear all;

N = 10;
inicio = 0;
step = 1;
time = 100;
s = 4; % box size

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

for i=inicio:step:time-1
    nomef = sprintf('./net_sample.time_%d.par',i);
    a = load(nomef);
    [N C] = size(a);
    [i N]

    for j=1:N
        p = plot3(a(j,1),a(j,2),a(j,3),'.');
        hold on;

        
        set(p, 'MarkerSize',50);

%         set(p,'Color',cores(a(j,5),:), 'MarkerSize',50);
%         set(pb,'Color',cores(a(j,4),:), 'MarkerSize',10);
%         if (a(j,7) ~= 0)
%             minatt = min(a(:,7))-5;
% %             set(p,'Color',cores(a(j,4),:), 'MarkerSize',((a(j,7)-minatt)*2)^1.3);
%             set(p,'Color',cores(a(j,4),:), 'MarkerSize',50);
%         else
%             set(p,'Color',cores(a(j,4),:), 'MarkerSize',50);
%         end;
    end;
    
%     nomef = sprintf('./Particle02/data128/net_m0.46_r1.time_%d.centroids',i);
%     b = load(nomef);    
% 
%     [nC att] = size(b);
%     
%     for j=1:nC
%         q = plot3(b(j,1),b(j,2),b(j,3),'.');
%         set(q, 'MarkerSize',100,'Color',cores(b(j,4),:));
%     end;
    box on;
    hold off;
%     axis([-10 10 -10 10 -10 10]);
    axis([-s s -s s -s s]);

     
    pause(0.01);
end;

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


xlabel('x_1','FontSize',16);

% Create ylabel
ylabel('x_2','FontSize',16);

% Create zlabel
zlabel('x_3','FontSize',16);

