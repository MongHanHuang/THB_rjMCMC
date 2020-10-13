%%  Mean model
close all

H=figure(1);clf;
imagesc(Xg(1,:)',Z(1,:),mask_Wmean);daspect([1 1 1]);
colormap(cc);colorbar;
caxis([250 5000]);
hold on;patch(polygonX,polygonY,'w');
plot(Topo(:,1),max(Topo(:,2)) - Topo(:,2)-dX/2,'k','linewidth',2);
plot(deep_ray(:,1),deep_ray(:,2),'w--','linewidth',2);
title(sprintf('Mean velocity (m/s) Profile %s',pronum),'Fontsize',14);
xlabel('Distance (m)');ylabel('Depth (m)')

% phone = unique(Master(:,4));
% for i = 1:length(phone)
%     for j = 1:length(Topo)
%         if (phone(i)-Topo(j,1)) == 0
%             plot(Topo(j,1),max(Topo(:,2)) - Topo(j,2),'ko'); %plot geophone location
%         end
%     end
% end

% Make profile
PT = ginput(1);
PT = round(PT/delta_X);

% PT = 240;
MEAN_STD = sqrt(Wvar/cnt); %mask_std;
% MEAN_STD = sqrt(Wvar/cnt);

PROF = mask_Wmean(:,PT(1));
% PROF = Wmean(:,PT(1));
PROF_STD = MEAN_STD(:,PT(1));
PROF_STD(PROF < 100) = [];
PROF(PROF < 100) = [];

plot([PT(1)*dX PT(1)*dX],[min(Z) max(Z)],'k--');

figure;
subplot(1,2,1);
hold on;
plot(PROF,(1:length(PROF))*delta_Z-delta_Z/2,'r','linewidth',3);
plot(PROF+PROF_STD*.5,(1:length(PROF))*delta_Z-delta_Z/2,'k--','linewidth',2);
plot(PROF-PROF_STD*.5,(1:length(PROF))*delta_Z-delta_Z/2,'k--','linewidth',2);
set(gca,'ydir','reverse');
xlabel('Velocity (m/s)');ylabel('Depth (m)');
title(['Vertical profile at ' num2str(PT(1)*dX) ' m']);
ylim([0 70]);xlim([0 5000]);
grid on;

subplot(1,2,2);
plot(gradient(PROF),(1:length(PROF))*delta_Z-delta_Z/2,'r','linewidth',3);
set(gca,'ydir','reverse');
xlabel('Velocity Gradient (m/s/m)');ylabel('Depth (m)');
title(['Vertical profile at ' num2str(PT(1)*dX) ' m']);
ylim([0 70]);
grid on;
