% Transdimensional Hierarchical Bayesian (THB) framework with reversible-jump
% Markov Chain Monte Carlo (rjMCMC) code for seismic refraction
%
% Descreiption:
%         Please find the user guide (THB2D MCMC User Guide v4.pdf) for more details
%
% Please cite this code as:
% Huang, M.-H., Hudson-Rasmussen, B., Burdick, S., Lekic, V., Nelson, M.D.,
%     Fauria, K.E., and Schmerr, N., (2020), Bayesian seismic refraction
%     inversion for critical zone science and near-surface applications,
%     Geopchemistry, Geophysics, Geosystems, 22, e2020GC009172.
%     https://doi.org/10.1029/2020GC009172
%
% Authors: Mong-Han Huang (mhhuang@umd.edu)
%          Scott Burdick (sburdick@wayne.edu)
%          Vedran Lekic (ved@umd.edu)
%          Berit Hudson-Rasmussen (hudsonb@umd.edu)
%
% Date: 11 Oct, 2020
%        2 Mar, 2021 (update input file names; enable FMM toolbox by Kroon 2021)
%        1 Apr, 2021 (update noise hyperparameter allowing noise to vary with distance
%       14 Sep, 2021 (enable different interpolation methods)

%%  Mean model
close all

MaxDep = 30;
Mean_Velocity = mask_Wmean; % mask_Wmean2;%mask_W_intp;

meter2ft = 3.28;

% Add toolbox location to matlab path
addpath(genpath('~/Documents/MATLAB/crameri_v1.05')) % for color scale

cc = [1 1 1; crameri('roma')];
cstd = [1 1 1; crameri('lajolla')];
ccc = [1 1 1;(parula(1024))];

%%% calculate mean velocity relative to ground surface
[y,x] = size(Wmean);
for i = 1:x
    id = ElevFMM(i);
    Wmean_ths(1:y-id+1,i) = Wmean(id:y,i);
    mask_Wmean_ths(1:y-id+1,i) = Mean_Velocity(id:y,i);
    mask_dWmean_ths(1:y-id+1,i) = Mean_Velocity(id:y,i);
end


H=figure;
clf;
subplot(2,1,1);
imagesc(Xg(1,:)',Z(1,:),Mean_Velocity);daspect([1 1 1]);
colormap(cc);colorbar;
caxis([250 4000]);
hold on;patch(polygonX,polygonY,'w');
% plot(Topo(:,1),max(Topo(:,2)) - Topo(:,2)-dX/2,'k','linewidth',2);
title(sprintf('Mean velocity (m/s) Profile %s',pronum),'Fontsize',14);
xlabel('Distance (m)');ylabel('Depth (m)')

subplot(2,1,2);
hold on
imagesc(Xg(1,:)',Z,mask_Wmean_ths);axis image;colormap(cc);colorbar;
caxis([250 4000]);xlabel('Distance (m)');ylabel('Depth (m)')
set(gca,'ydir','reverse');

%{
asdfasdfadsf
%}
% phone = unique(Master(:,4));
% for i = 1:length(phone)
%     for j = 1:length(Topo)
%         if (phone(i)-Topo(j,1)) == 0
%             plot(Topo(j,1),max(Topo(:,2)) - Topo(j,2),'ko'); %plot geophone location
%         end
%     end
% end

% Make profile
minX = min(Xg(1,:));
PT = ginput(1);  %distance
PT(2) = [];
PT = round((PT-minX)/delta_X);  %convert to index
skip = ElevFMM(PT); %should be index stored in ElevFMM, not distance 

% PT = 240;
MEAN_STD = sqrt(Wvar/cnt); %mask_std;
% MEAN_STD = sqrt(Wvar/cnt);

PROF = Mean_Velocity(:,PT(1));
% PROF = Wmean(:,PT(1));
PROF_STD = MEAN_STD(:,PT(1));
PROF_STD(1:skip) = [];
PROF(1:skip) = [];

subplot(2,1,1);
plot([PT(1)*dX+minX PT(1)*dX+minX],[min(Z) max(Z)],'k--');
subplot(2,1,2);
plot([PT(1)*dX+minX PT(1)*dX+minX],[min(Z) max(Z)],'k--');

figure;
subplot(1,2,1);
hold on;
plot(PROF,(1:length(PROF))*delta_Z-delta_Z/2,'r','linewidth',3);
plot(PROF+PROF_STD*.5,(1:length(PROF))*delta_Z-delta_Z/2,'k--','linewidth',2);
plot(PROF-PROF_STD*.5,(1:length(PROF))*delta_Z-delta_Z/2,'k--','linewidth',2);
set(gca,'ydir','reverse');
xlabel('Velocity (m/s)');ylabel('Depth (m)');
title(['Vertical profile at ' num2str(PT(1)*dX) ' m']);
ylim([0 MaxDep]);
xlim([0 5000]);
grid on;

subplot(1,2,2);
plot(gradient(PROF),(1:length(PROF))*delta_Z-delta_Z/2,'r','linewidth',3);
set(gca,'ydir','reverse');
xlabel('Velocity Gradient (m/s/m)');ylabel('Depth (m)');
title(['Vertical profile at ' num2str(PT(1)*dX) ' m']);
ylim([0 MaxDep]);
grid on;
