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

%%
clear
close all

load('data/062421_colesville_line1_AW.mat'); % the master file with source/receiver info

plot_TT(Master);

%%
function plot_TT(Master)
% Rearrange the structure of Master (arrival time, receiver loc - shot loc , shot loc, receiver loc)
Master(:,4) = Master(:,3); % receiver location
Master(:,3) = Master(:,2); % source location
Master(:,2) = Master(:,4) - Master(:,3); % distance between source and receiver

Tmm = Master(:,1);
figure;
hold on
k = 0;
for i = 1:length(Master)-1
    if Master(i,3) == Master(i+1,3)
        k = k+1;
        tmpMaster(k,1) = Master(i,1);
        tmpMaster(k,2) = Master(i,4);
        mmod(k,1) = Tmm(i);
    else
        k = k+1;
        tmpMaster(k,1) = Master(i,1);
        tmpMaster(k,2) = Master(i,4);
        mmod(k,1) = Tmm(i);
        plot(tmpMaster(:,2),mmod*1e3,'linewidth',1); %plot model
        clear tmpMaster mmod
        k = 0;
    end
end
k = k+1;
tmpMaster(k,1) = Master(i+1,1);
tmpMaster(k,2) = Master(i+1,4);
mmod(k,1) = Tmm(i+1);

plot(tmpMaster(:,2),mmod*1e3,'linewidth',1); %plot model
xlim([min(Master(:,4)) max(Master(:,4))])
xlabel('Distance (m)');ylabel('Arrival time (ms)');
set(gca,'fontsize',13);
end
