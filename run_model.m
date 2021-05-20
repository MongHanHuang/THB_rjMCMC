% Transdimensional Hierarchical Bayesian (THB) framework with reversible-jump 
% Markov Chain Monte Carlo (rjMCMC) code for seismic refraction
% 
% Descreiption:
%         Please find the user guide for more details
%         
% Please cite this code as: Huang, M.-H., Hudson-Rasmussen, B., Burdick, S., 
% Lekic, V., Nelson, M.D., Fauria, K.E., and Schmerr, N., (2020), Bayesian 
% seismic refraction inversion for critical zone science and near-surface 
% applications, submitted to Geochem. Geophys. Geosys.
%     
% Authors: Mong-Han Huang (mhhuang@umd.edu)
%          Scott Burdick (sburdick@wayne.edu)
%          Vedran Lekic (ved@umd.edu)
%          Berit Hudson-Rasmussen (hudsonb@umd.edu)
% 
% Date: 11 Oct, 2020
%        2 Mar, 2021 (update input file names; enable FMM toolbox by Kroon 2021)

%%
clear
close all

set(0,'defaultaxesfontsize',14)

% Add subdirectory to matlab path
addpath(genpath('./matcodes/'))

% Add fast marching toolbox location to matlab path
addpath(genpath('~/Documents/MATLAB/toolbox_fast_marching'))

clear Tmm ensemble

clearvars -except profiles topofiles Npro ipro

%% Model setup and Priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input files (p-wave arrival times and surface topography)
input_file = 'data/ModelA.mat'; % order in format: arrival time (sec), shot loc (m), receiver loc (m)
elevation_file = 'data/elevation_ModelA.txt'; % order in format: horizontal distance (m), elevation (m)

% Load Master data file
load (input_file);
Topo = load(elevation_file);

% Select how 'Filename' & 'Profile number' will be set:
Fname_Pronum_swtch = 1; % set to '1' to provide "File Name" & "Profile Number" here;
% Set 'Filename' & 'Profile number':
MODEL_FOLDER_NAME = 'ModelA_test'; % File Name; it will be under "model" folder
OUTPUT_FILE_NAME = 'result';  % output mat file name, saved under the model folder name

%%% Initial model setting %%%
% Choose depth range
maxZ = 70; % meters
minZ = 0; % meters

% Define grid size in X & Z dimensions
delta_X = 1; % m; horixontal delta distance between points (this value should be smaller or equal to geophone interval)
delta_Z = delta_X; % m; vertical (depth) distance between points (suggest to be the same as delta_X)

% Initial velocity and depth to top of layer
v0 = [400 800 1500 2000 2500 3000 4000 5000]'; % unit: meters per sec
zz0 = [0 10 20 30 40 50 60 maxZ]'; % depth of each layer in meters; need to be exactly the same number as v0
Ncol0 = 30; % hingeline number for the initial model geometry

%%% Prior model setting %%%
% Prior range for velocities and interface depths
prior.v1D = [300 5000]; % Lowest & Highest velocities allowed (m/s)
prior.h1D = [0 maxZ];
prior.Nuclei = 500; % Maximum number of nuclei (has to be > 4)
prior.n = [-6 1]; % Range for prior noise estimate; in log10 base; unit in log [ms]

% Set iteration and burn-in constraints
NumChain = 4; % number of chains you want to use (set this to be <= 4 if you are running this on a laptop)
maxiter = 1e5;  % Maximum number of model 'iterations' (ex. 500000)
hier = 1;  % Hierachical switch? Otherwise = 0 (suggest to keep it to 1)
datsav= 200;  % Save model to ensemble every this many 'iterations' (ex. 100)
burn_percent = 60; % Start doing stats after this many saved models; in percentage

%%% Proposal model parameters %%%
% "Proposal Sigmas"
% Standard deviations for proposing changes to model
psig.v1D = 300; % Proposal Sigma for changing velocity (m/s)
psig.z1D = 10;  % Proposal Sigma for moving single nucleus verticaly (in m; # determines how much distance it can move)
psig.x1D = 30;  % Proposal Sigma for moving single nucleus laterally (in m; # determines how much distance it can move)
psig.n = 0.1;   % Proposal Sigma for changing noise parameter (unit in natural log)

%% PROCESSING (No need to change below here) 
% Rearrange the structure of Master (arrival time, receiver loc - shot loc , shot loc, receiver loc)
Master(:,4) = Master(:,3);
Master(:,3) = Master(:,2);
Master(:,2) = Master(:,4) - Master(:,3);
[~,id0] = sort(Master(:,4)); % sort the table with increasing geophone location orders
Master = Master(id0,:);

% Calculate the actual burn-in number
total_save = maxiter/datsav;
burn = fix(total_save * NumChain * burn_percent / 100);

if Fname_Pronum_swtch == 1
    fname = MODEL_FOLDER_NAME; % Filename
    outname = OUTPUT_FILE_NAME;
    pronum = fname; % Profile number
end

if ~exist(['./models/' fname], 'dir')
    mkdir(['./models/' fname]);
end

% convert noise prior from base 10 log to natural log
prior.n = log(10.^(prior.n));

% include delta_X/Z to prior
prior.delta_X = delta_X;
prior.delta_Z = delta_Z;

save(['./models/' fname '/' fname '_2D_params.mat'])

T= Master(:,1);  % Traveltimes
Xsrc = Master(:,3); % Source Locations
Xrec = Master(:,4); % Receiver Locations

% Find entries with zero offset or no traveltime pick
% and remove them from data
nono = abs(Master(:,2))==0 | isnan(T) | T==0 ;
Xsrc(nono)=[];
Xrec(nono)=[];
T(nono)=[];

% Find unique source locations
[XsrcUnique,~,SrcNumber] = unique(Xsrc);
Nshot = length(XsrcUnique);

% Determine edges of model domain (colume 3 is source, 4 is receiver
% locations)
if min(Master(:,3)) < min(Master(:,4))
    minX = min(Master(:,3));
elseif min(Master(:,3)) >= min(Master(:,4))
    minX = min(Master(:,4));
else  %this shouldn't exist
    minX = min(XsrcUnique);
end

if max(Master(:,3)) > max(Master(:,4))
    maxX = max(Master(:,3));
elseif max(Master(:,3)) <= max(Master(:,4))
    maxX = max(Master(:,4));
else %this shouldn't exit
    maxX = max(XsrcUnique);
end

% Define size of domain
Xrange = maxX-minX;
Zrange = maxZ-minZ;

% Number of nodes in domain
Nx = Xrange/delta_X+1;
Nz = Zrange/delta_Z+1;

% Create grid
X = linspace(minX,maxX,Nx);
Z = linspace(minZ,maxZ,Nz);
[Xg,Zg] = meshgrid(X,Z);
prior.x1D = [minX,maxX];

save(['./models/' fname '/' fname '_2D_params.mat'],'Nx','Nz','X','Z','Xg','Zg','-append')

% Load topography
Topox = Topo(:,1);
Topoz = max(Topo(:,2))-Topo(:,2);

% Interpolate topography onto X grid
Elev = (interp1(Topox,Topoz,X,'spline'))';

% Find top of topography for FMM
ElevFMM = round((Nz-1)*(Elev-minZ)/Zrange + 1);

% Define source and receiver locations for FMM algorithm
XsrcFMM = (Nx-1)*(XsrcUnique-minX)/Xrange + 1;
XrecFMM = (Nx-1)*(Xrec-minX)/Xrange + 1;
ZsrcFMM = ElevFMM(fix(XsrcFMM));
ZrecFMM = ElevFMM(fix(XrecFMM));

save(['./models/' fname '/' fname '_2D_grid.mat'],'X','Z','Xg','Zg','ElevFMM')

% Initial model
model0.v1D = v0;
model0.xx = linspace(minX,maxX,Ncol0);
model0.zz = repmat(zz0,1,Ncol0);

[dimz , dimx] = size(model0.zz);

k = 0;
% Generate nuclei for the 4 corners
for i = [1,dimx]
    for j = [1,dimz]
        k = k+1;
        model.xx(k) = model0.xx(i);
        model.zz(k) = model0.zz(j,i);
        model.v1D(k) = model0.v1D(j);
    end
end    
% Generate nuclei for the rest of the model
for i = 2:dimx-1
    for j = 2:dimz-1
        k = k+1;
        model.xx(k) = model0.xx(i);
        model.zz(k) = model0.zz(j,i);
        model.v1D(k) = model0.v1D(j);
    end
end

% Create model based on layer parameters and elevation
W = griddata(model.xx,model.zz,model.v1D,Xg,Zg);
[dimz , dimx] = size(W);
for i = 1:dimx
    W(1:dimz<ElevFMM(i),i) = 0.01;
end

% Initialize modeled traveltime vector
Tm = zeros(size(T));

% FMM computation, iterate over shot
for ishot = 1:Nshot
    
    % Define source location
    start_point = [ZsrcFMM(ishot); XsrcFMM(ishot)];
    
    options.nb_iter_max = Inf;
    
    % Calculate traveltime field D
    D = perform_fast_marching(W, start_point, options);
    
    % Which data are from this shot?
    iind = find(SrcNumber==ishot);
    
    % Assign result to correct spot in modeled traveltime vector
    for ii = 1:length(iind)
        Tm(iind(ii),1) = (maxZ-minZ)*D(ZrecFMM(iind(ii)),XrecFMM(iind(ii)))';
    end
    
end

%%
Ndata=length(T); % Number of picks
save(['./models/' fname '/' fname '_2D_params.mat'],'T','Ndata','-append');

Nsig = 1; % Number of noise variables

% Initialize hierarchical covariance param. based on initial data fit
model.xsig = log(sqrt(sum((T-Tm).^2)./Ndata)); % mean misfit in log_e scale

% Initialize model success index
model.success = 0;

% data covariance vector
sigind = ones(Ndata,1);
Dsig = exp(model.xsig(sigind));

% Calculate error function and weighted error function
E=sum((T-Tm).^2./Dsig.^2);
E0=sum((T-Tm).^2);

% Initialize ensemble
ensemble(round(total_save)) = struct();

%% initialize parameters
% Initialize parameters for running mean and variance calc
Wmean = 0;
Wvar = 0;

% Initialize counters
cnt=1;
kept=zeros(7,4);
birthDeath = zeros(1,6);

ENS = cell(NumChain,1);
CNT = cell(NumChain,1);
BD = cell(NumChain,1);

%% run iteration with number of chains
tic
parfor i = 1:NumChain
    ENS{i} = ensemble;
    
    [ENS{i},CNT{i},BD{i}] = run_MCMC(i,maxiter,hier,prior,model,psig,sigind,...
        Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,...
        minZ,Tm,T,E,E0,cnt,kept,Dsig,datsav,ENS{i},birthDeath,delta_X);
    
end

t = toc;

disp(['Total CPU time: ' num2str(t/3600) ' hours']);

%% merge all of the chains back to ensemble
clear ensemble
k = 0;
for j = 1:CNT{1}-1
    for i = 1:NumChain
        k = k+1;
        ensemble(k) = ENS{i}(j);
    end
end

%Compile summary of the steps
BirthDeath = zeros(1,6);
for i = 1:NumChain
    BirthDeath(1) = BirthDeath(1) + BD{i}(1);
    BirthDeath(2) = BirthDeath(2) + BD{i}(2);    
    BirthDeath(3) = BirthDeath(3) + BD{i}(3);    
    BirthDeath(4) = BirthDeath(4) + BD{i}(4);    
    BirthDeath(5) = BirthDeath(5) + BD{i}(5);    
    BirthDeath(6) = BirthDeath(6) + BD{i}(6);    
end

% disp(['Birth: ' num2str(BirthDeath(1)) ' Death: ' num2str(BirthDeath(2))]);
disp(['Birth rate: ' num2str(BirthDeath(1)/maxiter/NumChain*100) ' %']);
disp(['Death rate: ' num2str(BirthDeath(2)/maxiter/NumChain*100) ' %']);
disp(['Noise rate: ' num2str(BirthDeath(3)/maxiter/NumChain*100) ' %']);
disp(['SwapV rate: ' num2str(BirthDeath(4)/maxiter/NumChain*100) ' %']);
disp(['Move rate: ' num2str(BirthDeath(5)/maxiter/NumChain*100) ' %']);
disp(['ChangeV rate: ' num2str(BirthDeath(6)/maxiter/NumChain*100) ' %']);

%% Save products
saveDataOne(['./models/' fname '/' fname '_Final_2D_ensemble.mat'],ensemble)

%% Save all of the products
save(['./models/' fname '/' outname '.mat']);

PlotEnsembleRef2D_verD

%% plot misfit time series
H=figure(1);clf;set(gca,'Fontsize',14);box on
for i = 1:NumChain
    tmp = [ENS{i}.E0]/Ndata; %misfit square
    NNEN = length(ENS{i});
    semilogx((1:NNEN)*datsav,sqrt(tmp)*1e3,'linewidth',1.5); %plot in log scale
    hold on;
end
xlabel('Iterations','Fontsize',14);
ylabel('RMSE misfit (ms)','Fontsize',14);

figfolder = ['./models/' fname '/figures/'];
p1=3.5;
set(H,'Units','Inches','Position',[1 1 1.5*p1 p1])
imfile = sprintf('%s/MulChain_%s_RawMisfit.pdf',figfolder,pronum);
set(gcf,'PaperPositionMode','auto')
set(gcf,'Units','Inches', 'PaperSize', [1.5*p1 p1]);
print(H,imfile,'-dpdf','-cmyk');

imfile = sprintf('%s/MulChain_%s_RawMisfit.fig',figfolder,pronum);
savefig(H,imfile,'compact');

%% subroutine running iterations
function [ensemble,cnt,birthDeath] = run_MCMC(ChainID,maxiter,hier,prior,model,psig,sigind,...
    Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,minZ,Tm,T,E,E0,...
    cnt,kept,Dsig,datsav,ensemble,birthDeath,delta_X)

LabelPercent = maxiter/100*5; %show status with 5% increment
tic

% Iterate THB algorithm!
for m = 1:maxiter
    
    % estimate current Nuclei
    Nuclei = length(model.v1D);
    
    % Choose an operation for updating model
    oper = RandomOperRef2D(hier,Nuclei,prior);
    
    % Update model
    model2 = UpdateMod2D(oper,model,psig,prior);
    if model2.success == 1
        
        % Update hierarchical noise parameter
        Dsig2 = exp(model2.xsig(sigind));
        
        if ~strcmp(oper(1:3),'noi')  %If not hierarchical step
            
            % Create new model
            F = scatteredInterpolant(model2.zz',model2.xx',model2.v1D');
            F.Method = 'linear';
            W2 = F(Zg,Xg);
            [dimz , dimx] = size(W2);
            for i = 1:dimx
                W2(1:dimz<ElevFMM(i),i) = 0.01;
            end
            W2(isnan(W2)) = 0.01; % this is unlikely
            W2(W2 < 0.01) = 0.01; % this is unlikely
            
            % FMM computation
            for ishot = 1:Nshot
                
                % Define source location
                start_point = [ZsrcFMM(ishot); XsrcFMM(ishot)];
                
                % Which data are from this shot?
                iind = find(SrcNumber==ishot);
                
                % Define receiver locations
                end_points = [ZrecFMM(iind)';XrecFMM(iind)'];
                
                % Find min and max X coordinates for this shot
                maxi = max([end_points(2,:) start_point(2,:)]);
                mini = min([end_points(2,:) start_point(2,:)]);
                
                % Define source and receiver locations on reduced domain
                sp2 = start_point;
                sp2(2,:) = sp2(2,:) - mini+1;
                
                ep2 = end_points;
                ep2(2,:) = ep2(2,:) - mini+1;
                               
                options.nb_iter_max = Inf;
                
                %%% Run FMM on reduced domain using Peyre (2020) FMM toolbox (uncomment the line below if you want to use this toolbox)
                D = perform_fast_marching(W2(:,mini:maxi), sp2, options);

                %%% Run FMM on reduced domain using Kroon (2021) FMM toolbox (uncomment the line below if you want to use this toolbox)
                % D = msfm2d(W2(:,mini:maxi), sp2, true, true)*delta_X/(maxZ-minZ));
                                                     
                % Assign result to traveltime vector
                for ii = 1:length(iind)
                    Tm2(iind(ii),1) = (maxZ-minZ)*D(ep2(1,ii),ep2(2,ii))';
                end
            end
        else
            
            % Otherwise keep same traveltimes.
            Tm2= Tm;
        end
        
        % Recalculate error function
        E2=sum((T-Tm2).^2./Dsig2.^2);
        E02=sum((T-Tm2).^2);
        
        dE=E2-E; % Change in error
        
        % Decide whether to accept or reject model
        keep = AcceptItRef(oper,dE,Dsig,Dsig2,Nuclei);
        
        % Update counter
        kept(OpNumRef(oper),2) = kept(OpNumRef(oper),2)+1;
        kept(OpNumRef(oper),4) = kept(OpNumRef(oper),4)+1;
        
        % If we accept the new model, update values
        if keep >= log(rand(1))            
            E=E2;
            E0=E02;
            Dsig=Dsig2;
            Tm=Tm2;
            model = model2;
            Nuclei = length(model.v1D);
            kept(OpNumRef(oper),1) = kept(OpNumRef(oper),1)+1;
            kept(OpNumRef(oper),3) = kept(OpNumRef(oper),3)+1;
            if strcmp(oper,'birthz')
                birthDeath(1) = birthDeath(1)+1;
            elseif strcmp(oper,'deathz')
                birthDeath(2) = birthDeath(2)+1;
            elseif strcmp(oper(1:3),'noi')
                birthDeath(3) = birthDeath(3)+1;
            elseif strcmp(oper,'swap')
                birthDeath(4) = birthDeath(4)+1;
            elseif strcmp(oper,'movez')
                birthDeath(5) = birthDeath(5)+1;
            elseif strcmp(oper(1:3),'cha')
                birthDeath(6) = birthDeath(6)+1;
            end            
        end
    end
    
    % Save model values and regular intervals
    if  mod(m,datsav)==0
        ensemble(cnt).v1D=model.v1D;
        ensemble(cnt).xx=model.xx;
        ensemble(cnt).zz=model.zz;
        ensemble(cnt).xsig=model.xsig;
        ensemble(cnt).E=E;
        ensemble(cnt).E0=E0;
        ensemble(cnt).Nuclei=Nuclei;
        
        cnt=cnt+1;
    end
    
    if mod(m,LabelPercent) == 0 && ChainID == 1
        t_2 = toc;
        work_percent = m/LabelPercent*5;
        disp(['Working... (' num2str(work_percent) '% done; ~' ...
            num2str(round((t_2)/60*(100-work_percent)/work_percent)) ' min(s) left)']);
    end
    
end
end
