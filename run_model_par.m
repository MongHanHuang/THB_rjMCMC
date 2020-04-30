clear
close all

set(0,'defaultaxesfontsize',14)

% Add subdirectory to matlab path
addpath(genpath('./matcodes/'))

% Add fast marching toolbox location to matlab path
addpath(genpath('~/Documents/MATLAB/toolbox_fast_marching'))

clear Tmm ensemble

%    clearvars -except ShtDir maxiter hier datsav burn profiles topofiles Npro ipro ...
%        maxNv maxDpth maxVel maxZ minZ delta_X delta_Z v0 zz0 Ncol0 prior psig ...
%        fname DsaveNo
clearvars -except profiles topofiles Npro ipro

%% Model setup and Priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input files (p-wave arrival times and surface topography)
input_file = 'data/Input_test.mat'; % order in format: arrival time (sec), shot loc - receiver loc (m), receiver loc (m), shot loc (m)
elevation_file = 'data/elevation_20191216.txt'; % order in format: horizontal distance (m), elevation (m)

% Load Master data file
load (input_file);
Topo = load(elevation_file);

% Select how 'Filename' & 'Profile number' will be set:
Fname_Pronum_swtch = 1; % set to '1' to provide "File Name" & "Profile Number" here;
% Set 'Filename' & 'Profile number':
MODEL_FOLDER_NAME = 'MH234_NS_par'; % File Name; it will be under "model" folder
OUTPUT_FILE_NAME = 'result';  % output mat file name, saved under the model folder name 

PrgrsDsp = 1; % Display progress stats & figures, and save progress figures every this many 'data save steps' (variable = "datsave")


%%% Initial model setting %%%
% Choose depth range
maxZ = 90; % meters
minZ = 0; % meters

% Define grid size in X & Z dimensions
delta_X = 3; % m; horixontal delta distance between points (orig: 0.5 m) 
delta_Z = 3; % m; vertical (depth) distance between points (orig: 0.5 m)

% Set if maximum or minimum X is equal to zero (assuming 0 = first receiver location)
% set to '1' if the 1st shot location is < reiceiver locations;
% set to '2' if the 1st shot location is >= reiceiver locations;
% set to '0' for others
MaxMinX = 1; %%% this could be very confusing!!!

% Initial velocity and depth to top of layer
v0 = [300 500 1000 3000 4000 5000]'; % unit: meters per sec
zz0 = [0 10 20 50 80 maxZ]'; % depth of each layer in meters; need to be exactly the same number as v0
Ncol0 = 50; % hingeline number


%%% Prior model setting %%%
% Prior range for velocities and interface depths
prior.v1D = [300 5000]; % Lowest & Highest velocities allowed 
prior.h1D = [0 maxZ];
prior.Nlay = 70; % Maximum number of layers (has to be > 4)
prior.n = [-6 0.7]; % Range for prior noise estimate; in log10 base; unit in [ms]
prior.smooth_X = 4; % horizontal smoothing; model node distance has to be > smooth*delta_X [m]

% Set iteration and burn-in constraints
NumChain = 10; % number of chains you want to use (set this to be <= 4 if you are running this on a laptop)
maxiter = 100000;  % Maximum number of model 'iterations' (ex. 500000)
hier = 1;  % Hierachical switch? Otherwise = 0
datsav= 500;  % Save model to ensemble every this many 'iterations' (ex. 100)
burn_percet = 70; % Start doing stats after this many saved models; in percentage
total_save = maxiter/datsav;
burn = fix(total_save * NumChain * burn_percet / 100);  


%%% Proposal model parameters %%%
% "Proposal Sigmas"
% Standard deviations for proposing changes to model
psig.v1D = 200;% Proposal Sigma for changing velocity (m/s)
psig.h1D = 5; % Proposal Sigma for moving single layer element verticaly, or entire hingepoint laterally (# determines how much distance it can move)
psig.d1D = 100; % Proposal Sigma for the velocity variation of the added layer (m/s)
psig.n = 0.1; % Proposal Sigma for changing noise parameter (unit in natural log)

%% PROCESSING (No need to change below here)
% Rearrange the structure of Master (arrival time, receiver loc - shot loc , shot loc, receiver loc)
Master(:,4) = Master(:,3);
Master(:,3) = Master(:,2);
Master(:,2) = Master(:,4) - Master(:,3);

DsaveNo = 1; 

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

% Determine edges of domain
if MaxMinX == 1
    minX = min(Master(:,3));
    maxX = max(Master(:,3));
elseif MaxMinX == 2
    minX = min(Master(:,4));
    maxX = max(Master(:,4));
else
    maxX = max(XsrcUnique);
    minX = min(XsrcUnique);
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

% This -append save command added these variables to the '2D_params.mat' file,
%    but caused a conflict in the 'PlotEnsembleRef2D_verC.m' script where the ensemble is also set to 'X'.
% 'PlotEnsembleRef2D_verD.m' provides a bandaid to de-conflict,
%    but needs a variable naming convention update (see 'PlotEnsembleRef2D_verD.m' file).
save(['./models/' fname '/' fname '_2D_params.mat'],'Nx','Nz','X','Z','Xg','Zg','-append')

% Load topography
Topox = Topo(:,1)+minX;
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
model.v1D = v0;
model.xx = linspace(minX,maxX,Ncol0);
model.zz = repmat(zz0,1,Ncol0);

% Create model based on layer parameters and elevation
W = CreateModel2D(model,Xg,Zg,ElevFMM);

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

% data covariance vector
sigind = ones(Ndata,1);
Dsig = exp(model.xsig(sigind));

Nlay = length(model.v1D); %Initialize number of layers
Ncol = length(model.xx);

% Calculate error function and weighted error function
E=sum((T-Tm).^2./Dsig.^2);
E0=sum((T-Tm).^2);

% Initialize ensemble
ensemble(round(total_save)) = struct();

%% initialize parameters
% Initialize parameters for running mean and variance calc
Wmean = 0;
Wvar =0;

% Initialize counters
cnt=1;
kept=zeros(5,4);
tic

ENS = cell(NumChain,1);
CNT = cell(NumChain,1);

%% run iteration with number of chains
parfor i = 1:NumChain
    ENS{i} = ensemble;
    
    [ENS{i},CNT{i}] = run_MCMC(i,maxiter,hier,prior,model,psig,sigind,...
        Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,minZ,Tm,T,E,E0,...
        cnt,kept,Dsig,datsav,DsaveNo,ENS{i},delta_X);
    
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

%% Save products
saveDataOne(['./models/' fname '/' fname '_Final_2D_ensemble.mat'],ensemble)

%% Plot all of the figures
PlotEnsembleRef2D_verD

%% Save all of the products
save(['./models/' fname '/' outname '.mat']);

%% plot misfit time series
H=figure(1);clf;set(gca,'Fontsize',14);box on
for i = 1:NumChain
    tmp = [ENS{i}.E0]/Ndata; %misfit square
    NNEN = length(ENS{i});
%     plot((1:NNEN)*datsav,sqrt(tmp)*1e3,'linewidth',1.5); %plot in linear scale
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

%% subroutine running iterations
function [ensemble,cnt] = run_MCMC(ChainID,maxiter,hier,prior,model,psig,sigind,...
    Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,minZ,Tm,T,E,E0,...
    cnt,kept,Dsig,datsav,DsaveNo,ensemble,delta_X)

% Iterate THB algorithm!
for m = 1:maxiter
    
    % estimate current Nlay Ncol
    Nlay = length(model.v1D);
    Ncol = length(model.xx);

    % Choose an operation for updating model
    oper = RandomOperRef2D(hier,Nlay,Ncol,prior);
    
    % Update model
    [model2,delv] = UpdateMod2D(oper,model,psig,prior,delta_X);
    
    % Update hierarchical noise parameter
    Dsig2 = exp(model2.xsig(sigind));   
    
    if ~strcmp(oper(1:3),'noi')  %If not hierarchical step
        
        % Create new model
        W2 = CreateModel2D(model2,Xg,Zg,ElevFMM);
              
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
            
            % Run FMM on reduced domain
            D = perform_fast_marching(W2(:,mini:maxi), sp2, options);
                        
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
    keep = AcceptItRef(oper,dE,psig,delv,prior,Dsig,Dsig2,Nlay,Ncol);
    
    % Update counter
    kept(OpNumRef(oper),2) = kept(OpNumRef(oper),2)+1;
    kept(OpNumRef(oper),4) = kept(OpNumRef(oper),4)+1;
    
    % If we accept the new model, update values
    if keep>=rand(1)
        
        E=E2;
        E0=E02;
        Dsig=Dsig2;
        Tm=Tm2;
        model = model2;
        Nlay = length(model.v1D);
        Ncol = length(model.xx);
        kept(OpNumRef(oper),1) = kept(OpNumRef(oper),1)+1;
        kept(OpNumRef(oper),3) = kept(OpNumRef(oper),3)+1;
        
    end
    
    % Save model values and regular intervals
    if  mod(m,datsav)==0
        ensemble(cnt).v1D=model.v1D;
        ensemble(cnt).xx=model.xx;
        ensemble(cnt).zz=model.zz;
        ensemble(cnt).xsig=model.xsig;
        ensemble(cnt).E=E;
        ensemble(cnt).E0=E0;
        ensemble(cnt).Nlay=Nlay;
        ensemble(cnt).Ncol=Ncol;
        
        cnt=cnt+1;
                      
        if  ChainID == 1
            disp (['Save # ',num2str(DsaveNo),' of ',num2str(maxiter/datsav)]);
            DsaveNo = DsaveNo + 1;
        end
    end
end

end
