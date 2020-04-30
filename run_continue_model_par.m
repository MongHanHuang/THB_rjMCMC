clear
close all
% clc

file = 'models/topo3_B/result.1.mat';  % need to update by user
load(file)

output_name = 'result.2';  % need to update by user

close all

%% Users can redefine model parameters (highest/lowest velocity, proposal parameters)
% Add subdirectory to matlab path
addpath(genpath('./matcodes/'))

% Add fast marching toolbox location to matlab path
addpath(genpath('~/Documents/MATLAB/toolbox_fast_marching'))  % need to update by user

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set iteration and burn-in constraints
% Set iteration and burn-in constraints
maxiter = 2e5;  % Maximum number of model 'iterations' (ex. 500000)
hier = 1;  % Hierachical switch?
%total_save = 100;
%datsav= maxiter/total_save;  % Save model to ensemble every this many 'iterations' (ex. 100)
datsav = 500;
total_save = maxiter/datsav;

PrgrsDsp = 1; % Display progress stats & figures, and save progress figures every this many 'data save steps' (variable = "datsave")

% reassign prior and proposal sigmas
prior.v1D = [250 5000]; % Lowest & Highest velocities allowed 
prior.smooth_X = 4; % horizontal smoothing; model node distance has to be > smooth*delta_X [m]

% "Proposal Sigmas"
% Standard deviations for proposing changes to model
psig.d1D = 200;% 200; % Proposal Sigma for changing velocity
psig.h1D = 5; % Proposal Sigma for moving single layer element verticaly, or entire hingepoint laterally
psig.v1D = 100; %100; % Proposal Sigma for adding or removing a layer
psig.n = 0.1; %0.1; % Proposal Sigma for changing noise parameter (unit in natural log)

%% run iteration with number of chains; do not need to change below here
% re-set cnt
tic
cnt = length(ENS{1})+1;

clear model Dsig Tm E E0

% initiate model parameters from the last iteration
for i = 1:NumChain
    model{i} = ENS{i}(length(ENS{i}));
    Dsig{i} = exp(model{i}.xsig(sigind));
    W = CreateModel2D(model{i},Xg,Zg,ElevFMM);    
    for ishot = 1:Nshot
        start_point = [ZsrcFMM(ishot); XsrcFMM(ishot)];
        iind = find(SrcNumber==ishot);
        end_points = [ZrecFMM(iind)';XrecFMM(iind)'];
        
        % Find min and max X coordinates for this shot
        maxi = max([end_points(2,:) start_point(2,:)]);
        mini = min([end_points(2,:) start_point(2,:)]);
        % Define source and receiver locations on reduced domain
        sp2 = start_point;
        sp2(2,:) = sp2(2,:) - mini+1;
        ep2 = end_points;
        ep2(2,:) = ep2(2,:) - mini+1;
        
        TMPoptions.nb_iter_max = Inf;
        
        % Run FMM on reduced domain
        D = perform_fast_marching(W(:,mini:maxi), sp2, TMPoptions);
                
        for ii = 1:length(iind)
            Tm{i}(iind(ii),1) = (maxZ-minZ)*D(ep2(1,ii),ep2(2,ii))';
        end
    end  
    % Recalculate error function
    E{i}=sum((T-Tm{i}).^2./Dsig{i}.^2);
    E0{i}=sum((T-Tm{i}).^2);

end

parfor i = 1:NumChain        

    [ENS{i},CNT{i}] = run_MCMC(i,maxiter,hier,prior,model{i},psig,sigind,...
        Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,minZ,Tm{i},T,E{i},E0{i},...
        cnt,kept,Dsig{i},datsav,DsaveNo,ENS{i},delta_X);
    
end

tCPU = toc;
disp(['CPU time: ' num2str(tCPU/3600) ' hours']);

%% merge all of the chains back to ensemble
tic
clear ensemble
k = 0;
for j = 1:length(ENS{1})
    for i = 1:NumChain
        k = k+1;
        ensemble(k) = ENS{i}(j);
    end
end

burn = fix(burnPercent * k);

%% plot misfit time series
H=figure(16);clf;set(gca,'Fontsize',14);box on
for i = 1:NumChain
    tmp = [ENS{i}.E0]/Ndata;
    NNEN = length(ENS{i});
%     plot((1:NNEN)*datsav,sqrt(tmp)*1e3,'linewidth',1.5); %plot in linear scale
    semilogx((1:NNEN)*datsav,sqrt(tmp)*1e3,'linewidth',1.5); %plot in log scale
    hold on;
end
xlabel('Iterations','Fontsize',14);
ylabel('RMSE misfit (ms)','Fontsize',14);
grid on;

figfolder = ['./models/' fname '/figures/'];
p1=3.5;
set(H,'Units','Inches','Position',[1 1 1.5*p1 p1])
imfile = sprintf('%s/MulChain_%s_RawMisfit.pdf',figfolder,pronum);
set(gcf,'PaperPositionMode','auto')
set(gcf,'Units','Inches', 'PaperSize', [1.5*p1 p1]);
print(H,imfile,'-dpdf','-cmyk');

set(gcf,'PaperPositionMode','auto')
set(H,'Units','Inches','Position',[1 1 1.5*p1 p1])
set(gcf,'Units','Inches', 'PaperSize', [1.5*p1 p1]);
print(H,imfile,'-dpdf','-cmyk');

imfile = sprintf('%s/MulChain_%s_RawMisfit.fig',figfolder,pronum);
savefig(H,imfile,'compact');

%%
tCPU = toc;
disp(['CPU time: ' num2str(tCPU/3600) ' hours']);

%%
saveDataOne(['./models/' fname '/' fname '_Final_2D_ensemble.mat'],ensemble)

%% Plot all of the figures
PlotEnsembleRef2D_verD

%% save all of the products
save(['./models/' fname '/' output_name '.mat']);

%% subroutine
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
    
    model2.xx = model2.xx/max(model2.xx)*max(Xg(1,:)); %rescale
    
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
%         W=W2;
        Tm=Tm2;
        model = model2;
        Nlay = length(model.v1D);
        Ncol = length(model.xx);
        kept(OpNumRef(oper),1) = kept(OpNumRef(oper),1)+1;
        kept(OpNumRef(oper),3) = kept(OpNumRef(oper),3)+1;
        
    end
    
    % Save model values and regular intervals
    if  mod(m,datsav)==0
        %             cnt=cnt+1;
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
