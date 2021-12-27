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
%        1 Apr, 2021 (update noise hyperparameter with noise varied with distance
%       14 Sep, 2021 (enable different interpolation methods; disable "swap")
%       19 Dec, 2021 (include raypath as an output in each iteration; add transparency to masked images)

%%  
clear
close all

% Add subdirectory to matlab path
addpath(genpath('./matcodes/'))

% Add fast marching toolbox location to matlab path
addpath(genpath('~/Documents/MATLAB/toolbox_fast_marching'))  % need to update by user

% Add fast marching toolbox location to matlab path (for Fast Marching by
% Dirk-Jan Kroon)
addpath(genpath('~/Documents/MATLAB/FastMarching_version3b'))

%%

MODEL_FOLDER_NAME = 'ModelA_test'; % File Name; it will be under "model" folder
Last_Model_File_Name = 'result';  % output mat file name, saved under the model folder name

%% Users can redefine model parameters (highest/lowest velocity, proposal parameters)

file = ['models/' MODEL_FOLDER_NAME '/' Last_Model_File_Name '.mat'];  % need to update by user
load(file)

%%%%%%% New output file name %%%%%%%%%%% 
new_output_name = 'result.1';  % need to update by user

% Set iteration and burn-in constraints
maxiter = 1e6;  % Maximum number of model 'iterations' (ex. 500000)

% reassign prior and proposal sigmas
prior.v1D = [250 6000]; % Lowest & Highest velocities allowed

burn_percent = 60; % Start doing stats after this many saved models; in percentage

% "Proposal Sigmas"
% Standard deviations for proposing changes to model
psig.v1D = 100; % Proposal Sigma for adding or removing a layer
psig.z1D = 1; % Proposal Sigma for moving single layer element verticaly, or entire hingepoint laterally
psig.x1D = 2; % Proposal Sigma for moving single layer element verticaly, or entire hingepoint laterally (in m; # determines how much distance it can move)
psig.n = 0.1; %0.1; % Proposal Sigma for changing noise parameter (unit in natural log)
psig.n2 = 0.01e-3; % Proposal Sigma for noise parameter slope with distance (unit in [ms/m])

%% run iteration with number of chains; do not need to change below here

output_name = new_output_name;

close all

tic
% reset cnt
cnt = length(ENS{1})+1;

clear model Dsig Tm E E0

% initiate model parameters from the last iteration
for i = 1:NumChain
    model{i} = ENS{i}(length(ENS{i}));
    Dsig{i} = exp(model{i}.xsig(sigind) .* model{i}.xsig2(sigind) .* Xsrc_rec);
    model{i}.success = 0;  % initialize success as 0
    F = scatteredInterpolant(model{i}.zz',model{i}.xx',model{i}.v1D');
    F.Method = interpMethod;
    W = F(Zg,Xg);
    [dimz , dimx] = size(W);
    for ii = 1:dimx
        W(1:dimz<ElevFMM(ii),ii) = 0.01;
    end
    W(W < 0.01) = 0.01;
    
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

% reset the B/D counting
BD = cell(NumChain,1);

% run iteration with number of chains
parfor i = 1:NumChain
    
    [ENS{i},CNT{i},BD{i},Ray{i}] = run_MCMC(i,maxiter,hier,prior,model{i},psig,...
        sigind,Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,...
        SrcNumber,maxZ,minZ,Tm{i},T,E{i},E0{i},cnt,kept,Dsig{i},datsav,...
        ENS{i},birthDeath,Xsrc_rec,delta_X,interpMethod,minX,Ray{i});
    
end

tCPU = toc;
disp(['CPU time: ' num2str(tCPU/3600) ' hours']);

%% merge all of the chains back to ensemble

clear ensemble
k = 0;
for j = 1:CNT{1}-1
    for i = 1:NumChain
        k = k+1;
        ensemble(k) = ENS{i}(j);
        raypath(k) = Ray{i}(j);
    end
end

BirthDeath = zeros(1,6);
for i = 1:NumChain
    BirthDeath(1) = BirthDeath(1) + BD{i}(1);
    BirthDeath(2) = BirthDeath(2) + BD{i}(2);
    BirthDeath(3) = BirthDeath(3) + BD{i}(3);
    BirthDeath(4) = BirthDeath(4) + BD{i}(4);
    BirthDeath(5) = BirthDeath(5) + BD{i}(5);
    BirthDeath(6) = BirthDeath(6) + BD{i}(6);
end

disp(['Birth: ' num2str(BirthDeath(1)) ' Death: ' num2str(BirthDeath(2))]);
disp(['Birth rate: ' num2str(BirthDeath(1)/maxiter/NumChain*100) ' %']);
disp(['Death rate: ' num2str(BirthDeath(2)/maxiter/NumChain*100) ' %']);
disp(['Noise rate: ' num2str(BirthDeath(3)/maxiter/NumChain*100) ' %']);
% disp(['Swap rate: ' num2str(BirthDeath(4)/maxiter/NumChain*100) ' %']);
disp(['Move rate: ' num2str(BirthDeath(5)/maxiter/NumChain*100) ' %']);
disp(['ChangeV rate: ' num2str(BirthDeath(6)/maxiter/NumChain*100) ' %']);

burn = fix(burn_percent/100 * k);

%% save all of the products
save(['./models/' fname '/' output_name '.mat'],'-v7.3');

%% Plot all of the figures
PlotEnsembleRef2D_verD

%% subroutine running iterations
function [ensemble,cnt,birthDeath,Ray] = run_MCMC(ChainID,maxiter,hier,prior,model,psig,sigind,...
    Xg,Zg,ElevFMM,Nshot,ZrecFMM,XrecFMM,ZsrcFMM,XsrcFMM,SrcNumber,maxZ,minZ,Tm,T,E,E0,...
    cnt,kept,Dsig,datsav,ensemble,birthDeath,Xsrc_rec,delta_X,interpMethod,minX,Ray)

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
        Dsig2 = exp(model2.xsig(sigind)) + model2.xsig2 .* Xsrc_rec; %updated Apr. 21, 2021
        
        if ~strcmp(oper(1:3),'noi')  %If not hierarchical step
            
            % Create new model
            F = scatteredInterpolant(model2.zz',model2.xx',model2.v1D');
            F.Method = interpMethod;
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
%                 D = msfm2d(W2(:,mini:maxi), sp2, false, false)*delta_X/(maxZ-minZ);
                                                     
                % Assign result to traveltime vector
                for ii = 1:length(iind)
                    Tm2(iind(ii),1) = (maxZ-minZ)*D(ep2(1,ii),ep2(2,ii))';
                end
            end
        else
            
            % Otherwise keep same traveltimes.
            Tm2= Tm;
        end
        
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
        ensemble(cnt).xsig2=model.xsig2;
        ensemble(cnt).E=E;
        ensemble(cnt).E0=E0;
        ensemble(cnt).Nuclei=Nuclei;
        
        % re-run fast marching and output raypath
        RayPath = [0 0];       
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
            D = perform_fast_marching(W2, start_point, options);
            gpath = compute_geodesic(D,end_points);
            
            for i = 1:length(gpath)
                [M,~] = size(RayPath);
                N = length(gpath{i});
                RayPath(M:M+N-1,1) = gpath{i}(2,:)*delta_X+minX;
                RayPath(M:M+N-1,2) = gpath{i}(1,:)*delta_X;
            end
        end

        [y,x] = size(W2);
        RayDen0 = zeros(y,x);
        
        MinX = min(RayPath(:,1))-delta_X;
        MinY = min(RayPath(:,2))-delta_X;
        
        for i = 1:length(RayPath)
            tmp_x = round((RayPath(i,1)-MinX)/delta_X);
            tmp_y = round((RayPath(i,2)-MinY)/delta_X);
            if tmp_x < x && tmp_y < y
                RayDen0(tmp_y,tmp_x) = RayDen0(tmp_y,tmp_x)+1;
            end
        end
        
        Ray(cnt).model = RayDen0; % raypath density of a model
        
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
