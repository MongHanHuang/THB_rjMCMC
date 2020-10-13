function model2 = UpdateMod2D(oper,model,psig,prior)
model.success = 0;
model2 = model;

if ~strcmp(oper(1:3),'noi')
    
    v2 = model.v1D;
    xx2 = round(model.xx/prior.delta_X);
    zz2 = round(model.zz/prior.delta_Z);

    N = length(v2);
    
    if strcmp(oper(1:3),'cha')   %CHANGE V
        %%
        nind = randi(N);        
        newv = v2(nind) + psig.v1D*randn(1);
        
        if  newv >= min(prior.v1D) && newv <= max(prior.v1D)
            v2(nind) = newv;
            
            model2.success = 1;
        else
            return
        end

    elseif strcmp(oper,'swap')       %swap v of 2 cells
        %%
        nind1 = randi(N);
        nind2 = randi(N);
        
        tmpv = v2(nind1);
        v2(nind1) = v2(nind2);
        v2(nind2) = tmpv;
               
        model2.success = 1;
        
    elseif strcmp(oper,'movez')     %MOVE CELL (X,Z)
        %%
        
        nind = 4 + randi(N-4); %randomly choose a nucleus# > 4 (#1-4 are the 4 corners)
        
        delvZ = (psig.z1D*randn(1))/prior.delta_Z; % set a perturbation in vertical (psig.z1D)
        delvX = (psig.x1D*randn(1))/prior.delta_X; % set a perturbation in horizontal (psig.x1D)        
        newX = round(xx2(nind) + delvX);
        newZ = round(zz2(nind) + delvZ);
        
        while sum(abs(find(newX == xx2 & newZ == zz2))) ~= 0 %re-propose until no re-occupying
            delvZ = (psig.z1D*randn(1))/prior.delta_Z; % set a perturbation in vertical (psig.z1D)
            delvX = (psig.x1D*randn(1))/prior.delta_X; % set a perturbation in horizontal (psig.x1D)
            newX = round(xx2(nind) + delvX);
            newZ = round(zz2(nind) + delvZ);
        end
        
        if newZ > min(prior.h1D)/prior.delta_Z && newZ < max(prior.h1D)/prior.delta_Z ...
                && newX > min(prior.x1D)/prior.delta_X && newX < max(prior.x1D)/prior.delta_X
            
            xx2(nind) = newX;
            zz2(nind) = newZ;
            
            model2.success = 1;
        else
            return
        end
        
    elseif strcmp(oper,'birthz')    %ADD CELL (X,Z)
        if N >= prior.Nuclei 
            return
        end
        
        newv = min(prior.v1D) + rand*range(prior.v1D); %proposed velocity   
        newX = round((min(prior.x1D) + rand*range(prior.x1D))/prior.delta_X); %proposed x loc
        newZ = round((min(prior.h1D) + rand*range(prior.h1D))/prior.delta_Z); %proposed z loc
        
        while sum(abs(find(newX == xx2 & newZ == zz2))) ~= 0 %re-propose until no re-occupying
            newX = round((min(prior.x1D) + rand*range(prior.x1D))/prior.delta_X); %proposed x loc
            newZ = round((min(prior.h1D) + rand*range(prior.h1D))/prior.delta_Z); %proposed z loc
        end
        
        xx2(N+1) = newX;
        zz2(N+1) = newZ;
        v2(N+1) = newv;
        
        model2.success = 1;
                       
        %%
    elseif strcmp(oper,'deathz')     %REMOVE CELL (X,Z)
    
        nind = 4 + randi(N-4); %randomly choose a nucleus# > 4 (#1-4 are the 4 corners)
        
        v2(nind) = [];
        zz2(nind) = [];
        xx2(nind) = [];
        
        model2.success = 1;
            
        %%
    else
        disp('Thats not a thing')
    end
    
    model2.v1D = v2;
    model2.xx = xx2*prior.delta_X;
    model2.zz = zz2*prior.delta_Z;
     
else    %CHANGE NOISE (root mean square error)
    
    Dsig = model.xsig;
    delv=psig.n*randn(1);
    
    if Dsig+delv>=prior.n(1) && Dsig+delv<=prior.n(2)
        Dsig2=Dsig+delv;
    else
        Dsig2 = Dsig;
    end
    
    model2.xsig = Dsig2;
    model2.success = 1;
    
end

end
