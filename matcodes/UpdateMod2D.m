function [model2,delv] = UpdateMod2D(oper,model,psig,prior,delta_X)
model2 = model;

if ~strcmp(oper(1:3),'noi') 
    
    v2 = model.v1D;
    xx2 = model.xx;   
    zz2 = model.zz;   
    
    Nz = length(v2);
    Nx = length(xx2);
    
    if strcmp(oper(1:3),'cha')   %CHANGE V
        %%
        nind = randi(length(v2));
        
        vtmp=[prior.v1D(1); v2; prior.v1D(2)];
               
        delv = psig.d1D*randn(1);
        delv = max([delv vtmp(nind)-vtmp(nind+1)]);
        delv = min([delv vtmp(nind+2)-vtmp(nind+1)]);
        
        v2(nind) = v2(nind) + delv;
                
    elseif strcmp(oper,'movez')     %MOVE LAYER (Z)
        
        nind = 1 + randi(Nz-2);
        nindx = randi(Nx);
        
        htmp=zz2(:,nindx);
        
        delv = psig.h1D*randn(1); %original code
        delv = max([delv htmp(nind-1) - htmp(nind)]);
        delv = min([delv htmp(nind+1) - htmp(nind)]);
        
        zz2(nind,nindx) = zz2(nind,nindx) + delv * .5 ; %make sure there is no overlap

%%
    elseif strcmp(oper,'birthz')    %ADD LAYER (Z)
        
        nind = 1 + randi(Nz-2);
        
        htmp = zz2;
        vtmp = v2;

        ioldv = nind; 

        delv = abs(psig.v1D*randn(1));        
        delv = max([delv vtmp(ioldv)-vtmp(ioldv+1)]);
        delv = min([delv vtmp(ioldv+1)-vtmp(ioldv)]);
        
        newv = vtmp(ioldv)+delv; % velocity of the new layer        
        
        for ix = 1:Nx % create topography of the new layer
            junk = rand*(htmp(nind+1,ix)-htmp(nind,ix)) + htmp(nind,ix);
            junk = max([min(prior.h1D) junk]);
            junk = min([max(prior.h1D) junk]);
            newh(1,ix) = junk;
        end

        if length(v2)<prior.Nlay 
            zz2 = [zz2; newh];
            v2 = [v2; newv];
            
            [~,ih]=sort(sum(zz2,2));
            zz2=zz2(ih,:);
            v2=v2(ih);
        else
            delv = 0;
        end 

        %%
    elseif strcmp(oper,'deathz')     %REMOVE LAYER (Z)
        
        nind = 1 + randi(length(v2)-2);
        
        delv = v2(nind-1) - v2(nind);
        
        v2(nind)=[];
        zz2(nind,:)=[];
      
%% 
    elseif strcmp(oper,'movex')     %MOVE HP (X)
        
        nindx = 1 + randi(Nx-2); % randomly choose an HP; make sure it is not at the L/R boundaries
                
        delv = psig.h1D*randn(1); % how much HP can move
        
        if delv>0 && (xx2(nindx+1) - (xx2(nindx) + delv)) > delta_X * prior.smooth_X % make sure not overlap & difference > delta_X (avoid sharp bd)
             xx2(nindx) = xx2(nindx) + delv;
        elseif delv<0 && ((xx2(nindx) + delv) - xx2(nindx-1)) > delta_X * prior.smooth_X % make sure not overlap & difference > delta_X
             xx2(nindx) = xx2(nindx) + delv;
        end

        
    else
        display('Thats not a thing')    
    end
        
    model2.v1D = v2;
    model2.xx = xx2;
    model2.zz = zz2;
        
else    %CHANGE NOISE (root mean square error)
   
    Dsig = model.xsig;
    delv=psig.n*randn(1);
     
    if Dsig+delv>=prior.n(1) && Dsig+delv<=prior.n(2)
        Dsig2=Dsig+delv;
    else
        Dsig2 = Dsig;
    end
    
    model2.xsig=Dsig2; 
    
end

end
