function keep = AcceptItRef(oper,dE,psig,delv,prior,Dsig,Dsig2,Nlay,Ncol)


if ~strcmp(oper,'noise') %strcmp(oper(end-1:end),'1D')
    
    
    range=prior.v1D(2)-prior.v1D(1);
    
    if strcmp(oper(1:3),'cha')
        
        X = 1;
        keep = min(1,X*exp(-(dE)/2));
        
        
    elseif strcmp(oper,'birthz')
        
        
        %X = (Nlay+1)/(Nlay+2);
        X = psig.v1D*sqrt(2*pi)/range*exp(delv^2/2/psig.v1D^2);
        keep = min(1,X*exp(-(dE)/2));
        
        
    elseif strcmp(oper,'deathz')
        
        %X = (Nlay+1)/(Nlay);
        X = range/psig.v1D/sqrt(2*pi)*exp(-delv^2/2./psig.v1D^2);
        keep = min(1,X*exp(-(dE)/2));
        
        
    elseif strcmp(oper,'movez')
        
        X = 1;
        keep = min(1,X*exp(-(dE)/2));
        
            
    elseif strcmp(oper,'birthx')
        
        
        X = (Ncol+1)/(Ncol+2);
        keep = min(1,X*exp(-(dE)/2));
        
        
    elseif strcmp(oper,'deathx')
        
        X = (Ncol+1)/(Ncol);
        keep = min(1,X*exp(-(dE)/2));
        
        
    elseif strcmp(oper,'movex')
        
        X = 1;
        keep = min(1,X*exp(-(dE)/2));
        
        
    end
        
elseif  strcmp(oper,'noise')
    
    X = sum(-log(Dsig2))-sum(-log(Dsig));
    keep = min(1,exp(X-(dE)/2));
    
else
    display('Nope, come on now')
end