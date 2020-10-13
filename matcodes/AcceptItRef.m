function keep = AcceptItRef(oper,dE,Dsig,Dsig2,Nuclei)

if ~strcmp(oper,'noise') %strcmp(oper(end-1:end),'1D')
    
    if strcmp(oper(1:3),'cha')
        
        keep = min(0, -dE/2);
        
    elseif strcmp(oper,'birthz')
        
        X = (Nuclei)/(Nuclei+1); 
        keep = min(0,log(X) -dE/2);
        
    elseif strcmp(oper,'deathz')
        
        X = (Nuclei)/(Nuclei-1); 
        keep = min(0,log(X) -dE/2);
        
    elseif strcmp(oper,'movez')

        keep = min(0, -dE/2);
                     
    elseif strcmp(oper,'swap')
        
        keep = min(0, -dE/2);
        
    end
    
elseif  strcmp(oper,'noise')
    
    X = sum(-log(Dsig2))-sum(-log(Dsig));
    keep = min(0,X -dE/2);
    
    
else
    disp('Nope, come on now')
end
