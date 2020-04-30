function oper = RandomOperRef2D(hier,Nlay,Ncol,prior)


if hier == 0
    
    nxt = randi(160);
    
    if nxt<21 && Nlay<prior.Nlay
        oper = 'birthz';
    elseif nxt<41 && Nlay>=4
        oper = 'deathz';
    elseif nxt<61
        oper = 'movez';
    elseif nxt<101
        oper = 'change1D';
    elseif nxt<121 && Ncol<prior.Ncol
        oper = 'birthx';
    elseif nxt<141 && Ncol>=4
        oper = 'deathx';
    elseif nxt<161
        oper = 'movex';
    end
    
elseif hier==1
    
    
    nxt = randi(120);
    
    if nxt<21 && Nlay<prior.Nlay
        oper = 'birthz';
    elseif nxt<41 && Nlay>=4 % make sure there are at least 4 lyr
        oper = 'deathz';
    elseif nxt<61
        oper = 'movez';
    elseif nxt<81
        oper = 'change1D';
    elseif nxt<101
        oper = 'movex';
    elseif nxt<121
        oper = 'noise';
    end
    
    
else
    disp('Is it hierarchical or not?')
end
