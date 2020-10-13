function oper = RandomOperRef2D(hier,Nuclei,prior)


if hier == 0
    
    nxt = randi(100);
    
    if nxt<21 && Nuclei<prior.Nuclei
        oper = 'birthz';
    elseif nxt<41 && Nuclei>=4
        oper = 'deathz';
    elseif nxt<61
        oper = 'movez';
    elseif nxt<81
        oper = 'change1D';
    elseif nxt<101 
        oper = 'swap';
    end
    
elseif hier==1
    
    
    nxt = randi(120);
    
    if nxt<21 && Nuclei<prior.Nuclei
        oper = 'birthz';
    elseif nxt<41 && Nuclei>5 
        oper = 'deathz';
    elseif nxt<61
        oper = 'movez';
    elseif nxt<81
        oper = 'change1D';
    elseif nxt<101
        oper = 'noise';
    elseif nxt<121
        oper = 'swap';
    end
    
    
else
    disp('Is it hierarchical or not?')
end
