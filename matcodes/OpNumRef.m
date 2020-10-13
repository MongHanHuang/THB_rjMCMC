function n=OpNumRef(oper)

if strcmp(oper(1:3),'cha')
    n=1;
elseif strcmp(oper(1:3),'bir')
    n=2;
elseif strcmp(oper(1:3),'dea')
    n=3;
elseif strcmp(oper(1:3),'mov')
    n = 4;
elseif strcmp(oper(1:3),'swa')
    n = 5;
elseif strcmp(oper,'noise')
    n = 6;
end

