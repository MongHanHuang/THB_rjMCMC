function W = CreateModel2D(model,Xg,Zg,Elev)
%%
AirV = 0.01;

xx = model.xx;
zz = model.zz;
v = model.v1D;

[xx,ixx]=unique(xx);
zz = zz(:,ixx);

Nzz = length(v);
Nxx = length(xx);

X = Xg(1,:);
Z = Zg(:,1);

Nx = length(X);
Nz = length(Z);


for ix = 1:Nxx-1
    xind(X>=xx(ix) & X<=xx(ix+1)) = ix;
end


for iz = 1:Nzz
    bd(iz,:) = zz(iz,xind) + (zz(iz,xind+1)-zz(iz,xind))./(xx(xind+1)-xx(xind)).*(X-xx(xind));
end


for izz = 1:Nzz-1
    for iz = 1:Nz
        zind(iz,Z(iz)>=bd(izz,:) & Z(iz)<=bd(izz+1,:)) = izz;
    end
end


for ix = 1:Nx
    W(:,ix) = v(zind(:,ix)) + (v(zind(:,ix)+1)-v(zind(:,ix)))./(bd(zind(:,ix)+1,ix)-bd(zind(:,ix),ix)).*(Z-bd(zind(:,ix),ix));
end


for ix = 1:Nx
    W(1:Nz<Elev(ix),ix) = AirV;
end

end
