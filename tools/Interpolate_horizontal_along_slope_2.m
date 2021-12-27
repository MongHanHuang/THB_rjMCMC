clear
close all

ele_file = '../data/elevation_BORR_20210811';
master_file = '../data/BORR_20210811_shallow';
TD = 201; % TD = max tape distance
Td = 0; % min tap distance
delta_X = .2; % desired model grid size (m)

ele = load([ele_file '.txt']);
new_delta_X = 1/delta_X;

if Td < 0
    originPt = Td;
    ele(:,1) = ele(:,1) - originPt;
    ele(ele(:,1)<0,:) = [];
    TD = TD - Td;
end

eleP = ele(1,:);
eleP(1,3) = round(ele(1,1)*new_delta_X)/new_delta_X;
eleP(1,2) = ele(1,2);
for i = 2:length(ele)
    eleP(i,1) = eleP(i-1,1)+sqrt((ele(i,1)-ele(i-1,1))^2+(ele(i,2)-ele(i-1,2))^2);  % true horizontal distance
    eleP(i,2) = ele(i,2); % elevation
    eleP(i,3) = ele(i,1); % distance along slope
end

%%

SF = TD/max(eleP(:,1)); %SF: stretching factor; TD = max tape distance

eleP(:,1) = eleP(:,1)*SF;
eleP(:,3) = eleP(:,3)*SF;

% figure;
% % plot(eleP(:,1),eleP(:,2));
% hold on;plot(eleP(:,3),eleP(:,2),'.');
% axis image;

xx = 0:delta_X:TD;
zz = interp1(eleP(:,1),eleP(:,2),xx,'spline');
hori_xx = interp1(eleP(:,1),eleP(:,3),xx,'spline');
hori_xx = round(hori_xx*new_delta_X)/new_delta_X;

elev(:,1) = xx;
elev(:,2) = zz;
elev(:,3) = hori_xx;

%%
%% Redo Elevation File?
load([master_file '.mat']);

% set everything relative to the first location on the tape in the Master
if Td < 0
    Master(:,2:3) = Master(:,2:3) - Td;
end

for i = 1:length(Master)
    id = find(Master(i,2)==elev(:,1));
    Master(i,2) = elev(id,3);
    id = find(Master(i,3)==elev(:,1));
    Master(i,3) = elev(id,3);
end

new_elev(:,1) = elev(:,3);
new_elev(:,2) = elev(:,2);

figure;
subplot(2,1,1);
plot(elev(:,1),elev(:,2),'o');axis image;
subplot(2,1,2);plot(new_elev(:,1),new_elev(:,2),'o');axis image;

save([master_file '_modified.mat'],'Master');
save([ele_file '_modified.txt'],'new_elev','-ascii');
