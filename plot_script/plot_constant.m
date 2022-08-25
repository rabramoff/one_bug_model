close all;
clear all;
clc;

%RZA updated these paths and deleted extra plots Aug 24, 22
mdir='./mat_files/constT/';
nfs={'one_box_deb_constT_Is1.314_Ic0.0876_dIs26.28_dIc1.752_plastic.mat',...
'one_box_deb_constT_Is1.314_Ic0.0876_dIs13.14_dIc0.876_plastic.mat'}; 

ax(1)=subplot(3,1,1);
ax(2)=subplot(3,1,2);
ax(3)=subplot(3,1,3);

load([mdir,nfs{1}]);
set(gcf,'CurrentAxes',ax(1));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.som]),2));
hold on;
set(gcf,'CurrentAxes',ax(2));
plot(TOUT_ctl(2:end)./365,diff(YOUT_ctl(:,vid.co2)));
hold on;
set(gcf,'CurrentAxes',ax(3));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.micb,vid.micc]),2));
hold on;

load([mdir,nfs{2}]);
set(gcf,'CurrentAxes',ax(1));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.som]),2),'r');
set(gcf,'CurrentAxes',ax(2));
plot(TOUT_ctl(2:end)./365,diff(YOUT_ctl(:,vid.co2)),'r');
set(gcf,'CurrentAxes',ax(3));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.micb,vid.micc]),2),'r');

set(ax,'Xlim',[100,200]);

