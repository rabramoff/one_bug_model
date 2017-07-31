close all;
clear all;
clc;

mdir='../mat_files/transtT/';
nfs={'one_box_deb_transtT_Is1.314_Ic0.0876_dIs0_dIc0.0876_plastic.mat',...% 2xdoc
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs1.314_dIc0_plastic.mat',... %2xsoc
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs1.314_dIc0.0876_plastic.mat'};        % 2xsoc & 2xdoc
figure;
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

load([mdir,nfs{3}]);
set(gcf,'CurrentAxes',ax(1));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.som]),2),'k');
set(gcf,'CurrentAxes',ax(2));
plot(TOUT_ctl(2:end)./365,diff(YOUT_ctl(:,vid.co2)),'k');
set(gcf,'CurrentAxes',ax(3));
plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.micb,vid.micc]),2),'k');

set(ax,'Xlim',[100,200]);



nfs={'one_box_deb_transtT_Is1.314_Ic0.0876_dIs1.314_dIc0.0876_plastic.mat',...    % 2xsoc & 2xdoc
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs6.57_dIc0.438_plastic.mat', ...     % 5x
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs13.14_dIc0.876_plastic.mat', ...    % 10x
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs19.71_dIc1.314_plastic.mat', ...    % 15x
'one_box_deb_transtT_Is1.314_Ic0.0876_dIs26.28_dIc1.752_plastic.mat'};        % 20x

fig=figure;
set(fig,'unit','normalized','position',[.1,.1,.6,.85]);
ax(1)=subplot(3,1,1);
ax(2)=subplot(3,1,2);
ax(3)=subplot(3,1,3);

linc={'b','g','r','k','c'};
for k = 1 : 5
    load([mdir,nfs{k}]);
    set(gcf,'CurrentAxes',ax(1));
    plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.som,vid.doc]),2),linc{k});
    hold on;
    set(gcf,'CurrentAxes',ax(2));
    plot(TOUT_ctl(2:end)./365,diff(YOUT_ctl(:,vid.co2)),linc{k});
    hold on;
    set(gcf,'CurrentAxes',ax(3));
    plot(TOUT_ctl./365,sum(YOUT_ctl(:,[vid.micb,vid.micc]),2),linc{k});
    hold on;
end
legend('2x','5x','10x','15x','20x');

set(ax,'Xlim',[100,150]);