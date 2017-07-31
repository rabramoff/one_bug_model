function one_box_deb_constT_2xInput_driver(Is,Ic)
%Ms: C protection equivalent of mineral surface, 500, 1000, 2000
%Tref: reference temperature, 290 or 280 K

%clc;
global vid;
global par;
opt=1;%plastic by default

[status,results]=system('pwd');
sstrs=strsplit(results,'/one_bug_model');

matfldir=[sstrs{1},'/one_bug_model/mat_files/constT2x/'];

system(['mkdir -p ', matfldir]);

%define vid
vid.som=1;
vid.doc=2;
vid.micc=3;
vid.micb=4;
vid.ee=5;
vid.co2=6;
vid.mss=7;
vid.fdoc=8;
vid.cout=9;


par=set_par_default();
par=par_scal(par,opt);
par.M1=5;
%frc_soc=.8;
par.input_som=Is./365;
par.input_doc=Ic./365;

iofile=[matfldir,'one_box_deb_constT2x_Is',num2str(Is),'_Ic',num2str(Ic),'_plastic.mat'];


%do spin up for 4000 years

[doc,micb,micc,ee,som,rco2,cue]=one_box_ss(par);

%define initial conditions
x0(vid.som)=som;
x0(vid.doc)=doc;
x0(vid.micc)=micc;
x0(vid.micb)=micb;
x0(vid.ee)=ee;
x0(vid.co2)=0;
x0(vid.mss)=par.M1;
x0(vid.fdoc)=0;
x0(vid.cout)=0;
%define the solver

xpar=[249.544170969785   5341.422691388677   5.617549086429];
TEa=set_Ea_default();

Tref=293.15;

par0=par;

T_fact0=fact_murphy(Tref,xpar);
Tfrz=293.15;


%do transient run 120 years
dt=1;
tend=365*200;


%introduce the temperature effect of mineral adsorption

A1=0;         %seasonal amplitude
A2=0;          %daily amplitude
w1=2.*pi/365;
w2=2.*pi;
kend=tend/dt;
flg=0;
TOUT_ctl=zeros(kend+1,1);
YOUT_ctl=zeros(kend+1,9);
YOUT_ctl(1,:)=x0;
TOUT_ctl(1,:)=0;
TEMP=zeros(kend,1);
for kk = 1 : kend
    t=(kk-0.5)*dt;
    xii=YOUT_ctl(kk,:);
    if(t>100*365)        
        par0.input_som=Is.*2./365;
        par0.input_doc=Ic.*2./365;
        if(flg==0)        
            xii=xii.*2;        
            xii(vid.mss)=xii(vid.mss)./2;
        end
        flg=1;
    else
        par0.input_som=Is./365;
        par0.input_doc=Ic./365;
    end
    temp=Tfrz+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);
    TEMP(kk)=temp;
    T_fact=fact_murphy(temp,xpar);
    fref=T_fact./T_fact0.*(temp./Tref);
    Tinv=1./temp-1./Tref;
    
    par=par_scal_temp(par0, TEa, fref, Tinv);
    
    YOUT_ctl(kk+1,:)=adptmbbks1(@one_box_deb,xii,TOUT_ctl(kk),dt);
    TOUT_ctl(kk+1)=TOUT_ctl(kk)+dt;
end



%plot(TEMP,diff(YOUT11(:,vid.co2))./dt,'.');
%hold on;
%plot(TEMP,diff(YOUT12(:,vid.co2))./dt,'r.');
%legend('Microbe+Mineral','Microbe only');
%set(gca,'FontSize',14);
%grid on;
%ylabel('g C day^-^1','FontSize',14);
%xlabel('Temperature (K)','FontSize',14);

save(iofile,'YOUT_ctl','TOUT_ctl','Is','Ic','vid');
end
