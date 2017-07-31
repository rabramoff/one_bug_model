function dxdt=one_box_deb_plastic(t,x)
%Jinyun Tang: jinyuntang@lbl.gov
%one box model, the model includes leaching,
%mineral surface and exoenzyme dynamics
%it is for plastic microbe, that with explicit reserve
global vid;
global par;


dxdt=zeros(size(x));


%enzymatic degradation of som
fsom=x(vid.ee)*x(vid.som)*par.vmax_ee/(par.ke_som+x(vid.som)+x(vid.ee)...
    +x(vid.mss)*par.ke_som/par.km_ee);

%transporter specific doc binding
fdoc=x(vid.doc)/(par.kb_doc+x(vid.doc)+...
    par.zb*x(vid.micb)+x(vid.mss)*par.kb_doc/par.km_doc);

%
ec=x(vid.micc)*x(vid.micb)/(x(vid.micb)^2+1.d-20);

%maximum specific enzyme production rate
pe_max=par.pro_ee0+par.pro_ee1*(1.0-fdoc)*fdoc;


fdoc=fdoc.*par.zb*x(vid.micb)*par.vmax_micb;

%metabolites flux
je=par.kappa_micb*ec;

%determine what type of growth can be supported.    
isgrw=deb_microbe_init(par.mr_micb,[par.gmax_micb,pe_max]',[],[par.Yld_micb,par.Yld_ee]',je,ec);

switch isgrw
    case -1
        %no growth
        %there is no flux for growth, but penalty for mortality
        gB0=0.0;
        pE=0.0;
        mr=je;        
    case 0
        %feasible growth
        %normalized maximum growth rate > 0
        gmax_vec=[1,pe_max./par.gmax_micb]';
        %normalized maintenance respiration 
        mr_micb=par.mr_micb./par.gmax_micb;
        %normalized reserve export
        je=je./par.gmax_micb;    
        %yield rate
        yld_vec=[par.Yld_micb,par.Yld_ee]';
    
        %solve for the growth rate
        gp=deb_microbe(mr_micb,gmax_vec,0,yld_vec,je,ec);
 
        %specific population growth
        gB0=gp(1).*par.gmax_micb;
        %specific enzyme production rate
        pE=gp(2).*par.gmax_micb;
        %specific maintenance rate
        mr=par.mr_micb;          
    case 1
        %maximum growth
        gB0=par.gmax_micb;
        %enzyme production
        pE=pe_max;
        %maitenance
        mr=par.mr_micb;
end




%compute specific mortality
decay_micb=par.decay_micb0*(1.0+par.mr_micb/((par.kappa_micb-gB0)*ec+...
    par.mr_micb)*par.decay_micb1).*x(vid.micb)./(par.micb0+x(vid.micb));

gB=gB0-decay_micb;                 %net growth rate


%fprintf('gp=%f,gmx=%f,gB0=%f,gB=%f\n',gp(1),par.gmax_micb,gB0,gB);
%fprintf('gmax=%f,gB0=%f,decay_micb=%f\n',par.gmax_micb,gB0,decay_micb);
%reserve pool increment
Ja=par.Yld_micc*fdoc;              %the gross carbon assimilation


%get the reserve export for metabolism
Jc0=(par.kappa_micb-gB0)*x(vid.micc);

%actual reserve export including mortality
Jc=Jc0+decay_micb.*x(vid.micc);

%export residual, currently is used to fuel excessive burning
residual=Jc0-(mr+pE/par.Yld_ee+gB0/par.Yld_micb)*x(vid.micb);               

%fprintf('residual=%e\n',residual);
%get the enzyme loss
dloss_ee=par.decay_ee*x(vid.ee);
%==========================================================================
%update state variables
dxdt(vid.som)=par.input_som-fsom+decay_micb*x(vid.micb)+...
    par.fee2som*dloss_ee;

dxdt(vid.doc)=par.input_doc+fsom-fdoc+decay_micb*x(vid.micc)+...
    (1.0-par.fee2som)*dloss_ee-par.q_leach*x(vid.doc);

dxdt(vid.micc)=Ja-Jc;

dxdt(vid.micb)=gB*x(vid.micb);

dxdt(vid.ee)=pE*x(vid.micb)-dloss_ee;

dxdt(vid.co2)=fdoc-Ja+(mr+pE*(1d0/par.Yld_ee-1d0)+...
    gB0*(1d0/par.Yld_micb-1d0))*x(vid.micb)+residual;

dxdt(vid.cout)=dxdt(vid.co2)+pE.*x(vid.micb);

dxdt(vid.fdoc)=fdoc;

%if(t>39.9*365)
%fprintf('pE=%e\n',pE);
%end

if((par.input_som+par.input_doc)>0)
    if(abs(sum(dxdt(1:vid.co2))./(par.input_som+par.input_doc)-1)>1d-10)
        fprintf('diff=%e\n',sum(dxdt(1:vid.co2))./(par.input_som+par.input_doc));
        error('deb bad\n');
    
    end
else
    if(abs(sum(dxdt(1:vid.co2)))>1d-10)
        fprintf('zero input diff=%e\n',sum(dxdt(1:vid.co2)));
        error('deb bad\n');        
    end
end
