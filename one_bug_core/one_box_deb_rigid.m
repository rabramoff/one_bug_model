function dxdt=one_box_deb_rigid(t,x)
%Jinyun Tang: jinyuntang@lbl.gov
%one box model, the model includes leaching,
%mineral surface and exoenzyme dynamics
%it is for rigid microbe with zero reserve
global vid;
global par;


dxdt=zeros(size(x));


%enzymatic degradation of som
fsom=x(vid.ee)*x(vid.som)*par.vmax_ee/(par.ke_som+x(vid.som)+x(vid.ee)...
    +x(vid.mss)*par.ke_som/par.km_ee);

%transporter specific doc binding
fdoc=x(vid.doc)/(par.kb_doc+x(vid.doc)+...
    par.zb*x(vid.micb)+x(vid.mss)*par.kb_doc/par.km_doc);

%maximum specific enzyme production rate
pe_max=par.pro_ee0+par.pro_ee1*(1.0-fdoc)*fdoc;

%cell specific doc binding
fdoc=fdoc.*par.zb*par.vmax_micb;

%specific gross metabolites flux to support maintenance, growth and enzyme
%production
je=fdoc*par.Yld_micc;

%actual doc putake
fdoc=fdoc*x(vid.micb);

if(je>par.mr_micb)
    %positive growth
    mr=par.mr_micb;
    %net specific flux to support non-maintenance processes
    jc=je-mr;
    %maximum production rate for different processes
    gmax_vec=[par.gmax_micb,pe_max]';
    %yield rate for different processes
    yld_vec=[par.Yld_micb,par.Yld_ee]';
    %potential production rate
    gp=1./(1./gmax_vec+1./(jc.*yld_vec));
    %scaling factor due to substrate limitation
    scal=min([jc/(sum(gp./yld_vec)),1]);
    %actual population growth rate
    gB0=gp(1)*scal;
    %actual enzyme production rate
    pE=gp(2)*scal;
    %
    
else
    %actual maintenance rate
    mr=je;
    %no population growth
    gB0=0.;
    %no enzyme production
    pE=0.;
end



%compute specific mortality
decay_micb=par.decay_micb0*(1.0+par.mr_micb/(je+par.mr_micb)...
    *par.decay_micb1).*x(vid.micb)./(par.micb0+x(vid.micb));

gB=gB0-decay_micb;                 %net growth rate

%fprintf('gp=%f,gmx=%f,gB0=%f,gB=%f\n',gp(1),par.gmax_micb,gB0,gB);
%fprintf('gmax=%f,gB0=%f,decay_micb=%f\n',par.gmax_micb,gB0,decay_micb);
%reserve pool increment
Ja=par.Yld_micc*fdoc;              %the gross carbon assimilation


%export residual, currently is used to fuel excessive burning
residual=Ja-(mr+pE/par.Yld_ee+gB0/par.Yld_micb)*x(vid.micb);               

%fprintf('residual=%e\n',residual);
%get the enzyme loss due to degradation
dloss_ee=par.decay_ee*x(vid.ee);
%==========================================================================
%update state variables
dxdt(vid.som)=par.input_som-fsom+decay_micb*x(vid.micb)+...
    par.fee2som*dloss_ee;

dxdt(vid.doc)=par.input_doc+fsom-fdoc+...
    (1.0-par.fee2som)*dloss_ee-par.q_leach*x(vid.doc);


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
    if(abs(sum(dxdt(1:vid.co2))./(par.input_som+par.input_doc)-1)>1d-8)
        fprintf('diff=%e,inc=%e,input=%e\n',sum(dxdt(1:vid.co2))./(par.input_som+par.input_doc),...
            sum(dxdt(1:vid.co2)),par.input_som+par.input_doc);
        error('deb bad\n');
    
    end
else
    if(abs(sum(dxdt(1:vid.co2)))>1d-10)
        fprintf('zero input diff=%e\n',sum(dxdt(1:vid.co2)));
        error('deb bad\n');        
    end
end
