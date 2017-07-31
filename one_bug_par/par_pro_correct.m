function par=par_pro_correct(par,micb, micc,doc)
%do back calculation of the maximum growth rate and maximum enzyme production rate
%based on steady state pool sizes

if(micc==0.0)
    %doc uptake
    Fc=par.zb.*doc.*par.vmax_micb./(par.kb_doc+doc...
        +par.zb.*micb+par.M1.*par.kb_doc./par.km_doc);
    jc=Fc.*par.Yld_micc-par.mr_micb;
    %maximum enzyme production rate
    par.pro_ee0=1./(1./par.pro_ee-1./(jc.*par.Yld_ee));
    %maximum population growth rate
    par.gmax_micb=1./(1./par.decay_mic-1./(jc.*par.Yld_micb));
else
    %carbon export flux
    jc=(par.kappa_micb-par.decay_mic).*micc./micb-par.mr_micb;
    %maximum enzyme production rate
    par.pro_ee0=1./(1./par.pro_ee-1./(jc.*par.Yld_ee));
    %maximum population growth rate
    par.gmax_micb=1./(1./par.decay_mic-1./(jc.*par.Yld_micb));
end

end