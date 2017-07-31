function par=set_par_default(scal)
if(nargin==0)
    scal=1;
end
par.input_som  =0.80;                                                      %input rate of somc         (gC/day)
par.input_doc  =0.20;                                                      %input rate of doc          (gC/day)
par.decay_mic  =1.5e-5.*scal;                                              %microbial death rate       1/day
par.decay_micb0=1.5d-5.*scal;                                              %microbial death rate       1/day
par.decay_micb1=0d0;
par.decay_ee   =1.e-5.*scal;                                               %enzyme decay rate          1/day
par.pro_ee     =1d-5.*scal;                                                %enzyme production rate     1/day/(g mic C)
par.gmax_micb  =1d-3.*scal;                                                %growth rate 1/day
par.pro_ee0    =5d-6.*scal;                                                %constintutive enzyme production rate     1/day
par.pro_ee1    =0d-6.*scal;                                                %inductive enzyme production rate  1/day
par.Yld_micb   =0.8;                                                       %growth efficiency of enzyme and microbes  (g mic C/g res C)
par.Yld_micc   =0.5;                                                       %assimilation efficiency from doc uptake   (g res C/g DOC C)
par.Yld_ee     =0.8;
par.mr_micb    =4d-4.*scal;                                                %microbial maintenance rate                (1/day)
par.kappa_micb =1d-1.*scal;                                                %reserve turnover rate                     (1/day)
par.q_leach    =0d-3;                                                      %leaching rate                             (1/day)
par.zb         =0.05;                                                      %scaling factor between transporter and microbial cell biomass
par.vmax_micb  =2d-2.*scal;                                                %maximum doc uptake rate                   (1/day)
par.vmax_ee    =2d-2.*scal;                                                %maximum som degradation rate              (1/day)
par.kb_doc     =1d0;                                                       %microbial doc affinity                    (g C/m3)
par.km_doc     =25d0;                                                      %adsorption surface doc affinity           (g C/m3)
par.ke_som     =20d1;                                                      %enzyme affinity to som                    (g ee c/m3)
par.km_ee      =5d1;                                                       %enzyme affinity for adsorptive surface    (g ee c/m3)
par.M1         =1d3;                                                       %abundance of mineral surface              (g C surface/m3)
par.fee2som    =0.2;                                                       %proportation of degraded exoenzyme into som (g som C/g ee C)
par.micb0      =1d-4;                                                      %half saturating microbial biomass for mortality (gC/m3), I have used an alternative value 1d-3 for sensitivity test
end