function par=par_scal_temp(par1, TEa, fref, Tinv)
%
%scale the temperature effects for different parameters
%
par=par1;


km_doc0=par.km_doc;
km_ee0=par.km_ee;
vmax_doc0=par.vmax_micb;
vmax_som0=par.vmax_ee;
mr0=par.mr_micb;
kb_doc0=par.kb_doc;
ke_som0=par.ke_som;
kappa0=par.kappa_micb;
decay_micb0=par.decay_micb0;

par.km_doc=km_doc0.*exp(-TEa.Ea_kmdoc.*Tinv);
par.km_ee=km_ee0.*exp(-TEa.Ea_kmee.*Tinv);
par.vmax_micb=vmax_doc0.*fref.*exp(-TEa.Ea_vbdoc.*Tinv);
par.vmax_ee=vmax_som0.*exp(-TEa.Ea_vesom.*Tinv);%*fref.
par.mr_micb=mr0.*exp(-TEa.Ea_mr.*Tinv);
par.kb_doc=kb_doc0.*exp(-TEa.Ea_kbdoc.*Tinv);
par.ke_som=ke_som0.*exp(-TEa.Ea_kesom.*Tinv);
par.kappa_micb=kappa0.*fref.*exp(-TEa.Ea_x.*Tinv);
par.decay_micb0=decay_micb0.*exp(-TEa.Ea_mr.*Tinv);

end