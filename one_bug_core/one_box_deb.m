function dxdt=one_box_deb(t,x)
%Jinyun Tang: jinyuntang@lbl.gov
%one box model, the model includes leaching,
%mineral surface and exoenzyme dynamics
global par;


if(par.kappa_micb>1d5)
   %rigid microbe
   dxdt=one_box_deb_rigid(t,x);
else
   %plastic microbe
   dxdt=one_box_deb_plastic(t,x);
end

end
