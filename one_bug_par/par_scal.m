function par1=par_scal(par,opt)
%scaling the rate parameters using the scaling factor searched from
%par_search
if(nargin==1)
    opt=1;
end
par1=par;
decay_mic0 =par.decay_mic;  
decay_ee0  =par.decay_ee;
pro_ee0    =par.pro_ee;
mr_micb0   =par.mr_micb;
kappa_micb0=par.kappa_micb;
vmax_micb0 =par.vmax_micb;
vmax_ee0   =par.vmax_ee;

decay_micb00=par.decay_micb0;
pro_ee00   =par.pro_ee0;
gmax_micb0 =par.gmax_micb;

%sc=[205.395497,953.561404,132.446849,10.732413,145.009704,458.397097,751.539950];
if(opt==0)
    %rigid microbe
    sc=[322.187769,320.863382,7.177276,2.681695,593.039365,412.434479,916.287907];
    %sc=[875.965050,611.427356,386.294367,57.653333,536.542,546.713826,120.666087];    
%sc=[707.504767,767.577053,141.304963,6.609155,377.900358,574.270844,589.337024];
%sc=[386.890771,568.956376,137.547903,62.477142,1.652079,594.278352,424.413762];
    sc=[875.965050,611.427356,386.294367,57.653333,536542d3,546.713826,120.666087];
else
    %plastic microbe
    sc=[875.965050,611.427356,386.294367,57.653333,0.536542,546.713826,120.666087];
end
par1.decay_mic = decay_mic0.*sc(1); 
par1.decay_ee  = decay_ee0.*sc(2);
par1.pro_ee    = pro_ee0.*sc(3);
par1.mr_micb   = mr_micb0.*sc(4);
par1.kappa_micb= kappa_micb0.*sc(5);
par1.vmax_micb = vmax_micb0.*sc(6);
par1.vmax_ee   = vmax_ee0.*sc(7);

par1.decay_micb0=decay_micb00.*sc(1);
par1.pro_ee0 = pro_ee00.*sc(3);
par1.gmax_micb=gmax_micb0.*sc(1);


end