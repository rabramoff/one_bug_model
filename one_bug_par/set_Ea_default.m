function TEa=set_Ea_default()
%
%set default activation energy for different processes
%created by Jinyun Tang, Jan 27, 2014


rgas=8.31446;                      %universal gas constant, [J/K/mol]

TEa.Ea_vbdoc=45d3/rgas;            %doc uptake        vmax      [K]
TEa.Ea_vesom=45d3/rgas;            %soc degradation,  vmax         [K]
TEa.Ea_kbdoc=15d3/rgas;            %doc addsorption to microbe, affinity parameter [K]
TEa.Ea_kesom=15d3/rgas;            %som adsorption to enzyme, affinity parameter [K]
TEa.Ea_mr=60d3/rgas;               %maintenance  [K],  0.625ev, Brown 2004
TEa.Ea_x=60d3/rgas;                %reserve export [K], 0.625ev, Brown 2004

TEa.Ea_kmdoc=10d3/rgas;            %doc adsorption to mineral    [K]
TEa.Ea_kmee=10d3/rgas;             %enzyme adsorption to mineral [K]

end