function one_bug_transT_runs(Tref,opt,Yld_x)
%
%Tref, reference temperature
%opt, 1(plastic)/0(rigid) microbe
%Yld_x, substrate yield rate

Ms=[500,1000,1500,2000];
Ems=[0,10d3,40d3];
Esc=[25d3,45d3];

for j1 = 1 : length(Ms)
    for j2 = 1 : length(Ems)
        for j3= 1: length(Esc)
            one_box_deb_transientT_driver(Ms(j1),Tref,Ems(j2),Esc(j3),Yld_x,opt);
        end
    end
end

end