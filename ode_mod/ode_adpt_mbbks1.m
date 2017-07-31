function [tout,yout]=ode_adpt_mbbks1(odefun,y0,dt,tspan)

nstep=fix((tspan(2)-tspan(1))/dt);
dtf=tspan(2)-tspan(1)-nstep.*dt;
if(dtf~=0.0)
    nstep=nstep+1;
end

tout=zeros(nstep+1,1);
yout=zeros(nstep+1,length(y0));

tout(1)=tspan(1);
yout(1,:)=y0;
for n = 1 : nstep-1
    yout(n+1,:)=adptmbbks1(odefun,yout(n,:),tout(n),dt);
    tout(n+1)=tout(n)+dt;   
    
end
%last step
if(dtf~=0.0)
    dt1=dt+dtf;
else
    dt1=dt;
end
n=nstep;
yout(n+1,:)=adptmbbks1(odefun,yout(n,:),tout(n),dt1);
tout(n+1)=tout(n)+dt1;        
end
