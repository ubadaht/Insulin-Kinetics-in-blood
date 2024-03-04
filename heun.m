function [t, y]=heun(funz,tspan,y0,h)

%% Solve an IVP ODE by Heun's explicit method
% INPUT:
%		funz	name of the file containing the f(x,y(x)) function
%		tspan   [t0,tf]	interval of interest for the solution
%		y0		Initial Condition
%		h		discretization step
% OUTPUT:
%		t		vector of discretization points in [t0,tf]
%		y		vector of solution values at points t 
%       FOR systems of m ODEs: each row contains a solution
%               y=[y_1(t);y_2(t);...;y_m(t)]
%%

t0=tspan(1);
tf=tspan(2);
t=t0:h:tf;
n=length(t);
y(1,:)=y0;
for i=1:n-1
   fi=funz(t(i),y(i,:));
   fip1=funz(t(i+1),y(i,:)+h*fi');
   y(i+1,:)=y(i,:)+0.5*h*(fi'+fip1');
end
