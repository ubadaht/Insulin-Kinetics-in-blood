function [t, y] = eul_imp(func,jacf,tspan,y0,h)
% 
%% Solve an IVP ODE by Euler's implicit method
% INPUT:
%       funz	IVP model function
%       jacf	Jacobian function
%       tspan [t0,tf]		integral interval
%       y0 		initial conditions
%       h 		discretization step
%
% OUTPUT
% t 	vector of discretization points in tspan interval
% y		vector of solution y(t) at t points
%
y0=y0';
t0=tspan(1);
tf=tspan(2);
t= t0:h:tf;
N=length(t);
I= eye(length(y0));
epsilon =1e-5; 
y(:,1)=y0;
for n=1:N-1
    y(:,n);
    y(:,n+1) =y(:,n);
 	A=I-h*feval(jacf,t(n+1),y(:,n+1));
 	b=y(:,n+1) -y(:,n) -h*feval(func,t(n+1),y(:,n+1));
 	u= A\b;
    y(:,n+1)=y(:,n+1) - u;
 	cont=1;
 	while(norm(u,'inf') >epsilon) & (cont<100)
  	 A=I-h*feval(jacf,t(n+1),y(:,n+1));
     b=y(:,n+1) -y(:,n) -h*feval(func,t(n+1),y(:,n+1));
     u= A\b;
  	 y(:,n+1)=y(:,n+1) - u;
  	 cont=cont+1;
 	end
 if cont ==100
   disp(['In' num2str(t(n)) ' iterations the method does not converge '])
 end
end
y=y';