function [x,y] = euler_exp(funz,tspan,y0,h)

t0=tspan(1);
tf=tspan(2);
x=t0:h:tf;
n=length(x);
y(1,:)=y0;
for i=1:n-1
   f=funz(x(i),y(i,:));
   y(i+1,:)=y(i,:)+h*f';
end
end