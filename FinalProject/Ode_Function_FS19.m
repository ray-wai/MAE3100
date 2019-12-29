%%Appendix F

%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Final Project ODE Function

%%This function just defines the ODE given in the project outline. It returns a vector 
%%that contains the velocity which is already given in the input yvec,
%%and the acceleration which is simply calculated from the ODE given in the
%%project outline. 

function [dr] = Ode_Function_FS19 (tvec, yvec)
mu=398716.708;
R=6378.137; %%given constants
J=0.0010826267;
x=yvec(1);
y=yvec(2); %%x, y, and z are the current position coordinates
z=yvec(3);
dx=yvec(4);
dy=yvec(5); %%dx, dy, dz are the current velocity vectors
dz=yvec(6);
r=sqrt((abs(x)).^2+(abs(y)).^2+(abs(z)).^2); %%distance equation

a=[(((-3.*J.*mu.*(R.^2).*x)./(2.*(r.^5))).*(1-((5.*z.^2)./(r.^2))));
   (((-3.*J.*mu.*(R.^2).*y)./(2.*(r.^5))).*(1-((5.*z.^2)./(r.^2)))); %%given 
   (((-3.*J.*mu.*(R.^2).*z)./(2.*(r.^5))).*(3-((5.*z.^2)./(r.^2))))];

d2x=((-mu.*x)./(r.^3))+a(1);
d2y=((-mu.*y)./(r.^3))+a(2); %%also given 
d2z=((-mu.*z)./(r.^3))+a(3);

dr=[dx;dy;dz;d2x;d2y;d2z]; %%this matrix will be multipled by step size to approximate the next set of position and velocity coordinates in the main program
end
