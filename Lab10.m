%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 10, 11/1/19

%%Using finite difference method to approximate derivatives

y=[0 32 58 78 92 100]; %%given distance in km 
x=[0 25 50 75 100 125]; %%given time in s
[dy d2y]=diff12(x,y); %%dy = velocity, d2y = acceleration
disp('  Time      Velocity   Acceleration')
disp([x' dy' d2y']) %%displaying results from function call

function [dy d2y]=diff12(x,y)
n=length(x); %%number of points
h=(x(n)-x(1))./(n-1); %%distance b/t points

%forward - O(h^2)
i=1; %%first data point has no previous points, so forward must be used
dy(i)=(-y(i+2)+4.*y(i+1)-3.*y(i))./(2.*h);
d2y(i)=(-y(i+3)+4.*y(i+2)-5.*y(i+1)+2.*y(i))./(h.^2);

%backward - O(h^2)
i=n; %%last data point has no further points, so backward must be used
dy(n)=(3.*y(i)-4.*y(i-1)+y(i-2))/(2.*h);
d2y(n)=(2.*y(i)-5.*y(i-1)+4.*y(i-2)-y(i-3))/(h.^2);

%centered - O(h^2)
for i=2:n-1 %%using centered for all other points
   dy(i)=(y(i+1)-y(i-1))/(2.*h);
   d2y(i)=(y(i+1)-2.*y(i)+y(i-1))/(h.^2);
end
end 