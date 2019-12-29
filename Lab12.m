%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 12, 11/14/19

%%Using Eulers and Runge-Kutta Methods to approximate solutions to an ODE

%%Eulers Method (First Order Runge-Kutta)
df=@(t) -2*t^3+12*t^2-20*t+8.5; %%given functions and constants
t=0:0.5:4;
h=0.5; 
x(1)=1; 

for i = 1 : length(t) - 1
   x(i+1) = x(i) + (df(t(i))).*h; %%Eulers method formula
end

fprintf('Eulers Method\n   Time(t)      x(t)\n'); %%displaying results to command window
disp([t' x']);
plot(t,x,'m');
xlabel('Time (t)'), ylabel('X(t)'); %%graphing results
hold on;

%%Ralstons Method (Second Order Runge-Kutta)
    %%a1 = 1/3, a2 = 2/3, p1 = q11 = 3/4 (given constants on lab doc)
y(1)=1;

for i = 1 : length(t) - 1
    k1 = df(t(i));
    k2 = df(t(i) + (3/4).*h); %%Ralstons method formulas
    y(i+1) = y(i) + ((1/3).*k1 + (2/3).*k2).*h; 
end

fprintf('\nRalstons Method\n   Time(t)      x(t)\n'); %%displaying and graphing results
disp([t' y']);
plot(t,y,'r');

%%Fourth Order Runge-Kutta Method
    %%p1 = p2 = 1/2, p3 = 1, a1 = a4 = 1/6, a2 = a3 = 1/3 (2/6) (given constants in class)
z(1)=1;

for i = 1 : length(t) - 1
    k1 = df(t(i));
    k2 = df(t(i) + (1/2).*h);
    k3 = df(t(i) + (1/2).*h); %%Runge-kutta 4th order formulas 
    k4 = df(t(i) + h);
    z(i+1) = z(i) + ((1/6).*k1 + (1/3).*k2 + (1/3).*k3 + (1/6).*k4).*h;
end

fprintf('\nRunge-Kutta 4th Order\n   Time(t)      x(t)\n'); %%displaying and graphing results
disp([t' z']);
plot(t,z,'g');

%%Exact Solution
f=@(t) -1/2.*t.^4+4.*t.^3-10.*t.^2+8.5.*t+1; %%solution to the ODE
s = f(t);
fprintf('\nExact Solution\n   Time(t)      x(t)\n'); %%displaying and graphing results
disp([t' s']);
t2=0:0.01:4;
plot(t2,f(t2),'b');
title('Using Eulers and Runge-Kutta Methods to Approximate Solutions to an ODE');
legend('Eulers','Ralstons','4th order RK','Exact Solution','Location','N');