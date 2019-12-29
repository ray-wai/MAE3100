%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, HW9, 11/20/19

%%%%% 22.2 %%%%%
dydx = @(x,y) (1+2.*x).*sqrt(y);
h=.25;
x=[0:.25:1];
y(1)=1;

%%Part A, exact Solution
u=[0:.01:1];
fx = ((u.^2 + u + 2)./2).^2;
figure(1);
plot(u,fx, 'm');
hold on;

% Euler's Method
for i=1:length(x)-1;
    y(i+1)=y(i)+dydx(x(i),y(i))*h;
end

uEuler=[x;y];
plot(x, uEuler(2,:),'b');

%%Huens w/o iteration (2nd order RK with a1 = a2 = 1/2, p1 = q11 = 1)
for i=1:length(x)-1;
    k1(i) = dydx(x(i),y(i));
    k2(i) = dydx(x(i)+h, y(i)+k1(i).*h);
    y(i+1)=y(i)+((1/2).*k1(i)+(1/2).*k2(i))*h;
end

uHuen=[x;y];
plot(x, uHuen(2,:),'r--');

%%Ralstons (2nd order RK with a1 = 1/4, a2 = 3/4, p1 = q11 = 2/3)
for i=1:length(x)-1;
    k1(i) = dydx(x(i),y(i));
    k2(i) = dydx(x(i)+(2/3).*h, y(i)+k1(i).*(2/3).*h);
    y(i+1)=y(i)+((1/4).*k1(i)+(3/4).*k2(i))*h;
end

uRalston=[x;y];
plot(x, uRalston(2,:),'g');

%%RK Fourth Order
for i=1:length(x)-1;
    k1(i) = dydx(x(i),y(i));
    k2(i) = dydx(x(i)+(1/2).*h, y(i)+k1(i).*(1/2).*h);
    k3(i) = dydx(x(i)+(1/2).*h, y(i)+k2(i).*(1/2).*h);
    k4(i) = dydx(x(i)+h, y(i)+k3(i).*h);
    
    y(i+1)=y(i)+(1/6).*(k1(i)+2.*k2(i)+2.*k3(i)+k4(i)).*h;
end

uRK=[x;y];
plot(x, uRK(2,:),'k');
legend('Exact', 'Euler', 'Huen', 'Ralston', 'RK 4th Order', 'location', 'northwest');
xlabel('X'), ylabel('Y');
title('22.2');
hold off;



%%%%% 22.7 %%%%%
t = [0,0.1,0.2];
y = [2,1.9819,1.9342]; %%simply plotting the data from our results for y and z
z = [4,3.7428,3.5061];
figure(2);
plot(t,y,'bo--');
hold on;
plot(t,z,'ro--');
legend('y', 'z');
xlabel('t'), ylabel('Value');
title('22.7');
hold off;



%%%%% 22.15 %%%%%
k = 20;
m = 20;
c = 5; %%underdamped
dx1 = @(x1,x2) x2; %%derivative of position is velocity
dx2 = @(x1,x2) -(k./m).*x1 - (c./m).*x2; %%derivative of velocity equation
x(1,1) = 1; %%x(1,...) represents position
x(2,1) = 0; %%x(2,...) represents velocity
h = 0.05;
t = [0:0.05:15];


%%since each k1, k2, ...kn is a vector and has two components, let
%%kn(1,...) be the first component, and kn(2,...) be the second component
for i=1:length(t)-1;
    k1(1,i) = dx1(x(1,i), x(2,i));
    k1(2,i) = dx2(x(1,i), x(2,i));
    k2(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(1,i));
    k2(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(2,i));
    k3(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(1,i));
    k3(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(2,i));
    k4(1,i) = dx1(x(1,i)+h, x(2,i)+h.*k3(1,i));
    k4(2,i) = dx2(x(1,i)+h, x(2,i)+h.*k3(2,i));
    x(1,i+1) = x(1,i) + (1/6).*(k1(1,i) + 2.*k2(1,i) + 2.*k3(1,i) + k4(1,i)).*h;
    x(2,i+1) = x(2,i) + (1/6).*(k1(2,i) + 2.*k2(2,i) + 2.*k3(2,i) + k4(2,i)).*h;
end

uc5pos=[t;x(1,:)];
%%uc5vel=[t;x(2,:)];
figure(3);
plot(t,uc5pos(2,:),'r');
hold on;
%%plot(t,uc5vel(2,:),'b'); %%not necessary to plot velocities, but it's cool to look at them

c = 40; %%critically damped
dx1 = @(x1,x2) x2; 
dx2 = @(x1,x2) -(k./m).*x1 - (c./m).*x2;
for i=1:length(t)-1;
    k1(1,i) = dx1(x(1,i), x(2,i));
    k1(2,i) = dx2(x(1,i), x(2,i));
    k2(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(1,i));
    k2(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(2,i));
    k3(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(1,i));
    k3(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(2,i));
    k4(1,i) = dx1(x(1,i)+h, x(2,i)+h.*k3(1,i));
    k4(2,i) = dx2(x(1,i)+h, x(2,i)+h.*k3(2,i));
    x(1,i+1) = x(1,i) + (1/6).*(k1(1,i) + 2.*k2(1,i) + 2.*k3(1,i) + k4(1,i)).*h;
    x(2,i+1) = x(2,i) + (1/6).*(k1(2,i) + 2.*k2(2,i) + 2.*k3(2,i) + k4(2,i)).*h;
end

uc40pos=[t;x(1,:)];
%%uc40vel=[t;x(2,:)];
plot(t,uc40pos(2,:),'b');
%%plot(t,uc40vel(2,:),'b--');

c = 200; %%overdamped
dx1 = @(x1,x2) x2; 
dx2 = @(x1,x2) -(k./m).*x1 - (c./m).*x2;
for i=1:length(t)-1;
    k1(1,i) = dx1(x(1,i), x(2,i));
    k1(2,i) = dx2(x(1,i), x(2,i));
    k2(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(1,i));
    k2(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k1(2,i));
    k3(1,i) = dx1(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(1,i));
    k3(2,i) = dx2(x(1,i)+(1/2).*h, x(2,i)+(1/2).*h.*k2(2,i));
    k4(1,i) = dx1(x(1,i)+h, x(2,i)+h.*k3(1,i));
    k4(2,i) = dx2(x(1,i)+h, x(2,i)+h.*k3(2,i));
    x(1,i+1) = x(1,i) + (1/6).*(k1(1,i) + 2.*k2(1,i) + 2.*k3(1,i) + k4(1,i)).*h;
    x(2,i+1) = x(2,i) + (1/6).*(k1(2,i) + 2.*k2(2,i) + 2.*k3(2,i) + k4(2,i)).*h;
end

uc200pos=[t;x(1,:)];
%%uc200vel=[t;x(2,:)];
plot(t,uc200pos(2,:),'g');
%%plot(t,uc200vel(2,:),'b-.');

legend('Underdamped (c=5)', 'Critically Damped (c=40)', 'Overdamped (c=200)', 'Location','Northeast');
title('22.15');
xlabel('Time (s)'), ylabel('Disp (m)');
hold off;


