%%Appendix A

%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Final Project Main Program
clear all;
clc;
format long;
warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale'); %%suppressing a warning message from our implementation of polyfit
load('pdata.mat'); %%loading in the distance and time data sets from the given file

%%*****PART 1.1*****%%
%%SEE REPORT FOR EXPLANATION OF THE DETERMINATION OF POLYNOMIAL DEGREES TESTED
%%3rd order regression
t3=[length(t) sum(t) sum(t.^2) sum(t.^3);
    sum(t) sum(t.^2) sum(t.^3) sum(t.^4); 
    sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5); %%creating the matrices of the normal equations
    sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6)];
r3=[sum(r); sum(t.*r) ; sum(t.^2.*r); sum(t.^3.*r)];

%%calling a gauss elimination function I wrote to solve for the coefficients of the polynomial
%%GaussEliminationFnct returns coefficients in increasing powers of x
coeff3 = GaussEliminationFnct_FS19(t3,r3); 
coeff3func = @(x) coeff3(4).*(x.^3)+coeff3(3).*(x.^2)+coeff3(2).*(x)+coeff3(1); 

%%4th order regression - see 3rd order comments
t4=[length(t) sum(t) sum(t.^2) sum(t.^3) sum(t.^4);
    sum(t) sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5); 
    sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6);
    sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6) sum(t.^7);
    sum(t.^4) sum(t.^5) sum(t.^6) sum(t.^7) sum(t.^8)];
r4=[sum(r); sum(t.*r) ; sum(t.^2.*r); sum(t.^3.*r); sum(t.^4.*r)];

coeff4 = GaussEliminationFnct_FS19(t4,r4);
coeff4func = @(x) coeff4(5).*(x.^4)+coeff4(4).*(x.^3)+coeff4(3).*(x.^2)+coeff4(2).*(x)+coeff4(1);

%%5th order regression - see 3rd order comments
t5=[length(t) sum(t) sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5);
    sum(t) sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6); 
    sum(t.^2) sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6) sum(t.^7);
    sum(t.^3) sum(t.^4) sum(t.^5) sum(t.^6) sum(t.^7) sum(t.^8);
    sum(t.^4) sum(t.^5) sum(t.^6) sum(t.^7) sum(t.^8) sum(t.^9);
    sum(t.^5) sum(t.^6) sum(t.^7) sum(t.^8) sum(t.^9) sum(t.^10);];
r5=[sum(r); sum(t.*r) ; sum(t.^2.*r); sum(t.^3.*r); sum(t.^4.*r); sum(t.^5.*r)];

coeff5 = GaussEliminationFnct_FS19(t5,r5);
coeff5func = @(x) coeff5(6).*(x.^5)+coeff5(5).*(x.^4)+coeff5(4).*(x.^3)+coeff5(3).*(x.^2)+coeff5(2).*(x)+coeff5(1);

%%*****PART 1.3*****%%
%%SEE REPORT FOR ANALYSIS
%%calling polyfit to solve for the coefficients of the polynomial for each order
%%polyfit returns coefficients in decreasing powers of x contrary to the increasing powers of x returned by the GaussEliminationFnct
pf3 = polyfit(t,r,3); 
pf3func = @(x) pf3(1).*(x.^3)+pf3(2).*(x.^2)+pf3(3).*(x)+pf3(4); 
pf4 = polyfit(t,r,4); 
pf4func = @(x) pf4(1).*(x.^4)+pf4(2).*(x.^3)+pf4(3).*(x.^2)+pf4(4).*(x)+pf4(5); 
pf5 = polyfit(t,r,5); 
pf5func = @(x) pf5(1).*(x.^5)+pf5(2).*(x.^4)+pf5(3).*(x.^3)+pf5(4).*(x.^2)+pf5(5).*(x)+pf5(6);
%%******************%%

%%*****PART 1.1 CONTINUED*****%%
%%calculation of errors to determine which polynomial of the three is the best fit
%%in this for loop, we are subtracting the distance values from the best fit polynomial
%%from the actual data points at that particular time. As per least squares regression, 
%%we sum the squares of each of these errors for each order of polynomial we have fitted
error3 = 0; error4 = 0; error5 = 0;
for a = 1 : length(t) 
    error3 = error3 + (abs(pf3func(t(a)) - r(a))).^2; 
    error4 = error4 + (abs(pf4func(t(a)) - r(a))).^2; 
    error5 = error5 + (abs(pf5func(t(a)) - r(a))).^2; 
end

fprintf("\nComparison of Squared Residuals for each Order of Polynomial in km^2 (1.1)\n3rd Order      4th Order      5th Order\n");
fprintf("%d   %d   %d\n\n", error3, error4, error5); %%displaying the results in the command window


%%*****PART 1.2*****%%
%%SEE REPORT


%%*****PART 1.3 CONTINUED*****%%
%%displaying functions from GaussEliminationFnct and polyfit to the command window to demonstrate they are equivalent
fprintf("Results of 5th Order Equations (1.3)\n");
fprintf("Gauss Elimination: %dfx^5 + %dfx^4 + %dfx^3 + %dfx^2 + %dfx + %df\n", coeff5(6), coeff5(5), coeff5(4), coeff5(3), coeff5(2), coeff5(1));
fprintf("Polyfit          : %dfx^5 + %dfx^4 + %dfx^3 + %dfx^2 + %dfx + %df\n\n", pf5(1), pf5(2), pf5(3), pf5(4), pf5(5), pf5(6));


%%*****PART 1.4*****%%
%%plotting distance (r) vs time (t) from the given pdata file
%%plotting the polynomials from both the GaussEliminationFnct and polyfit on the graph for all 3 polynomials
figure(1);
plot(t,r,'g'), grid on, hold on; 
xlabel('Time'), ylabel ('Distance'), title('Distance vs Time Part 1'); 
plot(t,coeff3func(t),'r'); 
plot(t,pf3func(t),'r--'); 
plot(t,coeff4func(t),'b'); 
plot(t,pf4func(t),'b--');
plot(t,coeff5func(t),'k');
plot(t,pf5func(t),'k--');


%%*****PART 1.5*****%%
%%Calling a 2 Golden Section Search functions I wrote to find the maximum and minimum timae values of the polynomial with the best fit. 
%%Displaying the man and min time and distance to the command window for user reference. 
tmax = GoldSecSearchMax_FS19(pf5func, 10000, 20000); %%The maximum is visibly inbetween the time interval [10000,20000]
tmin = GoldSecSearchMin_FS19(pf5func, 25000, 30000); %%The minimum is visibly inbetween the time interval [25000,30000]
plot(tmax,pf5func(tmax),'ro'); 
plot(tmin,pf5func(tmin),'bo'); 
legend('Data', '3rd Gauss', '3rd Polyfit', '4th Gauss', '4th Polyfit', '5th Gauss', '5th Polyfit', '5th Max', '5th Min');
hold off;

fprintf("Maximum and Minimum points using Golden Section Search (1.5)\n");
fprintf("Maximum Time: %d\nMaximum Distance: %d\n", tmax,pf5func(tmax)); 
fprintf("Minimum Time: %d\nMinimum Distance: %d\n", tmin,pf5func(tmin));


%%*****PART 2.1*****%%
%%SEE RK_FS19


%%*****PART 2.2*****%%
%%Computing the trajectory of the satellite at three times steps each with 3 different order of runge kutta. 
%%Throughout the rest of the program, the any vaiable with 2 numerical subscripts follows the following format: 
    %%For example: youtXY
        %%X = step size used (1, 2, or 3, where 1 is the largest stepsize)
        %%Y = order of RK used (1, 2, or 4)
tspan = [0, 32000];
y0 = [5000; 10000; 2100; -5.9925; 1.9254; 3.2456];
h1 =(tspan(2)-tspan(1))./50; %%given constants
h2 =(tspan(2)-tspan(1))./100;
h3 =(tspan(2)-tspan(1))./200;

rk = 1;
[~, yout11] = RK_FS19(@Ode_Function_FS19, tspan, y0, h1, rk); 
[~, yout21] = RK_FS19(@Ode_Function_FS19, tspan, y0, h2, rk); %%'~' since tout1, 2, and 3 are the same regardless of the order of rk used
[~, yout31] = RK_FS19(@Ode_Function_FS19, tspan, y0, h3, rk);    

rk = 2;
[~, yout12] = RK_FS19(@Ode_Function_FS19, tspan, y0, h1, rk); 
[~, yout22] = RK_FS19(@Ode_Function_FS19, tspan, y0, h2, rk); 
[~, yout32] = RK_FS19(@Ode_Function_FS19, tspan, y0, h3, rk); 

rk = 4;
[tout1, yout14] = RK_FS19(@Ode_Function_FS19, tspan, y0, h1, rk); 
[tout2, yout24] = RK_FS19(@Ode_Function_FS19, tspan, y0, h2, rk); %%tout is independent of the order of RK used, thus the subscript is which h is used
[tout3, yout34] = RK_FS19(@Ode_Function_FS19, tspan, y0, h3, rk); 


%%*****PART 2.3*****%%
%%plotting the result (trajectories) of RK for each step size, each on a seperate graph
figure(2), grid on, hold on;
plot3(yout11(:,1), yout11(:,2), yout11(:,3), 'k');
plot3(yout12(:,1), yout12(:,2), yout12(:,3), 'b');
plot3(yout14(:,1), yout14(:,2), yout14(:,3), 'r');
DrawEarth(6378.137);
view(45,45);
grid on;
legend('1st Order RK', '2nd Order RK', '4th Order RK');
xlabel('X Distance (km)'), ylabel('Y Distance (km)'), zlabel('Z Distance (km)');
title('Trajectories of the Satellite for Step Size (tf-to)/50 (Part 2.3)'), hold off;

figure(3), grid on, hold on;
plot3(yout21(:,1), yout21(:,2), yout21(:,3), 'k');
plot3(yout22(:,1), yout22(:,2), yout22(:,3), 'b');
plot3(yout24(:,1), yout24(:,2), yout24(:,3), 'r'); 
DrawEarth(6378.137);
view(45,45);
legend('1st Order RK', '2nd Order RK', '4th Order RK');
xlabel('X Distance (km)'), ylabel('Y Distance (km)'), zlabel('Z Distance (km)');
title('Trajectories of the Satellite for Step Size (tf-to)/100 (Part 2.3)'), hold off;

figure(4), grid on, hold on;
plot3(yout31(:,1), yout31(:,2), yout31(:,3), 'k');
plot3(yout32(:,1), yout32(:,2), yout32(:,3), 'b');
plot3(yout34(:,1), yout34(:,2), yout34(:,3), 'r');
DrawEarth(6378.137);
view(45,45);
legend('1st Order RK', '2nd Order RK', '4th Order RK');
xlabel('X Distance (km)'), ylabel('Y Distance (km)'), zlabel('Z Distance (km)');
title('Trajectories of the Satellite for Step Size (tf-to)/200 (Part 2.3)'), hold off;


%%*****PART 2.4*****%%
%%calculating distances the satellite is from earth for each step size and RK method used
dist11 = sqrt(abs(yout11(:,1).^2) + abs(yout11(:,2).^2) + abs(yout11(:,1).^2));
dist12 = sqrt(abs(yout12(:,1).^2) + abs(yout12(:,2).^2) + abs(yout12(:,1).^2));
dist14 = sqrt(abs(yout14(:,1).^2) + abs(yout14(:,2).^2) + abs(yout14(:,1).^2));

dist21 = sqrt(abs(yout21(:,1).^2) + abs(yout21(:,2).^2) + abs(yout21(:,1).^2));
dist22 = sqrt(abs(yout22(:,1).^2) + abs(yout22(:,2).^2) + abs(yout22(:,1).^2)); 
dist24 = sqrt(abs(yout24(:,1).^2) + abs(yout24(:,2).^2) + abs(yout24(:,1).^2));

dist31 = sqrt(abs(yout31(:,1).^2) + abs(yout31(:,2).^2) + abs(yout31(:,1).^2));
dist32 = sqrt(abs(yout32(:,1).^2) + abs(yout32(:,2).^2) + abs(yout32(:,1).^2));
dist34 = sqrt(abs(yout34(:,1).^2) + abs(yout34(:,2).^2) + abs(yout34(:,1).^2));

%%using the distances we calculated to plot each distance and RK method used on 3 graphs, 1 for each stepsize
figure(5), grid on, hold on;
plot(tout1, dist11, 'k');
plot(tout1, dist12, 'b');
plot(tout1, dist14, 'r');
legend('1st Order RK', '2nd Order RK', '4th Order RK', 'location', 'northwest');
xlabel('Time (Seconds)'), ylabel('Distance of Satellite from Earth (km)');
title('Distance of the Satellite from Earth for Step Size (tf-to)/50 (Part 2.4)'), hold off;

figure(6), grid on, hold on;
plot(tout2, dist21, 'k');
plot(tout2, dist22, 'b');
plot(tout2, dist24, 'r');
legend('1st Order RK', '2nd Order RK', '4th Order RK', 'location', 'northwest'); 
xlabel('Time (Seconds)'), ylabel('Distance of Satellite from Earth (km)');
title('Distance of the Satellite from Earth for Step Size (tf-to)/100 (Part 2.4)'), hold off;

figure(7), grid on, hold on;
plot(tout3, dist31, 'k');
plot(tout3, dist32, 'b');
plot(tout3, dist34, 'r');
legend('1st Order RK', '2nd Order RK', '4th Order RK', 'location', 'northwest');
xlabel('Time (Seconds)'), ylabel('Distance of Satellite from Earth (km)');
title('Distance of the Satellite from Earth for Step Size (tf-to)/200 (Part 2.4)'), hold off;


%%*****Part 2.5*****%%
%%Computing the trajectory of the satellite using ode45
[~, yode1] = ode45(@Ode_Function_FS19, tout1, y0); %%tout1, 2, and 3 are used because they contain the correct time steps  
[~, yode2] = ode45(@Ode_Function_FS19, tout2, y0); %%vectors for each step size; we can just reuse them from before
[~, yode3] = ode45(@Ode_Function_FS19, tout3, y0); 

%%calculating the satellite distances for each of the ode45 solutions
distode1 = sqrt(abs(yode1(:,1).^2) + abs(yode1(:,2).^2) + abs(yode1(:,1).^2));
distode2 = sqrt(abs(yode2(:,1).^2) + abs(yode2(:,2).^2) + abs(yode2(:,1).^2)); 
distode3 = sqrt(abs(yode3(:,1).^2) + abs(yode3(:,2).^2) + abs(yode3(:,1).^2));

%%relative approx error for each stepsize and all orders of RK
rae11 = abs((distode1 - dist11)./(distode1)).*100;
rae12 = abs((distode1 - dist12)./(distode1)).*100; 
rae14 = abs((distode1 - dist14)./(distode1)).*100;
rae21 = abs((distode2 - dist21)./(distode2)).*100;
rae22 = abs((distode2 - dist22)./(distode2)).*100; 
rae24 = abs((distode2 - dist24)./(distode2)).*100;
rae31 = abs((distode3 - dist31)./(distode3)).*100;
rae32 = abs((distode3 - dist32)./(distode3)).*100; 
rae34 = abs((distode3 - dist34)./(distode3)).*100;

%%plotting the errors between ode45 and each order of rk for all three stepsizes
figure(8), hold on;
subplot(3,1,1), plot(tout1, rae11), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/50, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout1, rae12), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/50, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout1, rae14), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/50, RK Order 4 (2.5)');
hold off;

figure(9), hold on;
subplot(3,1,1), plot(tout2, rae21), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/100, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout2, rae22), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/100, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout2, rae24), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/100, RK Order 4 (2.5)');
hold off;

figure(10), hold on;
subplot(3,1,1), plot(tout3, rae31), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/200, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout3, rae32), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/200, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout3, rae34), grid on;
xlabel('Time (s)'), ylabel('Rel Approx Err (%)'), title('Step Size h = (tf-to)/200, RK Order 4 (2.5)');
hold off;

%%true error for each stepsize and all orders of RK
te11 = abs((distode1 - dist11));
te12 = abs((distode1 - dist12));
te14 = abs((distode1 - dist14));
te21 = abs((distode2 - dist21));
te22 = abs((distode2 - dist22));
te24 = abs((distode2 - dist24));
te31 = abs((distode3 - dist31));
te32 = abs((distode3 - dist32));
te34 = abs((distode3 - dist34));

%%plotting the errors between ode45 and each order of rk for all three stepsizes
figure(11), hold on;
subplot(3,1,1), plot(tout1, te11), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/50, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout1, te12), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/50, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout1, te14), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/50, RK Order 4 (2.5)');
hold off;

figure(12), hold on;
subplot(3,1,1), plot(tout2, te21), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/100, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout2, te22), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/100, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout2, te24), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/100, RK Order 4 (2.5)');
hold off;

figure(13), hold on;
subplot(3,1,1), plot(tout3, te31), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/200, RK Order 1 (2.5)');
subplot(3,1,2), plot(tout3, te32), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/200, RK Order 2 (2.5)');
subplot(3,1,3), plot(tout3, te34), grid on;
xlabel('Time (s)'), ylabel('Absolute Error (km)'), title('Step Size h = (tf-to)/200, RK Order 4 (2.5)');
hold off;