%%Raymond Waidmann
%%MAE 3100 Lab 1 8/30/19
%%Using Eulers method and different step sizes to approximate an ODE then
    %%plotting the approximations and errors

g = 9.81;
cd = 0.15;
m = 68.1;
stepsize = [1, 2, 4, 8]; %%an array of Step Sizes
steps = (24./stepsize) + 1; %%number of Steps for a given step size

v = zeros(4, 25);

for i = 1:4; %%looping throgh the rows
    v(i, 1) = 0;
    for j = 1:((steps(i))-1); %%looping through the colums, ending at appropriate steps
        v(i, j+1) = v(i,j) + ((g-((cd./m).*(v(i,j)).^2)).*stepsize(i)); 
    end; %%setting the values of 4x25 matrix v to the eulers approx for all 4 step sizes
end;

a = 0:24; %%range of time values for step size 1
disp('Time (s)     Velocity (m/s)'); %%displaying step size 1 results
disp([a', v(1,:)']);

b = 0:2:24; %%range of time values for step size 2
disp('Time (s)     Velocity (m/s)'); %%displaying step size 2 results
disp([b', v(2,1:13)']);

c = 0:4:24; %%range of time values for step size 4
disp('Time (s)     Velocity (m/s)'); %%displaying step size 4 results
disp([c', v(3,1:7)']);

d = 0:8:24; %%range of time values for step size 8
disp('Time (s)     Velocity (m/s)'); %%displaying step size 8 results
disp([d', v(4,1:4)']);

truevalue = zeros(1, 25);
for k = 2:25; %%using a for loop to initialize the true values from the given equation
    truevalue(k) = (sqrt((g.*m)./cd).*tanh(sqrt((g.*cd)/m).*(k-1)));
end;

%%graphing true value, and all 4 euler approximations on the same graph
figure(1); 
plot(a, truevalue(1,:), 'ko--');
hold on;
plot(a, v(1,:), 'ro-');
plot(b, v(2,1:13), 'yo-');
plot(c, v(3,1:7), 'go-');
plot(d, v(4,1:4), 'bo-');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('True Velocity', 'Step Size = 1', 'Step Size = 2', 'Step Size = 4', 'Step Size = 8', 'Location', 'southeast');
title('Velocity vs Time');
hold off;

trueerror = zeros(1,4); %%trueerror = true value - approx
reltrueerror = zeros(1,4); %% RTE = |trueerror/truevalue| *100%
relapproxerror = zeros(1,4); %%RAE = |(present approx - previous approx)/present approx| *100%

%%initializing all 3 error types in one for loop (error @24 seconds)
for m = 1:4;
    trueerror(m) = truevalue(25) - v(m,(steps(m)));
    reltrueerror(m) = abs((trueerror(m)/truevalue(25)))*100;
    relapproxerror(m) = abs(((v(m,steps(m)))-(v(m,(steps(m)-1))))/v(m,(steps(m))))*100;
end;

%%graphing relative true error and relative approx error on the same graph
figure(2);
plot(stepsize, reltrueerror(1,:), 'ko-');
hold on;
plot(stepsize, relapproxerror(1,:), 'bo-');
xlabel('Step size (s)');
ylabel('Percent Error (%)');
title('Relative True Error and Relative Approximation Error');
legend('Relative True Error', 'Relative Approximation Error', 'Location', 'northwest');
hold off;
