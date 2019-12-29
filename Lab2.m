%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 2, 9/6/19

%%Using Newton-Raphsons method and modified secant method to find the largest real root of a function

%%initializing the given function and derivative
f = @(x) x.^3 - 6.*x.^2 + 11.*x - 6.1; 
df = @(x) 3.*x.^2 - 12.*x + 11;
y = 0:0.01:4; %%range of values that includes all three roots

%%plotting the function and using the graphical method
plot(y,f(y));
hold on;
grid on;
a = fzero(f,0);
plot(a,f(a),'ro');
b = fzero(f,2); %%a, b, and c represent the three roots
plot(b,f(b),'ro');
c = fzero(f,4);
plot(c,f(c),'ro');
xlabel('x');
ylabel('f(x)');
title('Using the Graphical Method to Find the Roots');
hold off;

%%using newton raphson method w/ initial quess of 3.5. terminating program when relative approx error < 1%
relapproxerr = zeros(1,5); 
relapproxerr(1) = 10; %%an arbitrary number greater than 1 so the while loop executes
z = zeros(1,5); 
z(1) = 3.5; %%initial guess
i = 1;
while (relapproxerr(i) > 1)
    i = i + 1;
    z(i) = z(i-1) - ((f(z(i-1)))./(df(z(i-1)))); %%finding the next prediction using Newt-Raps
    relapproxerr(i) = abs((z(i) - z(i-1))./z(i)).*100; %%RAE = |(present approx-previous approx)/present approx|*100%
end

%%displaying the results of newt-raps method
disp('Newton - Raphsons Method:'); %%displaying step size 4 results
disp([z(1:i)]);
fprintf("The final guess for the largest real root with initial guess \nx = 3.5 using Newton-Raphsons");
fprintf(" Method is: %f and has a\nrelative approximation error of %f%%.", z(i), relapproxerr(i));
fprintf(" It took %d iterations\nto get to this approximation.\n\n", i);

%%using modified secant method with initial guess 3.5. terminating program when relative approx error < 1%
relapproxerr2 = zeros(1,5);
relapproxerr2(1) = 10;
z2 = zeros(1,5);
z2(1) = 3.5;
j = 1;
delta = 0.01; %%given
while (relapproxerr2(j) > 1)
    j = j + 1;
    z2(j) = z2(j-1) - (((f(z2(j-1))).*(delta.*z2(j-1)))./((f((z2(j-1)).*(1+delta)))-(f(z2(j-1))))); %%modified-sec
    relapproxerr2(j) = abs((z2(j) - z2(j-1))./z2(j)).*100;
end

%%displaying the results of modified secant method
disp('Modified Secant Method:'); %%displaying step size 4 results
disp([z2(1:j)]);
fprintf("The final guess for the largest real root with initial guess \nx = 3.5 using the Modified");
fprintf(" Secant Method is: %f and has a\nrelative approximation error of %f%%. It took ", z2(j), relapproxerr2(j));
fprintf("%d iterations\nto get to this approximation. \n\n(The actual value of the root is %f).\n\n", j, c);
