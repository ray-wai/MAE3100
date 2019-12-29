%%Raymond Waidmann
%%MAE 3100 Lab 4 9/20/19
%%Using MatLab to perform a Golden Section Search

f = @(x) (4.*x) - (1.8.*(x.^2)) + (1.2.*(x.^3)) - (0.3.*(x.^4));
phi = ((1+sqrt(5))./2); %%golden ratio
es = 0.1; %%stop condition for RelAppErr

xl = zeros(1,20); %%lower bound
xl(1) = -2;
xu = zeros(1,20); %%upper bound
xu(1) = 4;
xA = zeros(1,20); %%intermediate point A
xB = zeros(1,20); %%intermediate point B
xOpt = zeros(1,20); %%array of the optimum for each iteration
fxOpt = zeros(1,20);
d = zeros(1,20);  %%distance between xl->xA and x2->xB
ea = zeros(1,20); %%relative approx error
ea(1) = 100; %%arbitrary number greater than es
i = 1; %%a counter for the iterations

while(ea(i) > es)
    d(i) = (phi-1).*(xu(i) - xl(i));
    xA(i) = xl(i) + d(i); %%calculating the two intermediate points
    xB(i) = xu(i) - d(i);
    
    if (-f(xA(i)) < -f(xB(i))) %%minus needed since we are finding Max
        xl(i+1) = xB(i); %%new lower bound
        xu(i+1) = xu(i); %%upper bound stays the same
        xOpt(i) = xA(i); %%optimum guess
    else
        xl(i+1) = xl(i); %%lower bound stays the same
        xu(i+1) = xA(i); %%new upper bound
        xOpt(i) = xB(i); %%optimum guess
    end
    ea(i+1) = (2-phi).*abs((xu(i+1)-xl(i+1))/(xOpt(i))).*100; %%RelAppErr
    fxOpt(i) = (f(xOpt(i))); %%calculating the function value @xOpt
    i = i+1; %%incrementing the counter
end
iterations = 1:(i-1);

%%displaying the results to the screen
fprintf("\nIteration No.         Xopt      f(Xopt)           Ea\n");
disp([iterations(1,:)', xOpt(1:(i-1))', fxOpt(1:(i-1))', ea(2:i)']);

%%showing the solution graphically and circling the optimum
x = [-2:0.01:4]; 
plot(x, f(x));
hold on;
plot(xOpt(i-1), fxOpt(i-1), 'ro');
title('Finding the Optimum of f(x) = 4x - 1.8x^2 + 1.2x^3 - 0.3x^4');
ylabel('f(x)');
xlabel('x');
hold off;