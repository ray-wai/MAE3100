%%Raymond Waidmann
%%MAE 3100 Lab 3 9/13/19

%%This program uses functions roots and fzero to find the roots of a polynomial
    %%Plot the polynomial and find the roots graphically (uses fzero)
    %%Use Newt-Raps method to find one of the roots with 2 different initial guesses
    %%Plot the root estimate vs iteration for both initial guesses
    
format short g; %%formating statement for good scientific notation output

%%Creating a column vector p that includes the coefficients of the
    %%polynomial; the roots function finds the roots of p.
    %%We use the results of the roots function in order to have good
    %%initial guesses for the individual root assignment in the fzero
    %%statements below. This allows us to create 3 variables, each one
    %%equal to one of the roots. 
p = [0.5, -4, 6, -2];
rts = roots(p)

%%Plotting the polynomial 
figure(1);
x = [-1:0.01:7]; %%range that includes all three roots
f = @(x) ((0.5).*x.^3) - (4.*x.^2) + (6.*x) - 2; %%anonymous functions
df = @(x) ((1.5).*x.^2) - (8.*x) + 6;
plot(x,f(x));
hold on;
grid on;

%%Here we are using initial guesses 0, 2, and 6 based on the return value
    %%of the roots function. By using fzero, we can specifically assign
    %%each root to a different variable and this allows us to plot each
    %%root on the graph with a red circle in the plot statements below.
root1 = fzero(f, 0);
root2 = fzero(f, 2); %%root1, 2, and 3 are the same roots that roots(p) finds
root3 = fzero(f, 6);
fprintf("The roots using the fzero function are %.4f, %.4f, and %.4f\n\n", root1, root2, root3);

plot(root1, f(root1), 'ro');
plot(root2, f(root2), 'ro'); %%circling the roots on the graph
plot(root3, f(root3), 'ro');

title('Using the Graphical Method to Find the Roots');
xlabel('x');
ylabel('f(x)');
hold off;

%%Newton Rapshsons Method with Xo = 4.5
relapproxerr = zeros(1,15); %%an array of the error for each step in Newt-Raps Method
relapproxerr(1) = 100; %%arbitrary number > 0.01
z = zeros(1,15); %%an array of guesses of the reoot for each step in Newt-Raps Method
z(1) = 4.5; %%initial guess
i = 1; 
while (relapproxerr(i) > 0.1) %%this loop executes Newton Rapsons method and calculates Ea for each iteration
    i = i+1;
    z(i) = z(i-1) - ((f(z(i-1)))./(df(z(i-1)))); %%Newt_Raps Method formula
    relapproxerr(i) = abs((z(i) - z(i-1))./(z(i))).*100; %%Rel Approx Error Formula
end    

iterations = [0:(i-1)]; %%a 1D matrix that has all the iteration numbers 
                        %%from 0 to however many times the while loop executed
disp('Newt_Raps Method with Initial Guess x = 4.5');
disp('Iteration No.           xi           ea');
disp([iterations(1,:)', z(1:i)', relapproxerr(1:i)']); %%displaying the results in a table

%%Newton Raps with Xo = 4.43
%%See comments on code for Xo = 4.5, code is identical except different
    %%variable names and initial values
relapproxerr2 = zeros(1,30);
relapproxerr2(1) = 100; 
z2 = zeros(1,30);
z2(1) = 4.43; %%This is the only line that is different, different Xo.
i2 = 1;
while (relapproxerr2(i2) > 0.1)
    i2 = i2+1;
    z2(i2) = z2(i2-1) - ((f(z2(i2-1)))./(df(z2(i2-1))));
    relapproxerr2(i2) = abs((z2(i2) - z2(i2-1))./(z2(i2))).*100;
end    

iterations2 = [0:1:(i2-1)];
disp('Newt_Raps Method with Initial Guess x = 4.43');
disp('Iteration No.           xi           ea');
disp([iterations2(1,:)', z2(1:i2)', relapproxerr2(1:i2)']);

%%Plotting Results for Xo = 4.5
figure(2);
plot(iterations(1,:), z(1:i), 'ro-');
hold on;
xlabel('Iteration No.');
ylabel('Estimate of the Root');
title('Iteration No. vs Estimate of the Root for Xo = 4.5');
grid on;
hold off;

%%Plotting Results for Xo = 4.43
figure(3);
plot(iterations2(1,:), z2(1:i2), 'bo-');
hold on;
xlabel('Iteration No.');
ylabel('Estimate of the Root');
title('Iteration No. vs Estimate of the Root for Xo = 4.43');
grid on;
hold off;

%%After completing this lab, it can be seen that using Newton Raphsons
%%method with an initial Xo = 4.43 does not correctly find the closest
%%root. This is because at X = 4.43, the derivative of the graph is less
%%than zero, but very small in magnitude:

%%df @(4.43) = ((1.5).*(4.43).^2) - (8.*(4.43)) + 6 = -0.00265

%%As we can see from the table and due to the above calculation, the
%%guess of the location of the root after just one iteration is -3939.1.
%%Since we are so far away from any of the three roots - which are all
%%located inbetween x = [0, 7] - there is really no chance that we can
%%find any of the roots in few iterations with a new initial guess so far away. 
%%This proves to be true as even after significantly more iterations than the
%%initial guess at Xo = 4.5, we predict a final value of the root to be
%%0.47457 which is a root, but not the one closest to 4.43. With an initial
%%guess of 4.43, we would expect to find the root at 6.1563; since the sign 
%%of df @(4.43) is negative, this prevents us from finding the root at 
%%6.1563 since we would need a positive slope at an initial guess less 
%%than the value of the root to approach the correct value, which is not
%%the case here.

%%The reason that the initial guess at Xo = 4.5 is able to correctly find a
%%root is because the slope is small in magnitude, but not nearly as small
%%as df @(4.43) and is positive:

%%df @(4.5) = ((1.5).*(4.5).^2) - (8.*(4.5)) + 6 = 0.375

%%As a result of this much larger and positive derivative, the guess 
%%of the root after just one iteration is X = 32.333. Although this is 
%%still guite far from the largest root of 6.1563, it is significantly 
%%close enough for us to approach the value of the root after
%%just 10 total iterations. The sign of df @(4.5) is also positive,
%%allowing us to approach to root we want to find at 6.1563 
%%since 4.5 is less than 6.1563; since the initial guess is less than the
%%value of the root, we need an initially positive slope to approach the
%%root on the first iteration as is the case here. 

%%These results can also be seen on figures 2 and 3; both figures initially
%%spike after one iteration, but the magnitude of the spike is much smaller on
%%figure 2 than it is on figure 3. 
