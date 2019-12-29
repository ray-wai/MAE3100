%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 7, 10/11/19

%%Using multiple linear regression to fit a plane to 3D data
%%Given data
x1 = [0 1 1 2 2 3 3 4 4]';
x2 = [0 1 2 1 2 1 2 1 2]';
y = [15.1 17.9 12.7 25.6 20.5 35.1 29.7 45.4 40.2]';

%%Using necessary equations to solve
A = [length(x1), sum(x1), sum(x2); sum(x1), sum(x1.^2), sum(x1.*x2); sum(x2), sum(x1.*x2), sum(x2.^2)];
B = [sum(y); sum(y.*x1); sum(y.*x2)];
C = A\B;
fprintf("\nThe result from necessary equations is: y = %.4f + %.4fx1 + %.4fx2", C(1), C(2), C(3));

D = (C(1)) + (C(2).*5) + (C(3).*3); %%finding the value of y when x1 = 5, x2 = 3
fprintf("\nThe value of y when x1 = 5 and x2 = 3 is: %.4f\n", D);

%%graphing the results (Using coefficients found in necessary equations)
%%This is not required, but I did this as part of the prelab and it is a
    %%good visual representation of the solution to both parts 1 and 2
scatter3(x1, x2, y, 'ro');
hold on;
plot3(5, 3, D, 'bo'); %%graphing the point x1 = 5, x2 = 3, y = D
[a b] = meshgrid(0:0.1:5); %%defining a range of values to be plotted in the x1 x2 plane
c = ((C(1)) + (C(2).*a) + (C(3).*b)); %%from necessary equations
surf(a, b, c); %%creating the surface
colorbar;
xlabel("X1");
ylabel("X2");
zlabel("Y");
title("Multiple Linear Regression");