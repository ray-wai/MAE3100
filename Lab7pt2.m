%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 7 part 2, 10/12/19

%%creating a polyfit (general LSR) function and fitting a cubic polynomial
    %%with given data

format short g;

%%fitting the given data to a cubic model
x1 = [3; 4; 5; 7; 8; 9; 11; 12];
y1 = [1.6; 3.4; 4.4; 3.5; 2.2; 2.7; 3.8; 4.6];
m1 = 3;

a1 = super_cool_glsr(x1, y1, m1);

scatter(x1, y1);
hold on;
x2 = [2:0.1:13];
fn = a1(1) + a1(2).*x2 + a1(3).*(x2.^2) + a1(4).*(x2.^3);
fprintf("\nThe resulting cubic fit from my function is; y = %.4f + %.4fx + %.4fx^2 + %.4fx^3", a1(1), a1(2), a1(3), a1(4));

a2 = polyfit(x1, y1, m1);
fprintf("\nThe resulting cubic fit from polyfit is; y = %.4f + %.4fx + %.4fx^2 + %.4fx^3\n\n", a2(4), a2(3), a2(2), a2(1));

plot(x2, fn);
legend('Data Points', 'Fit Curve', 'Location', 'Southeast');
title("Cubic fit of Given Data using a General Least Squares Regression Function");
xlabel("X");
ylabel("Y");

%%super_cool_glsr:
    %%This function is passed x and y data and a polynomial size
    %%The data is placed into the general Z matrix and necessary operations
        %%are performed
    %%The matrix returned contains coefficients in INCREASING powers of X
function [a] = super_cool_glsr(X, Y, M)
    N=length(X);
    Z=zeros(N,M+1);
    Z(:, 1) = ones; %%X.^0 = 1
    
    for i=2:M+1
        Z(:,i) = X.^(i-1); %%Z's columns are in increasing powers of X
    end
    a = inv(Z'*Z)*Z'*Y;
end