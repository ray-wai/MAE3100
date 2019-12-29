%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 8, 10/18/19

%%Using Newtons-Interpolating Polynomial to curve fit given data and
    %%graphing the results.
    
%%given data
x = [1920 1930 1940 1950 1960 1970 1980 1990 2000];
y = [106.45 123.07 132.10 151.27 180.67 204.05 227.23 246.46 281.42];

%%PART A (2nd Order)
x1 = [1970 1980 1990];
y1 = [204.05 227.23 246.46];

%%creating the divided difference table for quadratic polynomial
n = length(x1);
if length(y1) ~= n, error('X and Y must be same length'); end
b1 = zeros(n,n);
b1(:,1) = y1(:); %%first column is given y values
for j = 2:n
    for i = 1:n-j+1
        b1(i,j) = (b1(i+1, j-1)-b1(i, j-1))/(x1(i+j-1)-x1(i));
    end
end

%%interpolation at a point
xt = 1; %%xt is (x-x1)*(x-x2)*...*(x-xn-1)
poi = 2000; %%point of interest
yint = b1(1,1); %%yint is running value of coeff*xt
for j = 1:n-1
    xt = xt*(poi-x1(j)); %%number is the x-value of interest (2000)
    yint = yint + b1(1,j+1) * xt;
end
fprintf("\nThe funtion value using Newtons 2nd Order Polynomial at X = %.4f is %.4f\n", poi, yint);

%%PART B (3rd Order)
x2 = [1960 1970 1980 1990];
y2 = [180.67 204.05 227.23 246.46];

%%creating the divided difference table for cubic polynomial
%%See comments in part A
n = length(x2);
if length(y2) ~= n, error('X and Y must be same length'); end
b2 = zeros(n,n);
b2(:,1) = y2(:); 
for j = 2:n
    for i = 1:n-j+1
        b2(i,j) = (b2(i+1, j-1)-b2(i, j-1))/(x2(i+j-1)-x2(i));
    end
end

%%interpolation at a point
%%See comments in part A
xt = 1; 
%%POI is still 2000 from previous part
yint2 = b2(1,1); 
for j = 1:n-1
    xt = xt*(poi-x2(j)); 
    yint2 = yint2 + b2(1,j+1) * xt;
end
fprintf("The funtion value using Newtons 3rd Order Polynomial at X = %.4f is %.4f\n", poi, yint2);

%%graphing the result
plot(x,y, 'ro'); %%plotting data points
hold on;
range = 1910:1:2010;
f = @(a) b1(1,1) + b1(1,2).*(a-x1(1)) + b1(1,3).*(a-x1(1)).*(a-x1(2));
f2 = @(a) b2(1,1) + b2(1,2).*(a-x2(1)) + b2(1,3).*(a-x2(1)).*(a-x2(2)) + b2(1,4).*(a-x2(1)).*(a-x2(2)).*(a-x2(3));
plot(range, f(range),'b'); %%plotting Newtons Polynomial 2nd order
plot(poi, f(poi), 'bo'); %%plotting newtons value @ POI 2nd order
plot(range, f2(range),'g'); %%plotting Newtons Polynomial 3rd order
plot(poi, f2(poi), 'go'); %%plotting newtons value @ POI 3rd order
xlabel("X");
ylabel("f(X)");
title("Using Netwons-Interpolating Polynomial to Estimate a Function");
legend('Given Data', 'Newtons 2nd Order Polynomial', 'Newtons 2nd Order @X = 2000','Newtons 3rd Order Polynomial', 'Newtons 3rd Order @X = 2000', 'Location', 'southeast');
fprintf("The true value at X = %.4f is %.4f\n", poi, y(9));
TRE = abs((y(9) - f(poi))./(y(9))).*100; %%true relative error 2nd order
fprintf("True relative error for 2nd Order = %.4f percent\n", TRE);
TRE2 = abs((y(9) - f2(poi))./(y(9))).*100; %%true relative error 3rd order
fprintf("True relative error for 3rd Order = %.4f percent\n\n", TRE2);