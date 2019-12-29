%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 8, 10/18/19

%%Using Newtons-Interpolating Polynomial to curve fit given data and
    %%graphing the results.
    
%%PART B (3rd Order)
x2 = [3 4 2.5 5];
y2 = [7 3 6.5 1];

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
poi = 3.4;
yint2 = b2(1,1); 
for j = 1:n-1
    xt = xt*(poi-x2(j)); 
    yint2 = yint2 + b2(1,j+1) * xt;
end
fprintf("The funtion value using Newtons 3rd Order Polynomial at X = %.4f is %.4f\n", poi, yint2);

%%graphing the result
plot(x2,y2, 'ro'); %%plotting data points
hold on;
range = 2:0.1:6;
f2 = @(a) b2(1,1) + b2(1,2).*(a-x2(1)) + b2(1,3).*(a-x2(1)).*(a-x2(2)) + b2(1,4).*(a-x2(1)).*(a-x2(2)).*(a-x2(3));
plot(range, f2(range),'g'); %%plotting Newtons Polynomial 3rd order
plot(poi, f2(poi), 'go'); %%plotting newtons value @ POI 3rd order
xlabel("X");
ylabel("f(X)");
title("Using Netwons-Interpolating Polynomial to Estimate a Function");
legend('Given Data', 'Newtons 3rd Order Polynomial', 'Newtons 3rd Order @X = 3.4');
