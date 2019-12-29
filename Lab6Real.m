%%Raymond Waidmann
%%MAE 3100, Lab 6, 10/3/19
%%Using Newton Raphsons method for multiple equations

f1 = @(x,y) (x^2)+(y^2)-5; %%given functions
f2 = @(x,y) y-(x^2)+1; 

f1x = @(x,y) (2*x); %%partial derivatives w respect to x and y
f1y = @(x,y) (2*y);
f2x = @(x,y) (-2*x);
f2y = @(x,y) 1;

%%solving using equations from 12.12
size = 100;
RAEx2 = zeros(1,size); RAEx2(1) = 100; %%Rel approx error arrays
RAEy2 = zeros(1,size); RAEy2(1) = 100;
x2 = zeros(1,size); x2(1) = 1.5; %%arrays that store all the values of x and y
y2 = zeros(1,size); y2(1) = 1.5;

k = 1;

while RAEx2(k)>=0.01 || RAEy2(k)>=0.01
    J = [f1x(x2(k), y2(k)), f1y(x2(k), y2(k)); f2x(x2(k), y2(k)), f2y(x2(k), y2(k))];
    x2(k+1) = x2(k) - (((f1(x2(k),y2(k))).*(f2y(x2(k), y2(k)))) - ((f2(x2(k),y2(k))).*(f1y(x2(k), y2(k)))))./(det(J));
    y2(k+1) = y2(k) - (((f2(x2(k),y2(k))).*(f1x(x2(k), y2(k)))) - ((f1(x2(k),y2(k))).*(f2x(x2(k), y2(k)))))./(det(J));
    RAEx2(k+1) = abs((x2(k+1) - x2(k))./(x2(k+1))).*100;
    RAEy2(k+1) = abs((y2(k+1) - y2(k))./(y2(k+1))).*100;
    k = k+1;
end

%%dislaying the results to the screen
fprintf("\n");
disp('Results of Newt-Raps Method using equations 12.12');
disp(' Iteration         X         Y      Ea(x)     Ea(Y)');
disp([(0:(k-1))', x2(1:k)', y2(1:k)', RAEx2(1:k)', RAEy2(1:k)']);

%%solving using x(i+1) = x(i) - inv(J)*f (12.17)
RAEx = zeros(1,size); RAEx(1) = 100; %%Rel approx error arrays
RAEy = zeros(1,size); RAEy(1) = 100;
x = zeros(1,size); x(1) = 1.5; %%arrays that store all the values of x and y
y = zeros(1,size); y(1) = 1.5;

i = 1;

while RAEx(i)>=0.01 || RAEy(i)>=0.01
    J = [f1x(x(i), y(i)), f1y(x(i), y(i)); f2x(x(i), y(i)), f2y(x(i), y(i))]; %%definition of J matrix
    f = [f1(x(i), y(i)); f2(x(i), y(i))]; %%array of f1 and f2 at the current iteration
    vars = [x(i); y(i)]; %%array of x and y at the current iteration

    tempvars = vars - (inv(J)*f); %%Newt-Raps Formula to find xnew and ynew
    x(i+1) = tempvars(1);
    y(i+1) = tempvars(2);

    RAEx(i+1) = abs((x(i+1) - x(i))./(x(i+1))).*100;
    RAEy(i+1) = abs((y(i+1) - y(i))./(y(i+1))).*100;
    i = i+1;
end

%%dislaying the results to the screen
disp('Results of Newt-Raps Method using equation 12.17');
disp(' Iteration         X         Y      Ea(x)     Ea(Y)');
disp([(0:(i-1))', x(1:i)', y(1:i)', RAEx(1:i)', RAEy(1:i)']);

ezplot(f1);
hold on;
ezplot(f2);
title('Finding the Solution Graphically');
xlabel('X');
ylabel('Y');
plot(x2(k), y2(k), 'ro');
plot(x(i), y(i), 'b*');
legend('Ellipse (f1): x^2 = 5 - y^2', 'Parabola (f2): y = x^2 -1', '12.12 Solution', '12.17 Solution', 'Location', 'southeast');