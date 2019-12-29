%%Raymond Waidmann
%%MAE 3100, Lab 5 Assignment 2, 9/30/19
%%Using Gauss-Seidel method to solve a 3x3 system of equations

a = [10 2 -1; -3 -6 2; 1 2 5];
b = [27; -61.3; -21.2];

x1 = zeros(1,10);
x2 = zeros(1,10);
x3 = zeros(1,10);

RAE = 100; %%arbitrary number > 5 so the loop executes at least once
i = 1; %%a counter for number of iterations

while RAE > 5
    x1(i+1) = (b(1,1) - (a(1,2).*x2(i)) - (a(1,3)*x3(i)))/(a(1,1));
    x2(i+1) = (b(2,1) - (a(2,1).*x1(i+1)) - (a(2,3)*x3(i)))/(a(2,2));
    x3(i+1) = (b(3,1) - (a(3,1).*x1(i+1)) - (a(3,2)*x2(i+1)))/(a(3,3));
    
    RAE1=abs((x1(i+1)-x1(i))/x1(i+1)).*100;
    RAE2=abs((x2(i+1)-x2(i))/x2(i+1)).*100;
    RAE3=abs((x3(i+1)-x3(i))/x3(i+1)).*100;
    
    c = [RAE1 RAE2 RAE3];
    RAE = max(c);
    i = i + 1;
end

fprintf("\nThe value for x1 is %f and has a RAE of %f percent\n", x1(i), RAE1);
fprintf("The value for x2 is %f and has a RAE of %f percent\n", x2(i), RAE2);
fprintf("The value for x3 is %f and has a RAE of %f percent\n", x3(i), RAE3);
fprintf("These results were found after %d iterations\n\n", (i-1));