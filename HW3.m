%%Raymond Waidmann
%%MAE 3100, 18157816, rcw5k2
%%HW3 9/22/2019

%%6.4
f =@(x) (7.*sin(x).*exp(-x))-1;
x = (0:0.1:5);
plot(x,f(x));
hold on;
grid on;
title('f(x)=7*sin(x)*exp(-x)-1');
xlabel('x');
ylabel('f(x)');
root1 = fzero(f, 0);
plot(root1, f(root1), 'ro');
fprintf("\n6.4\nThe root is approximately at x = 0.2\n");
hold off;

%%6.14
g =@(y) (y./(1-y)).*(sqrt(6./(2+y)))-0.05;
root2 = fzero(g,0);
fprintf("\n6.14\nThe value that satisfies the equation is x = %f\n\n", root2);