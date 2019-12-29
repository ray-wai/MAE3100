%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 9, 10/25/19

%%Using composite trapezoidal rule to estimate an integral

x = [0 400 800 1200 1600 2000];
p = [4 3.95 3.8 3.6 3.41 3.3];
a = [100 103 110 120 133 150];
PxA = zeros(1,length(x));

for i = 1:length(x)
    PxA(i) = p(i).*a(i); %%an array of p times a (f(x))
end

midsum = 0;
for j = 2:(length(x) - 1)
    midsum = midsum + PxA(j); %%summing f(x2) to f(x(n-1))
end

h = (x(length(x)) - x(1))./(length(x)-1); %%width of each trapezoid
m = ((x(length(x)) - x(1)).*(PxA(1) + (2.*midsum) + PxA(length(x))))./(2.*(length(x)-1)); %%comp trap formula
fprintf("\nThe value of the mass using composite trapezoidal rule is %d grams\n\n", m);