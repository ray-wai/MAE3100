%%Raymond Waidmann 
%%MAE 3100, Lab 5, 9/27/19
%%Using Guass Elimination to solve a system of equations

%%given parameters
a=[5 2 -1 1;4 -2 3 1;3 1 2 -1;2 3 5 7];
b=[8;4;-15;8];

c = [a b];
n = 4; %%size of matrix (5x5)

%%forward elimination
%%outer for loop (i) is for the column, inner for loop (j) is for rows
%%for a given column, divide the coeff's of a row by the corresponding
    %%coeff of the pivot equation for that iteration (pivot eq = the
    %%equation in row i) to find the multiplier
%%using the multiplier, subtract (multiplier*pivot row coeff's) from the
    %%current value in the non pivot rows. 
%%at the end of these nested for loops, triangular form is achieved. 
for i=1 : n-1
    for j= i+1 : n
        multiplier = c(j,i)./c(i,i);  %%multiplier for the row
        c(j,i:n+1) = c(j,i:n+1) - multiplier*c(i,i:n+1); %%row = current value - multiplier*row i
    end
end

%%back substitution:
    %%create an array "x" to hold the solutions
    %%with each loop, multiply the known x's by the coeff's and subtract from
        %%the value on the right hand side of the eq. Then divide by the coeff
        %%of the unknown X
    %%note: c(k,k+1:n)*x(k+1:n) in loop 1 = 0
x = zeros(n,1);
for k=(n) : -1 : 1
    x(k) = (c(k,n+1) - (c(k,k+1:n)*x(k+1:n)))./c(k,k); 
end

%%displaying the results
disp ("Using Gauss Elimination to Solve");
fprintf("    x1       x2        x3         x4\n");
disp(x');

%%solving and displaying with inverse function
y = inv(a)*b;
disp("Using inv() to solve");
fprintf("    x1       x2        x3         x4\n");
disp(y');


