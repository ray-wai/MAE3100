%%Appendix C

%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Final Project Gauss Elimination Function 
%%Using Guass Elimination to solve a system of equations

%%This function accepts the matrices of the normal equations and returns 
%%a matrix containing the solutions to that system of equations to the
%%calling program. 


function [x] = GaussEliminationFnct_FS19(a,b)
c = [a b];
n = length(b); %%size of matrix

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
end


