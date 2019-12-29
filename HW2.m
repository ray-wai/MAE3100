%%Raymond Waidmann
%%18157816 rcw5k2
%%MAE 3100 HW2 9/4/2019

%%iterating to find a root using the false position method

xl = zeros(1,3);
xu = zeros(1,3);
xr = zeros(1,3);
xlxr = zeros(1,3);
iteration = zeros(1,3);
xl(1) = 0.5;
xu(1) = 2;

%%a for loop that loops 3 tiumes and finds xr
for i = 1:3
    iteration(i) = i;
    xr(i) = xu(i)-(((log(xu(i).^2)-0.7).*(xl(i)-xu(i)))/((log(xl(i).^2)-0.7)-(log(xu(i).^2)-0.7)));
    xlxr(i) = (log(xl(i).^2)-0.7).*(log(xr(i).^2)-0.7); %%f(xl)*f(xr)
    if i == 3 %%break out of the for loop on the last iteration
        break;
    end
    if xlxr(i) < 0 %% if root between xl,xr
        xl(i+1) = xl(i);
        xu(i+1) = xr(i);
    else %%if root between xr,xu
        xl(i+1) = xr(i);
        xu(i+1) = xu(i);
    end 
end

%%a 3x3 matrix of function values, by column: xl,xu,xr
fx = zeros(3,3);
for j = 1:3
    fx(1,j) = (log(xl(j).^2)-0.7);
    fx(2,j) = (log(xu(j).^2)-0.7);
    fx(3,j) = (log(xr(j).^2)-0.7);
end

relapproxerror = zeros(1,3);
relapproxerror(1) = 0;
for k = 2:3
    relapproxerror(k) = abs((xr(k)-xr(k-1))/(xr(k))).*100;
end

%%displaying the results
disp(' Iteration    xl       f(xl)      xu        f(xu)     xr        f(xr)    xl*xr      RelApproxError(%)'); %%displaying step size 1 results
disp([iteration(:), xl(:), fx(1,:)', xu(:), fx(2,:)', xr(:), fx(3,:)', xlxr(:), relapproxerror(:)]);
    
    

