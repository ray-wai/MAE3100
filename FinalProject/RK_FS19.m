%%Appendix G

%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Final Project Runge Kutta Function

%%*****PART 2.1*****%%

%%this function intakes an ode function, time interval, initial values,
%%time step, and function order and returns a column vector of times and a
%%matrix containing the solutions at each corresponding time. This function
%%only supports 1st, 2nd, and 4th order Runge Kutta

%%Definition of variables:
%%tout: a vector of times at which y is evaluated 
%%yout: a matrix where each row contains y and y' for each x, y, z evaluated at the correspinding tout
%%ode_function_FS19: another function where the derivatives are computed
%%tspan: a vector containing initial and final times (the interval)
%%y0: a vector of initial values of the vector y
%%h: time step
%%rk: the order of Runge Kutta to be used

function [tout,yout] = RK_FS19(Ode_Function_FS19,tspan,y0,h,rk)
tout=[tspan(1):h:tspan(2)]';
yout(1,:)=y0;

if rk==1
    for i=1:length(tout)-1
        yout(i+1,:)=yout(i,:)+(Ode_Function_FS19(tout(i),yout(i,:)).*h)';
    end
     
elseif rk==2
    for i=1:length(tout)-1
        k1=(Ode_Function_FS19(tout(i),yout(i,:)))';
        k2=(Ode_Function_FS19(tout(i)+(3/4).*h,yout(i,:)+(3/4).*k1.*h))';
        yout(i+1,:)=yout(i,:)+((1/3).*k1+(2/3).*k2).*h;
    end
    
elseif rk==4
    for i=1:length(tout)-1
        K1=(Ode_Function_FS19(tout(i),yout(i,:)))';
        K2=(Ode_Function_FS19(tout(i)+(1/2).*h,yout(i,:)+(1/2).*K1.*h))';
        K3=(Ode_Function_FS19(tout(i)+(1/2).*h,yout(i,:)+(1/2).*K2.*h))';
        K4=(Ode_Function_FS19(tout(i)+h,yout(i,:)+K3.*h))';
        yout(i+1,:)=yout(i,:)+(1/6).*(K1+2.*K2+2.*K3+K4).*h;
    end
end
end