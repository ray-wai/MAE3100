%%Appendix E

%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Final Project Golden Section Search for Min Function

%%This function accepts an anonymous function, and lower and upper bounds.
%%It then uses golden section search to find the input value where the
%%minimum is located and returns that value to the calling program. 

function [xOpt] = GoldSecSearchMin_FS19(f,xl,xu)
phi = ((1+sqrt(5))./2); %%golden ratio
es = 0.1; %%stop condition for RelAppErr
xA = 0;
xB = 0;
xOpt = 0;
d = 0;
ea = 100;

while(ea > es)
    d = (phi-1).*(xu - xl);
    xA = xl + d; %%calculating the two intermediate points
    xB = xu - d;
    
    if (f(xA) < f(xB)) %%positive needed since we are finding min
        xl = xB; %%new lower bound
        xu = xu; %%upper bound stays the same
        xOpt = xA; %%optimum guess
    else
        xl = xl; %%lower bound stays the same
        xu = xA; %%new upper bound
        xOpt = xB; %%optimum guess
    end
    ea = (2-phi).*abs((xu-xl)/(xOpt)).*100; %%RelAppErr
end
end