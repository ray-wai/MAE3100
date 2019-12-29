%%Raymond Waidmann
%%18157816, rcw5k2
%%MAE 3100, Lab 11, 11/7/19

%%Predefined functions
func=@(t,y) t.^3-2.*y;
h=0.4;
t=0:h:2;
y(1)=1;

%%Eulers Method
for i=2:length(t)
    y(i)=y(i-1)+func(t(i-1),y(i-1)).*h;
end

plot(t,y);
hold on;
xlabel('t');
ylabel('y');
title('Eulers and Heuns Methods to Solve an ODE');
grid on;

%%Heuns Method
y2(1)=1;
for j=2:length(t)
    ea=1;
    y2(j)=y2(j-1)+func(t(j-1),y2(j-1)).*h;
    while ea>0.1
        temp=y2(j);
        y2(j)=y2(j-1)+((func(t(j-1),y2(j-1))+func(t(j),temp))/2).*h;
        ea=abs((y2(j)-temp)./(y2(j))).*100;
    end
end

plot(t,y2);
legend('Eulers Method','Heuns Method','location','northwest');
hold off;