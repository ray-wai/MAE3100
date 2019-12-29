%%Raymond Waidmann
%%18157816 rcw5k2
%%MAE 3100 HW5

%%14.4
x=[2 4 6 7 10 11 14 17 20];
y=[4 5 6 5 8 8 6 9 12];
figure(1);
plot(x,y,'ro'); %%plotting the data points
hold on;
a=[0:21];
z = (744./1211).*a; %%hard coded function we calculated
plot(a,z);
axis([0 21 0 13]);
title('14.4');
xlabel('X');
ylabel('Y');
hold off;

%%14.1
v =[10 20 30 40 50 60 70 80];
f =[25 70 380 550 610 1220 830 1450];
figure(2);
plot(v,f,'ro'); %%plotting the data points
hold on;
b=[0:90];
m = (-21253./119) + (13547./840).*b + (1061./28560).*(b.^2); %%hard coded function we calculated
plot(b,m);
title('14.1 Force vs. Velocity');
xlabel('Velocity');
ylabel('Force');
hold off;

