%dynamics
%ppark


clear 
close all
clc

%---Constants---%
mA=2;       %kg
lA=0.3;     %m
mB=1.5;     %kg
lB=0.25;    %m
C2=0.4;     %rad
g=9.81;     %m/s^2
w2=13;       %rad/s
%dtheta1(0) and theta2(0) won't be used since it's zero

%---Timespan---%
tspan = [0 10]; %0 - 10 sec

%---Initial Conditions---%
The1_0=0.8; %rad
dThe1_0=0;  %rad/s
y0=[The1_0;dThe1_0];

%---ODE function---%
[t,y] = ode45(@(t,y) odefun(t,y,mA,lA,mB,lB,C2,g,w2), tspan, y0);

%---Plotting---%

figure(1)
clf
plot(t,y(:,1))
xlabel('Time (seconds)')
ylabel('Angle (radians)')
title('Part B')

function  [dydt]=odefun(t,y,mA,lA,mB,lB,C2,g,w2)

    the2=w2*C2*cos(w2*t);
    ddthe2=-w2^2*C2*cos(w2*t);
    
    y1=y(1);
    y2=y(2);
    dy1=y2;
    dy2=((((-mA*lA*g)/2)*sin(y1) - mB*lA*g*sin(y1) - ((mB*lB*g)/2)*sin(y1+the2)) / ((mB*lB^2)/12 - (mA*lA^2)/12)) - ddthe2;
    
    dydt=[dy1;dy2];

end
