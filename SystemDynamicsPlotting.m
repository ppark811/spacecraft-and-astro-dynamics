%system dynamics
%ppark

clear 
close all
clc

%---Constants---%
mp = 2;     %kg
k = 18;     %N/m
lo = 0.5;   %m
g = 9.81;   %m/s/s (gravity
%Distance value - > since ddx=0, the x value won't be used

%---Timespan---%
tspan = [0 40]; %0 - 40 sec

%---Initial Conditions---%
l_0 = 0.75; %m
dl_0 = 0;   %m/s
l0 = [l_0;dl_0];

%---ODE function---%
[t,y] = ode45(@(t,y) odefun(t,y,mp,k,lo,g), tspan, l0);

%---Plotting---%

figure(1)
clf
plot(t,y(:,1))
xlabel('Time (seconds)')
ylabel('Angle (radians)')

%Had to handjam characters a->o in this next part
title('Part(o)')

function [dydt]=odefun(t,y,mp,k,lo,g)%,theta,dtheta)

vectheta=[0, 0.5*t ,t ,1.5*t ,2*t ,2.5*t ,2.9*t ,2.99*t ,2.999*t ,3.0*t ...
    ,3.001*t ,3.01*t ,0.03*t^2 ,0.04*t^2 ,0.05*t^2];
vecdtheta=[0,0.5,1,1.5,2,2.5,2.9,2.99,2.999,3.0,3.001,3.01,0.06*t,0.08*t...
    ,0.1*t];

%for i = 1:length(vectheta) --- Tried to make a loop to make this as 
%                               seemless as possible, but couldn't find 
%                               a way to hold dydt as a variable that 
%                               would held onto? 
%                               global variable i wont hold into a function 

    %theta=vectheta(i);
    %dtheta=vecdtheta(i);
%end

%had to hand-jam values from 1->15
    theta=vectheta(15);
    dtheta=vecdtheta(15);

    y1=y(1);
    y2=y(2);
    dy1=y2;
    dy2=(k*lo - mp*g*cos(theta) - (k - mp*dtheta^2)*y1);
    dydt=[dy1;dy2];

end
