%dynamics
%ppark

clear all
close all
clc

%Initial Conditions
x1_0=0.2;  %rad
x2_0=0.3;  %rad
dx1_0=0.4; %rad/s
dx2_0=0.5; %rad/s
x0=[x1_0;x2_0;dx1_0;dx2_0];

%constants
m1=2;   %kg
m2=3;   %kg
l1=0.8; %m
l2=0.6; %m
g=9.81; %m/s^2 (gravity)

%Time span
tspan=[0 20]; % 0 - 20 sec

%ODE function
[t,y]=ode45(@(t,x) odefun(t,x, m1,m2,l1,l2,g), tspan, x0);

%Plotting
figure(1)
clf
plot(t,y(:,1))
xlabel('Time (seconds)')
ylabel('Angle (radians)')
title('MAE 511 Problem 2 \theta_1(t)')

figure(2)
clf
plot(t,y(:,2))
xlabel('Time (seconds)')
ylabel('Angle (radians)')
title('MAE 511 Problem 2 \theta_2(t)')
