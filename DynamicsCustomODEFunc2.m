%dynamics
%ppark

clear 
close all
clc

%---Constants---%
%m = 0.155;    %kg <- mass left out due to diff eqns not having a use for
%it
a = 0.008;    %m
b = 0.068;    %m
c = 0.149;    %m
g = 9.81;     %m/s^2

%---Timespan---%
tspan = [0 2];  %0 - 2 sec

%---Initial Conditions---%
q0=[0;0;0;2;0;9.81;0;0;0;4*pi;0.1;0.1;0;0];


%---ODE function---%
[t,q] = ode45(@(t,q) odefun(t,q,a,b,c,g), tspan, q0);

    %-z vs x-%
figure(1);clf
plot(q(:,14),q(:,13))
title('z-x value')
xlabel('z-value (m)')
ylabel('x-value (m)')




function [dqdt]=odefun(t,q,a,b,c,g)
    
    Ixx=c^2 + a^2;
    Iyy=b^2 + a^2;
    Izz=c^2 + b^2;

    q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);q6=q(6);
    q7=q(7);q8=q(8);q9=q(9);q10=q(10);q11=q(11);q12=q(12);
    q13=q(13);q14=q(14);

    dq1=q4-q11*q3+q12*q2;                   %x
    dq2=q5-q12*q1+q10*q3;                   %y
    dq3=q6-q10*q2+q11*q1;                   %z
    dq4=-q11*q6+q12*q5;                     %u
    dq5=-q12*q4+q10*q6;                     %v
    dq6=-g-q10*q5+q11*q4;                   %w
    dq7=sec(q8)*(q12*cos(q7)+q11*sin(q7));  %phi
    dq8=q11*cos(q7)-q12*sin(q7);            %theta
    dq9=q10+q11*sin(q7)*tan(q8)+q12*cos(q7)*tan(q8); %psi
    dq10=(1/Ixx)*((Iyy-Izz)*q11*q12);       %wx
    dq11=(1/Iyy)*((Izz-Ixx)*q12*q10);       %wy
    dq12=(1/Izz)*((Iyy-Ixx)*q10*q11);       %wz

    dq13=q3*(sin(q7)*sin(q9)+cos(q7)*cos(q9)*sin(q8)) ...
        - q2*(cos(q7)*sin(q9)-cos(q9)*sin(q7)*sin(q8)) ...
        + q1*cos(q9)*cos(q8);               %xo
    dq14=q3*cos(q7)*cos(q8)-q1*sin(q8)+q2*cos(q8)*sin(q7);  %zo

    dqdt=[dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8;dq9;dq10;dq11;dq12;dq13;dq14];

end
