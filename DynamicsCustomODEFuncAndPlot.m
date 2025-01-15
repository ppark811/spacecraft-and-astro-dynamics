%system dynamics
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
%q0=[0;0;0;2;0;9.81;0;0;0;4*pi;0.1;0.1];part='Part i'; %part i
q0=[0;0;0;2;0;9.81;0;0;0;0.1;4*pi;0.1];part='Part ii'; %part ii

%---ODE function---%
[t,q] = ode45(@(t,q) odefun(t,q,a,b,c,g), tspan, q0);

    %-x-%
figure(1);clf
plot(t,q(:,1))
title(part,' x value')
xlabel('time (s)')
ylabel('distance (m)')

    %-y-%
figure(2);clf
plot(t,q(:,2))
title(part,' y value')
xlabel('time (s)')
ylabel('distance (m)')

    %-z-%
figure(3);clf
plot(t,q(:,3))
title(part,' z value')
xlabel('time (s)')
ylabel('distance (m)')

    %-phi-%
figure(7);clf
plot(t,q(:,7))
title(part,' \phi value')
xlabel('time (s)')
ylabel('position (rad)')

    %-theta-%
figure(8);clf
plot(t,q(:,8))
title(part,' \theta value')
xlabel('time (s)')
ylabel('position (rad)')

    %-psi-%
figure(9);clf
plot(t,q(:,9))
title(part,' \psi value')
xlabel('time (s)')
ylabel('position (rad)')

    %-wx-%
figure(10);clf
plot(t,q(:,10))
title(part,' \omega_x value')
xlabel('time (s)')
ylabel('Speed (rad/s)')

    %-wy-%
figure(11);clf
plot(t,q(:,11))
title(part,' \omega_y value')
xlabel('time (s)')
ylabel('Speed (rad/s)')

    %-wz-%
figure(12);clf
plot(t,q(:,12))
title(part,' \omega_z value')
xlabel('time (s)')
ylabel('Speed (rad/s)')



function [dqdt]=odefun(t,q,a,b,c,g)
    
    Ixx=c^2 + a^2;
    Iyy=b^2 + a^2;
    Izz=c^2 + b^2;

    q1=q(1);q2=q(2);q3=q(3);q4=q(4);q5=q(5);q6=q(6);
    q7=q(7);q8=q(8);q9=q(9);q10=q(10);q11=q(11);q12=q(12);

%%%%% Was working on this part to flip function to gather angles, but could
%%%%% not figure out how it was supposed to be written... 

%     for i=0:length(t)-1
%         LIA1 = ismembertol(pi/2,q8,0.0872665);
%         LIA2 = ismembertol(-pi/2,q8,0.0872665);
%         LIA3 = ismembertol(3*pi/2,q8,0.0872665);
%         LIA4 = ismembertol(-3*pi/2,q8,0.0872665);
%         test=norm(double([LIA1, LIA2, LIA3, LIA4]));
%         if test~=0
%             
% 
%         else
%             dq8=q11*cos(q7)-q12*sin(q7);                    %theta
%             dq7=sec(q8)*(q12*cos(q7)+q11*sin(q7));          %phi
%             dq9=q10+q11*sin(q7)*tan(q8)+q12*cos(q7)*tan(q8);%psi
%             dq10=(1/Ixx)*((Iyy-Izz)*q11*q12);               %wx
%             dq11=(1/Iyy)*((Izz-Ixx)*q12*q10);               %wy
%             dq12=(1/Izz)*((Iyy-Ixx)*q10*q11);               %wz
%         end
%     end

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

    dqdt=[dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8;dq9;dq10;dq11;dq12];

end
