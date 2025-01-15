%MAE 789
%Project - Lagrangian
%Paul Park

clear 
close all
clc

% Physical components
g = 9.81;       % m/s^2

% Pendulum components (already coded into equations)
L = 0.0263;     % m
m = 0.032;      % kg
r1 = -0.0737;   % m
r2 = 0.706;     % m
r3 = 0.0584;    % m

%---Timespan---%
tspan = [0 0.4]; %0 - 0.4 sec
%---Initial Conditions---%
% Order:
% q1, q2, dq1, dq2
% v1, v2, v3
% dv1,dv2,dv3
% phi, theta, psi
% dphi, dtheta, dpsi
%case one has acceleration pulse as 20
q0 = [0.1;0.1;0;0;0;0;0;20;0;0;0;0;0;0;0;0]; 
%---ODE function call---%
[t,q] = ode45(@(t,q) odefun(t,q), tspan, q0);
%---Plotting---%

 %q1 angle
figure(1)
plot(t,q(:,3))
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Angle of q1 vs Time')

 %q2 angle
figure(2)
plot(t,q(:,4))
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Angle of q2 vs Time')

%---ODE function---%
function [dqdt] = odefun(t,q)

    %---Creating time dependent variables:---%
   %q1         q2
    q1 = q(1); q2 = q(2);                   
    
   %dq1 (u1)   dq2 (u2)
    q3 = q(3); q4 = q(4);                   
    
   %v1         v2         v3
    q5 = q(5); q6 = q(6); q7 = q(7);    

   %dv1        dv2        dv3 
    q8 = q(8); q9 = q(9); q10 = q(10);      

   %phi          theta        psi
    q11 = q(11); q12 = q(12); q13 = q(13); 

   %dphi         dtheta       dpsi  
    q14 = q(14); q15 = q(15); q16 = q(16);  
    %---Creating 1st order ODEs---%
    dq1 = q3;
    dq2 = q4;
    dq3 = (-4*csc(q2)/(4+csc(q2)))*((2*cos(q2) + 0.25-cot(q2))*q3*q4 + ...
        q5*(-19.011*cos(q1)*q3 - 19.011*cot(q2)*sin(q1)*q4));
    dq4 = (q5*(76.045*cos(q2)*sin(q1)*q3 + 76.045*cos(q1)*sin(q2)*q4) + ...
        q4*((-0.75*cos(q1) + 0.75*cos((3*q1) + 0.75*sin(q1) + 0.75*sin(3*q1))* ...
        q3)))/(4+0.25*cos(q1)^3 + 0.75*sin(q1) - 0.75*cos(q1)^2*sin(q1) + ...
        0.25*sin(q1)^3 + cos(q1)*(0.75-0.75*sin(q1)^2));
    dq5 = q8;
    dq6 = q9;
    dq7 = q10;
    dq8 = 0;

    dq9 = 0;
    dq10 = 0;
    dq11 = q14;
    dq12 = q15;
    dq13 = q16;
    dq14 = 0;
    dq15 = 0;
    dq16 = 0;

    dqdt=[dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8;dq9;dq10;dq11;dq12;dq13;dq14;dq15;dq16];
end
