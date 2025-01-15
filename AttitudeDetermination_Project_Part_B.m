%AttitudeDetermination
%ppark
%Project Part B


clear;clc
format long

%---Earth and orbital parameters---%
Alt = 5*10^5;
me=5.9736*10^24;
G = 6.6743*10^-11;
Re = 6.378137*10^6;
R=Re+Alt;
Om = sqrt((G*me)/R^3);
OrbPer = ((2*pi)/Om);

%---Initial Conditions and testing timeframe---%

tspan = [0 1.6*OrbPer];
question = input('Enter desired angle in degrees: ');
phi_0 = question*pi/180; %Deg
theta_0 = question*pi/180; %Deg
psi_0 = question*pi/180; %Deg
dpsi_0 = 0;
dtheta_0 = 0;
dphi_0 = 0;
ang_0 = [phi_0  theta_0  psi_0 dphi_0 dtheta_0 dpsi_0];

[t,y]= ode45(@(t,y) odefun_NL(t,y,Om), tspan, ang_0);
[s,x]= ode45(@(s,x) odefun_L(s,x,Om), tspan, ang_0);


%---Plotting---%

    %-Phi-%
figure(1);clf
subplot(3,1,1)
plot(t/OrbPer,rad2deg(y(:,1)))
hold on
grid on
plot(s/OrbPer,rad2deg(x(:,1)))
ylabel('\phi (deg)')

legend({'nonlinear','linear'},'Location','southeast')
title([num2str(question),char(176), ' initial condition on \phi, \theta, and \psi, neglecting GG'])

%-theta-%
subplot(3,1,2)
plot(t/OrbPer,rad2deg(y(:,2)))
hold on
grid on
plot(s/OrbPer,rad2deg(x(:,2)))
ylabel('\theta (deg)')


    %-psi-%
subplot(3,1,3)
plot(t/OrbPer,rad2deg(y(:,3)))
hold on
grid on
plot(s/OrbPer,rad2deg(x(:,3)))
ylabel('\psi (deg)')
xlabel('time (period)')




%---ODE functions---%
function [dydt]=odefun_NL(tspan,y,Om) % Non-linear only equation
    
    y1=y(1); %phi
    y2=y(2); %theta
    y3=y(3); %psi
    y4=y(4); %dphi
    y5=y(5); %dtheta
    y6=y(6); %dpsi

    dy1=y4;
    dy2=y5;
    dy3=y6;
    dy6=(-(3/4)*Om^2*sin(2*y2)*sin(y1) + ...
        (3/4)*Om^2*cos(y3)^2*sin(2*y2)*sin(y1) - (3/4)*Om^2*sin(2*y2)*sin(y1)*sin(y3)^2 - ...
        Om^2*cos(y2)*cos(y1)*sin(2*y3) + ...
        0.5*Om^2*cos(y2)*sin(y1)*sin(2*y3)*tan(y1) + Om*(2*cos(y1)*cos(y3)*sin(y2)+ ...
        4*sin(y1)*sin(y3)+2*cos(2*y2)*sin(y1)*sin(y3))*y6 + 0.5*Om*(...
        10*cos(y3)*sin(y2)*sin(y1)-8*cos(y1)*sin(y3)+2*cos(2*y2)*cos(y1)*sin(y3))*tan(y1)*y6 + ...
        1.5*sin(2*y2)*sin(y1)*y6^2 - 5*tan(y1)*y4*(Om*cos(y3)*sin(y1)-Om*cos(y1)*sin(y2)*sin(y3)+...
        cos(y2)*cos(y1)*y6) - cos(y1)*y5*(-6*Om*cos(y2)*sin(y3)-2*y4-6*sin(y2)*y6) + ...
        0.5*sin(y1)*tan(y1)*y5*(6*Om*cos(y2)*sin(y3)+10*y4+6*sin(y2)*y6) - ...
        2*y4*(Om*cos(y1)*cos(y3)+Om*sin(y2)*sin(y1)*sin(y3)-cos(y2)*sin(y1)*y6)) / ((4*cos(y2)*cos(y1)+4*cos(y2)*sin(y1)*tan(y1)));  

    dy4=(1/24)*(-4*Om^2*cos(y1)*sin(y1) - 4*(Om^2)*cos(y2)^2*cos(y1)*sin(y1) - ...
        12*(Om^2)*cos(y1)*(cos(y3)^2)*sin(y1) + 4*Om^2*cos(y2)^2*cos(y1)*cos(y3)^2*sin(y1) + ...
        4*Om^2*cos(y1)*sin(y2)^2*sin(y1) - 4*Om^2*cos(y1)*cos(y3)^2*sin(y2)^2*sin(y1) + ...
        12*Om^2*cos(y1)*sin(y1)*sin(y3)^2 - 4*Om^2*cos(y2)^2*cos(y1)*sin(y1)*sin(y3)^2 + ...
        4*Om^2*cos(y1)*sin(y2)^2*sin(y1)*sin(y3)^2 + 8*Om^2*cos(y1)^2*sin(y2)*sin(2*y3) - ...
        8*Om^2*sin(y2)*sin(y1)^2*sin(2*y3) + 16*Om*cos(y3)*sin(2*y1)*y5 - 24*Om*sin(y2)*sin(y3)*y5 - ...
        16*Om*cos(2*y1)*sin(y2)*sin(y3)*y5 - 8*sin(2*y1)*y5^2 + 24*Om*cos(y2)*cos(y3)*y6 - ...
        16*Om*cos(y2)*cos(2*y1)*cos(y3)*y6 - 8*Om*sin(2*y2)*sin(2*y1)*sin(y3)*y6 + ...
        24*cos(y2)*y5*y6 + 16*cos(y2)*cos(2*y1)*y5*y6 + 8*cos(y2)^2*sin(2*y1)*y6^2 + 24*sin(y2)*dy6);

    dy5=(1/8)*sec(y1)* (...
        -0.5*Om^2*cos(y1)*sin(2*y2) + 0.5*Om^2*cos(y1)*cos(y3)^2*sin(2*y2) - 0.5*Om^2*cos(y1)*sin(2*y2)*sin(y3)^2 + ...
        Om^2*cos(y2)*sin(y1)*sin(2*y3) + Om*(10*cos(y3)*sin(y2)*sin(y1)-8*cos(y1)*sin(y3)+2*cos(2*y2)*cos(y1)*sin(y3))*y6 + ...
        cos(y1)*sin(2*y2)*y6^2 - 10*y4*(Om*cos(y3)*sin(y1)-Om*cos(y1)*sin(y2)*sin(y3)+cos(y2)*cos(y1)*y6) - ...
        sin(y1)*y5*(-6*Om*cos(y2)*sin(y3)-10*y4-6*sin(y2)*y6) - 8*cos(y2)*sin(y1)*dy6);

    dydt = [dy1;dy2;dy3;dy4;dy5;dy6];
end
function [dxdt]=odefun_L(tspan,x,Om) % linear only equation
   
    x1=x(1); %phi
    x2=x(2); %theta
    x3=x(3); %psi
    x4=x(4); %dphi
    x5=x(5); %dtheta
    x6=x(6); %dpsi

    dx1=x4;
    dx2=x5;
    dx3=x6;
    dx4=(Om/3)*(-2*Om*x1 + x6);
    dx5=0;
    dx6=(-Om/2)*(2*Om*x3 + x4);
    dxdt = [dx1;dx2;dx3;dx4;dx5;dx6];
end
