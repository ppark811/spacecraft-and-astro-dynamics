%ppark
%astrodynamics

clear;clc

%-----Given-----%

mu=3.986005e5;

%-----Initial Conditions-----&

a=input('enter semi-major axis: ');
e=input('enter eccentricity: ');

%-----Using givens and equations-----%

TP=2*pi*sqrt(a^3/mu);                % Equation for time period
rp=a*(1-e);                          % Equation for position at perigee
vp=sqrt(mu*((2/rp)-(1/a)));          % Equation for velocity at perigee     

%-----Euler eqn setup-----%

h=1;
r=[rp;0];                            % Setting initial point
v=[0;vp];                            % Setting initial velocity
ac=@(mu,r) -(mu/(norm(r)^3))*r;      % Function for acceleration

if a>10000                           % Making sure that the orbit is fully closed by increasing step
    step=a*(e/4.555);                % 4.555 is an arbritary constant found through testing
else
    step=0;                          % Making it so that if a is less than 10000
end                                  % the step size is zero to ensure no more than one revolution

%-----Loop to create vectors to make points along the orbit-----%

for t=1:TP+step
   r(:,t+h)=r(:,t)+v(:,t);           % Building matrix, with each column being a point
   v(:,t+h)=v(:,t)+ac(mu,r(:,t));    % Vector matrix being built to add into "r" vector   
end
 
%-----Seperating the "r" vector to make x and y-----%

x=r(1,:);                            % First row is x-axis
y=r(2,:);                            % Second row is y-axis

%-----Plotting data-----%

figure(1);clf                       
plot(x,y); grid on; hold on          % The euler step and ode45 graphs will 
    title('Orbital Path')            % be on top of each other to compare easier
    xlabel('x-axis (km)')
    ylabel('y-axis (km)')

%-------------------------------------------------------------------------------------------------%

%-----ODE45 graph-----%

% a, e, Rp ,vp, mu, TP, will remain constant. New variables will override old
% variables or new variable names will be created

y0 = [rp; 0; 0; vp];                 % initial condition will be put into 4x1 matrix 
                                     % so ode45 can work properly
%------ODE45 setup-----%

options = odeset('MaxStep', 3);      % An option to increase the step size within the ode45 function
[t, y] = ode45(@accel, [0, TP],y0,options);     % The ode45 will take from 0 to TP where the y0 is 
                                                % added on and step size is determined
%-----Showing integration-----%

difference = y(size(y),:) - y(1,:);   % Showing integration

fprintf('dx: %16.8f m\n', difference(1));       % Printing results onto command window
fprintf('dy: %16.8f m\n', difference(2));
fprintf('dvx: %16.8f m/sec\n', difference(3));
fprintf('dvy: %16.8f m/sec\n', difference(4));

%-----Plot for ODE45-----%

plot(y(:,1), y(:,2),'r'); grid on
   legend('Euler step', 'ODE45')

%-----Derivative function-----%

function a = accel(t, y)

    mu=3.986005e5;                  % Setting variables
    r = [y(1); y(2)];
    r3 = norm(r)^3;

    g = - mu/r3 * r;                % acceleration from deriving velocity
 
    a = [y(3); y(4); g(1); g(2)];   % y(3) and y(4) are the velocity terms and first derivatives
end
