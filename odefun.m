%MAE 511 
%ODE function for Homework1_problem2.m, HW set 1 - Problem 2
%9.7.2022
%Paul Park

function [dxdt]=odefun(t,x, m1,m2,l1,l2,g)
   
    x1=x(1);
    x2=x(2);
    x3=x(3);
    x4=x(4);
    
    dx1=x3;
    dx2=x4;
    dx3=(m2*l1*x3^2*sin(x2-x1)*cos(x2-x1) + m2*g*sin(x2)*cos(x2-x1) ... 
        + m2*l2*x4^2*sin(x2-x1) - (m1+m2)*g*sin(x1)) / ((m1+m2)*l1 ...
        - m2*l1*(cos(x2-x1)^2));
    dx4=(-l1*dx3*cos(x2-x1) - l1*x3^2*sin(x2-x1) - g*sin(x2)) / l2;
    dxdt=[dx1;dx2;dx3;dx4];

end 