%dynamics
%ppark

clear;clc

mA=1; mB=2; mC=3;
lA=4; lB=5; lC=6;

m1=0.5; m2=0.6; m3=0.7;
m4=0.8; m5=0.9; m6=1;

mT=mA+mB+mC+m1+m2+m3+m4+m5+m6;

l1=1; l2=2; l3=3;
l4=4;l5=5;l6=6;
l7=3.5;

rcmp=[(m1*l1-m2*l2)/mT;
    (m3*l3-m4*l4+m6*l7)/mT;
    (m5*l5-m6*l6)/mT]

Ipsys = [(mT/12)*(l1-l2)^2+mT*(rcmp(2)^2+rcmp(3)^2) , -mT*rcmp(1)*rcmp(2) , -mT*rcmp(1)*rcmp(3);
    -mT*rcmp(1)*rcmp(2) , (mT/12)*(l3-l4+l7)^2+mT*(rcmp(1)^2+rcmp(3)^2) , -mT*rcmp(2)*rcmp(3);
    -mT*rcmp(1)*rcmp(3) , -mT*rcmp(2)*rcmp(3) , (mT/12)*(l5-l6)^2+mT*((rcmp(1)^2)*rcmp(2)^2)]

