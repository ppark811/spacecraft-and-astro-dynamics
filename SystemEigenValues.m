%system dynamics
%ppark

clear 
close all
clc

%===================%
%-----Problem 1-----%
%===================%

IcmB=[3 -0.5 -0.7;
     -0.5 4 -0.9;
     -0.7 -0.9 8];

[V,D] = eig(IcmB);

E1=V(:,1);
E2=V(:,2);
E3=-V(:,3); %This is set negative bc cross(e1,e2)=-e3. Setting this 
            %negative is okay, as the eigenvector remains an eigenvector
I1=D(1,1);
I2=D(2,2);
I3=D(3,3);

%==================%
%----Problem 2c----%
%==================%

IcA=[9.4367 0 0;
     0 1.7692 0.0216;
     0 0.0216 0.0875];

[v,d] = eig(IcA);
e2=v(:,1);
e1=v(:,2);
e3=v(:,3);

i2=d(1,1);
i1=d(2,2);
i3=d(3,3);
