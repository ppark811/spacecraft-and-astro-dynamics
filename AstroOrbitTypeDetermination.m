%ppark
%astrodynamics

clear;clc

%-----Givens-----%

fprintf('In units of kms and seconds, in column form: \n')
fprintf(' \n')
r=input('enter radius vector: ');
v=input('enter velocity vector: ');
mu=3.986012e5; %km^3/s^2

%using inputs so all cases can be analyzed without
%changing script

%-----Calculations-----%

%these values needed for finding the different COEs

h=cross(r,v);           %angular momentum

rmag=norm(r);           %radius magnitude
vmag=norm(v);           %velocity magnitude
hmag=norm(h);           %angular momentum magnitude

n=[-h(2) h(1) 0];       %node vector
nmag=norm(n);           %node magnitude

e=(1/mu)*(((vmag^2)-(mu/rmag))*r-(dot(r,v)*v)); %eccentricity vector
emag=norm(e);                                   %eccentricity magnitude

E=((vmag^2)/2)-(mu/rmag);
a=-mu/(2*E);

%-----Finding angles-----%

i=acosd(h(3)/hmag);                  %Inclination

if n(2)>0                           %Right ascention of the ascending node
    RAAN=acosd(n(1)/nmag);
else
    RAAN=360-acosd(n(1)/nmag);
end

if e(3)<0                            %Argument of perigee
    ARGP=360-acosd(dot(n,e)/(nmag*emag));
else
    ARGP=acosd(dot(n,e)/(nmag*emag));    
end

if dot(r,v)>0                        %True anomoly
    TA=acosd(dot(r,e)/(rmag*emag));
else                                   
    TA=360-acosd(dot(r,e)/(rmag*emag));
end

%-----Description of orbit-----%

if emag<0.001
    disp('This orbit is circular')
elseif abs(1-emag)<0.001
    disp('This orbit is parabolic')
elseif emag>1
    disp('This orbit is hyperbolic')
else
    disp('This orbit is elliptical')
end

if i<0.001
    disp('This orbit is equatorial')
else
    disp('This orbit is not equatorial')
end



if i<0.001
    longper=RAAN+ARGP;
    COE=[a;emag;i;longper];
    text=sprintf(' a=%skm \n e=%s \n i=%s° \n longper=%s°', COE);
else
    COE=[a;emag;i;RAAN;ARGP;TA];                              
    text=sprintf(' a=%skm \n e=%s \n i=%s° \n RAAN=%s° \n ARGP=%s° \n TA=%s°', COE);  
end
fprintf('The COEs of this orbit are:\n')
disp(text)
fprintf('\n')
