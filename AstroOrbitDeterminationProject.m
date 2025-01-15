%ppark
%astrodynamics

clear;clc

%-----Case 1-----%

    %-Givens-%

rho = 2121.418;             %km
rhodot = -3.3204;       
Az = 18.9606;               %degree
Azdot = -0.001353153;       %rad/sec
El = 35.3507;               %degree
Eldot = 0.0035547;          %rad/sec
L = 77.7;   
theta = (1.3*(15))-68.6 ;   %LST degrees
GST = 0130;                 %GST
    
RE = 6378.145;
e = 0.08182;
h = 0;
w = [0;0;7.29211e-5];       %rad/sec

    %-Position and Velocity-%
    
rhovec = [-rho*cosd(El)*cosd(Az) ; rho*cosd(El)*sind(Az) ; rho*sind(El)];
rhodotvec = [-rhodot*cosd(El)*cosd(Az)+rho*sind(El)*Eldot*cosd(Az)+...
    rho*cosd(El)*sind(Az)*(Azdot) ; 
    rhodot*cosd(El)*sind(Az)-rho*sind(El)*(Eldot)*sind(Az)+...
    rho*cosd(El)*cos(Az)*Azdot ; 
    rhodot*sind(El)+rho*cosd(El)*Eldot];

T = [sind(L)*cosd(theta) -sind(theta) cosd(L)*cosd(theta);
    sind(L)*sind(theta) cosd(theta) cosd(L)*sind(theta);
    -cosd(L) 0 sind(L)];

R = [ ((RE/(sqrt(1-(e^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((RE/sqrt(1-(e^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((RE*(1-e^2))/sqrt(1-(e^2)*(sind(L)^2)))+h)*sind(L)];

rhoijk = T*rhovec;
rhovecijk = T*rhodotvec;

RIJK = rhoijk+R;
VIJK = rhovecijk+cross(w,RIJK);

    %-COE-%
    
mu=3.986012e5; %km^3/s^2
r=RIJK;
v=VIJK;

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
        
disp('        Case 1 values '); fprintf('\n')
if emag<0.001
    disp('This orbit is circular')
elseif 0.009>emag && emag>0.001
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
        %-----Alternate COE-----%

if r(3)<0        
    ARGLAT=acosd(dot(n,r)/(nmag*rmag));
else
    ARGLAT=360-acosd(dot(n,r)/(nmag*rmag));
end

if e(2)<0
    LONGPER=norm(acosd(e/emag));
else
    LONGPER=360-norm(acosd(e/emag));
end

if r(2)<0
    TRULON=acosd(r(1)/rmag);
else
    TRULON=360-acosd(r(1)/rmag);
end

            %-----Final values and displaying them-----%
            
CASE1 = [RIJK VIJK];
T2 = array2table(CASE1,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
COE1 = [a;emag;i;RAAN;ARGP;TA];
ACOE1=[ARGLAT;LONGPER;TRULON];
COORDS1=[L;GST;theta];
disp(T2)
fprintf('\n COE Values: \n')
fprintf('a = %f, e = %f, i = %f, RAAN = %f, ARGP = %f, TA = %f \n \n',COE1)
fprintf('Alternate COE \n')
fprintf('ARGLAT = %f, LONGPER = %f, TRULON = %f, \n \n',ACOE1)
fprintf('Coordinates: \n')
fprintf('Longitude = %f, GST = %f, LST = %f \n',COORDS1)
fprintf('==============================================================\n')
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----Case 2-----%

    %-Givens-%

rho = 37390;                  %km
rhodot = -4.914;       
Az = 148.25;                  %degree
Azdot = 0.00065537;           %rad/sec
El = 45.32;                   %degree
Eldot = -0.00084229;          %degree/sec
L = 34.15;   
theta = (6.3*(15))-118.15;    %LST
GST = 0630;                   %GST

RE = 6378.145;
e = 0.08182;
h = 0;
w = [0;0;7.29211e-5];   %rad/sec

    %-Position and Velocity-%
    
rhovec = [-rho*cosd(El)*cosd(Az) ; rho*cosd(El)*sind(Az) ; rho*sind(El)];
rhodotvec = [-rhodot*cosd(El)*cosd(Az)+rho*sind(El)*Eldot*cosd(Az)+...
    rho*cosd(El)*sind(Az)*(Azdot) ; 
    rhodot*cosd(El)*sind(Az)-rho*sind(El)*(Eldot)*sind(Az)+...
    rho*cosd(El)*cos(Az)*Azdot ; 
    rhodot*sind(El)+rho*cosd(El)*Eldot];

T = [sind(L)*cosd(theta) -sind(theta) cosd(L)*cosd(theta);
    sind(L)*sind(theta) cosd(theta) cosd(L)*sind(theta);
    -cosd(L) 0 sind(L)];

R = [ ((RE/(sqrt(1-(e^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((RE/sqrt(1-(e^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((RE*(1-e^2))/sqrt(1-(e^2)*(sind(L)^2)))+h)*sind(L)];

rhoijk = T*rhovec;
rhovecijk = T*rhodotvec;

RIJK = rhoijk+R;
VIJK = rhovecijk+cross(w,RIJK);

    %-COE-%
    
mu=3.986012e5; %km^3/s^2
r=RIJK;
v=VIJK;

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
        
disp('        Case 2 values '); fprintf('\n')
if emag<0.001
    disp('This orbit is circular')
elseif 0.009>emag && emag>0.001
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
        %-----Alternate COE-----%

if r(3)<0        
    ARGLAT=acosd(dot(n,r)/(nmag*rmag));
else
    ARGLAT=360-acosd(dot(n,r)/(nmag*rmag));
end

if e(2)<0
    LONGPER=norm(acosd(e/emag));
else
    LONGPER=360-norm(acosd(e/emag));
end

if r(2)<0
    TRULON=acosd(r(1)/rmag);
else
    TRULON=360-acosd(r(1)/rmag);
end

            %-----Final values and displaying them-----%
            
CASE2 = [RIJK VIJK];
T2 = array2table(CASE2,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
COE2 = [a;emag;i;RAAN;ARGP;TA];
ACOE2=[ARGLAT;LONGPER;TRULON];
COORDS2=[L;GST;theta];
disp(T2)
fprintf('\n COE Values: \n')
fprintf('a = %f, e = %f, i = %f, RAAN = %f, ARGP = %f, TA = %f \n \n',COE2)
fprintf('Alternate COE \n')
fprintf('ARGLAT = %f, LONGPER = %f, TRULON = %f, \n \n',ACOE2)
fprintf('Coordinates: \n')
fprintf('Longitude = %f, GST = %f, LST = %f \n',COORDS2)
