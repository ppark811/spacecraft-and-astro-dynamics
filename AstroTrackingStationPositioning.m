%ppark
%astrodynamics

clear;clc

format long 

%----------Case 1----------%

    %-----Initial-----%

    %-----Givens-----%

    %-Tracking site-%
    
GST = 265;                      %5:40pm
Long = 64.7;                    %degrees west
L = 76.53;                      %degrees latitude
h = 1.207;                      %km

    %-Satellite data-%
    
rho = 3000;                     %km
rhodot = 6;                     %km/s
Az = 7.5;                       %degrees
Azdot = 0.0174533;              %rad/sec
El = 85;                        %degrees
Eldot = 0.000174533;            %rad/sec
theta = ((17+(40/60))*15)-Long; %LST degrees

TOF = 1.21e6;                   %seconds (2 weeks)

RE = 6378.145;                  %radius of earth km
e = 0.08182;                    %eccentricity of earth
w = [0;0;7.292115856e-5];       %rad/sec

    %-Position and Velocity-%
    
rhovec = [-rho*cosd(El)*cosd(Az) ; rho*cosd(El)*sind(Az) ; rho*sind(El)];
rhodotvec = [-rhodot*cosd(El)*cosd(Az)+rho*sind(El)*Eldot*cosd(Az)+...
             rho*cosd(El)*sind(Az)*(Azdot) ; 
             rhodot*cosd(El)*sind(Az)-rho*sind(El)*(Eldot)*sind(Az)+...
             rho*cosd(El)*cosd(Az)*Azdot ; 
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
    
mu=3.986012e5;                  %km^3/s^2
r=RIJK;
v=VIJK;

h=cross(r,v);                   %angular momentum

rmag=norm(r);                   %radius magnitude
vmag=norm(v);                   %velocity magnitude
hmag=norm(h);                   %angular momentum magnitude

n=[-h(2) h(1) 0];               %node vector
nmag=norm(n);                   %node magnitude

e=(1/mu)*(((vmag^2)-(mu/rmag))*r-(dot(r,v)*v)); %eccentricity vector
emag=norm(e);                                   %eccentricity magnitude

Energy=((vmag^2)/2)-(mu/rmag);
a=-mu/(2*Energy);

        %-----Finding angles-----%

i=acosd(h(3)/hmag);                 %Inclination

if n(2)>0                           %Right ascention of the ascending node
    RAAN=acosd(n(1)/nmag);
else
    RAAN=360-acosd(n(1)/nmag);
end

if e(3)<0                           %Argument of perigee
    ARGP=360-acosd(dot(n,e)/(nmag*emag));
else
    ARGP=acosd(dot(n,e)/(nmag*emag));    
end

if dot(r,v)>0                       %True anomoly
    TA=acosd(dot(r,e)/(rmag*emag));
else                                   
    TA=360-acosd(dot(r,e)/(rmag*emag));
end

        %-----Description of orbit-----%
        
disp('        Case 1 values '); fprintf('\n Initial Points: \n')
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

        %-----displaying initial values-----%

fprintf('\n')
CASE1 = [RIJK VIJK];
T1 = array2table(CASE1,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
COE1 = [a;emag;i;RAAN;ARGP;TA];
ACOE1=[ARGLAT;LONGPER;TRULON];
COORDS1=[L;Long;GST;theta];
SATEPOS1=[rho;Az;El];
SATEVEL1=[rhodot;Azdot;Eldot];
disp(T1)
fprintf('\n   COE Values of Spacecraft: \n')
fprintf('a = %f, e = %f, i = %f, RAAN = %f, ARGP = %f, TA = %f \n \n',COE1)
fprintf('   Alternate COE \n')
fprintf('ARGLAT = %f, LONGPER = %f, TRULON = %f, \n \n',ACOE1)
fprintf('   Coordinates of tracking site: \n')
fprintf('Latitude = %f°, Longitude = %f° west, GST = %f°, LST = %f° \n\n',COORDS1)
fprintf('   Satellite Location: \n')
fprintf('Range = %fkm, Azimuth = %f°, Elevation = %f° \n\n', SATEPOS1)
fprintf('   Satellite Movement: \n')
fprintf('Range rate of change = %fkm/s, Azimuth rate of change = %frad/sec, Elevation rate of change = %frad/sec \n\n',SATEVEL1)
fprintf('--------------------------------------------------------------\n')

    %-----Final-----%
    
n = sqrt(mu/a^3);                                 %mean anomaly
E = zeros(1,1000000);

diff = 1000000; it=2;
while abs(diff)>0.001 && it<1000001
    E(it) = E(it-1) + ((n*TOF)-(E(it-1)-emag*sind(E(it-1))))/(1-(emag*cosd(E(it-1))));
    
        if it < 1000000               %this if else statement is to 
        diff = abs(E(it+1)-E(it));    %get the E vector to be exactly
        else                          %1 milliion elements long
        diff = 0;
        end
        
    it = it+1;
end

Ef = E(1000000)*(180/pi);                              %back into degrees
TAf = acosd((cosd(Ef)-emag)/(1-(emag*cosd(Ef))));      %Final true anomaly

rfmag = (a*(1-emag^2))/(1+emag*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                   %in PQW frame
vPQWf = sqrt(mu/(a*(1-emag^2))) * [ -sind(TAf) ; emag+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAAN)*cosd(ARGP)-sind(RAAN)*sind(ARGP)*cosd(i) -cosd(RAAN)*sind(ARGP)-sind(RAAN)*cosd(ARGP)*cosd(i) sind(RAAN)*sind(i);   
           sind(RAAN)*cosd(ARGP)+cosd(RAAN)*sind(ARGP)*cosd(i) -sind(RAAN)*sind(ARGP)+cosd(RAAN)*cosd(ARGP)*cosd(i) -cosd(RAAN)*sind(i);
           sind(ARGP)*sind(i) cosd(ARGP)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

    %-----Displaying Final results-----%

fprintf('\n Final Points: \n\n')
CASE1f = [RIJKf VIJKf];
T1f = array2table(CASE1f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T1f)
SATEPOS1f = [TAf;rhomagf;Azf;Elf];
fprintf('\n Position of tracking site constant')
fprintf('\n TA = %f°, Range = %fkm, Azimuth = %f°, Elevation = %f° \n\n', SATEPOS1f)
fprintf('=================================================================\n\n')

%----------Case 2----------%

    %-----Initial-----%

    %-----Givens-----%

    %-Tracking site-%
    
GST = 0;                        %10:00pm
Long = 104.54;                  %degrees west
L = 38.8;                       %degrees latitude
h = 1.915;                      %km

    %-Satellite data-%
    
rho = 2121.4180;                %km
rhodot = -3.32040;              %km/s
Az = 350;                       %degrees
Azdot = -0.07653*(pi/180);      %rad/sec
El = 35.3507;                   %degrees
Eldot = 0.20367*(pi/180);       %rad/sec
theta = ((10+(00/60))*15)-Long; %LST degrees

TOF = 345600;                   %seconds (4 days)

RE = 6378.145;                  %radius of earth km
e = 0.08182;                    %eccentricity of earth
w = [0;0;7.292115856e-5];       %rad/sec

    %-Position and Velocity-%
    
rhovec = [-rho*cosd(El)*cosd(Az) ; rho*cosd(El)*sind(Az) ; rho*sind(El)];
rhodotvec = [-rhodot*cosd(El)*cosd(Az)+rho*sind(El)*Eldot*cosd(Az)+...
             rho*cosd(El)*sind(Az)*(Azdot) ; 
             rhodot*cosd(El)*sind(Az)-rho*sind(El)*(Eldot)*sind(Az)+...
             rho*cosd(El)*cosd(Az)*Azdot ; 
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
    
mu=3.986012e5;                  %km^3/s^2
r=RIJK;
v=VIJK;

h=cross(r,v);                   %angular momentum

rmag=norm(r);                   %radius magnitude
vmag=norm(v);                   %velocity magnitude
hmag=norm(h);                   %angular momentum magnitude

n=[-h(2) h(1) 0];               %node vector
nmag=norm(n);                   %node magnitude

e=(1/mu)*(((vmag^2)-(mu/rmag))*r-(dot(r,v)*v)); %eccentricity vector
emag=norm(e);                                   %eccentricity magnitude

Energy=((vmag^2)/2)-(mu/rmag);
a=-mu/(2*Energy);

        %-----Finding angles-----%

i=acosd(h(3)/hmag);                 %Inclination

if n(2)>0                           %Right ascention of the ascending node
    RAAN=acosd(n(1)/nmag);
else
    RAAN=360-acosd(n(1)/nmag);
end

if e(3)<0                           %Argument of perigee
    ARGP=360-acosd(dot(n,e)/(nmag*emag));
else
    ARGP=acosd(dot(n,e)/(nmag*emag));    
end

if dot(r,v)>0                       %True anomoly
    TA=acosd(dot(r,e)/(rmag*emag));
else                                   
    TA=360-acosd(dot(r,e)/(rmag*emag));
end

        %-----Description of orbit-----%
        
disp('        Case 2 values '); fprintf('\n Initial Points: \n')
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

        %-----displaying initial values-----%

fprintf('\n')
CASE2 = [RIJK VIJK];
T2 = array2table(CASE2,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
COE2 = [a;emag;i;RAAN;ARGP;TA];
ACOE2=[ARGLAT;LONGPER;TRULON];
COORDS2=[L;Long;GST;theta];
SATEPOS2=[rho;Az;El];
SATEVEL2=[rhodot;Azdot;Eldot];
disp(T2)
fprintf('\n   COE Values of Spacecraft: \n')
fprintf('a = %f, e = %f, i = %f, RAAN = %f, ARGP = %f, TA = %f \n \n',COE2)
fprintf('   Alternate COE \n')
fprintf('ARGLAT = %f, LONGPER = %f, TRULON = %f, \n \n',ACOE2)
fprintf('   Coordinates of tracking site: \n')
fprintf('Latitude = %f°, Longitude = %f° west, GST = %f°, LST = %f° \n\n',COORDS2)
fprintf('   Satellite Location: \n')
fprintf('Range = %fkm, Azimuth = %f°, Elevation = %f° \n\n', SATEPOS2)
fprintf('   Satellite Movement: \n')
fprintf('Range rate of change = %fkm/s, Azimuth rate of change = %frad/sec, Elevation rate of change = %frad/sec \n\n',SATEVEL2)
fprintf('--------------------------------------------------------------\n')

    %-----Final-----%
    
n = sqrt(mu/a^3);                                 %mean anomaly
Ei = acosd((emag-cosd(TA))/(1+emag*cosd(TA)));    %initial eccentric anomaly
Ei = Ei*(pi/180);                                 %in radians
Mi = Ei-emag*sind(Ei);                            %initial mean motion
Mf = Mi + (n*TOF) - (2*pi);                       %final mean motion

E = zeros(1,1000000);

diff = 1000000; it=2;
while abs(diff)>0.001 && it<1000001
    E(it) = E(it-1) + ((n*TOF)-(E(it-1)-emag*sind(E(it-1))))/(1-(emag*cosd(E(it-1))));
    
        if it < 1000000               %this if else statement is to 
        diff = abs(E(it+1)-E(it));    %get the E vector to be exactly
        else                          %1 milliion elements long
        diff = 0;
        end
        
    it = it+1;
end

Ef = E(1000000)*(180/pi);                              %back into degrees
TAf = acosd((cosd(Ef)-emag)/(1-(emag*cosd(Ef))));      %Final true anomaly

rfmag = (a*(1-emag^2))/(1+emag*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                   %in PQW frame
vPQWf = sqrt(mu/(a*(1-emag^2))) * [ -sind(TAf) ; emag+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAAN)*cosd(ARGP)-sind(RAAN)*sind(ARGP)*cosd(i) -cosd(RAAN)*sind(ARGP)-sind(RAAN)*cosd(ARGP)*cosd(i) sind(RAAN)*sind(i);   
           sind(RAAN)*cosd(ARGP)+cosd(RAAN)*sind(ARGP)*cosd(i) -sind(RAAN)*sind(ARGP)+cosd(RAAN)*cosd(ARGP)*cosd(i) -cosd(RAAN)*sind(i);
           sind(ARGP)*sind(i) cosd(ARGP)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

    %-----Displaying Final results-----%

fprintf('\n Final Points: \n\n')
CASE2f = [RIJKf VIJKf];
T2f = array2table(CASE2f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T2f)
SATEPOS2f = [TAf;rhomagf;Azf;Elf];
fprintf('\n TA = %f, Range = %fkm, Azimuth = %f°, Elevation = %f° \n\n', SATEPOS2f)

    %-----finding optimal site location-----%
    
    
L = 76;                         %degrees latitude north
Long = 68;                      %degrees longditude west
theta = ((10+(00/60))*15)-Long; %LST degrees
h = 0.077;                      %km elevation
e = 0.08182;                    %eccentricity of earth

R = [ ((RE/(sqrt(1-(e^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((RE/sqrt(1-(e^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((RE*(1-e^2))/sqrt(1-(e^2)*(sind(L)^2)))+h)*sind(L)];

rhof = RIJKf - R;               %distance from site to satellite
rhomagf = norm(rhof);           %range
Elf = asind(rhof(3)/rhomagf);   %elevtion from site to satellite
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %azimuth from site to satellite
sitefinal = [rhomagf;Elf;Azf];

    %-----Displaying site data-----%
    
fprintf('Suitable tracking station to find satellite\n')
fprintf('is Thule Tracking Station (TTS), Thule Air Base Greenland\n')
fprintf('Range = %fkm, Elevation = %f°, Azimuth = %f° \n\n',sitefinal)
