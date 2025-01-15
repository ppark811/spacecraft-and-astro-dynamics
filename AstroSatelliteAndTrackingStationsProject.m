%ppark
%astrodynamics

clear;clc
format long

mu = 3.986004418e5;             %km^3/s^2
J2 = 1.08262668e-3;             %Zonal harmonic coeffcient
Re = 6378.137;                  %radius of earth (km)
e = 0.08182;                    %eccentricity of earth



%----------Part a----------%




%-----Initial points-----%

ndot = 0.00003156*(2*pi)/(86400^2);%drag rate (rad/sec)
i = 51.6434;                    %inclination (degrees)
RAANi = 348.9165;               %right ascention of the ascending node (degrees)
ei = 0.0006992;                 %ecentricity
ARGPi = 60.4905;                %argument of perigee (degrees)
M = 51.3661*pi/180;             %mean anomaly (radians)
n = (15.54103949*(2*pi))/86400; %mean motion (radians/s)
TOF = 60*60*24*14;              %time of flight (seconds)

ai = (mu/n^2)^(1/3);            %semi-major axis (km)
po = ai*(1-ei^2);               %semi-latus rectum (km)

%-True anomaly loop-%

Mi = M;
Ei = M;

for it = 1:10
    M = Ei - ei*sin(Ei);
    Ei = Ei + (Mi-M)/(1-ei*cos(Ei));
end

TAi = (180/pi)*acos((cos(Ei)-ei)/(1-(ei*cos(Ei))));      %initial true anomaly

%---Rate of change---%

adot = -(2*ndot*ai)/(3*n);                                          %km/s
edot = -(2/3)*(1-ei)*(ndot/n);                                      %rad/s
nbar = (1+(3/2)*J2*((Re/po)^2)*sqrt(1-ei^2)*(1-(3/2)*sind(i)^2))*n; %rad/s
RAANdot = (-(3/2)*J2*(Re/po)^2*cosd(i))*nbar;                       %rad/s
ARGPdot = ((3/2)*J2*(Re/po)^2*(2-(5/2)*sind(i)^2))*nbar;            %rad/s

%-----Final points-----%

af = ai + adot*TOF;             %semi-major axis (km)
ef = ei + edot*TOF;             %eccentricity
RAANf = RAANi + (RAANdot*TOF)*(180/pi);    %Right ascention of ascending node 
ARGPf = ARGPi + (ARGPdot*TOF)*(180/pi);    %Argument of perigee

Mf = Mi + n*TOF + ndot*(TOF^2);
Ef = Mf;
for it = 1:10
    M = Ef - ef*sin(Ef);
    Ef = Ef + (Mf-M)/(1-ei*cos(Ef));
end

TAf = (180/pi)*acos((cos(Ef)-ef)/(1-(ef*cos(Ef))));      %final true anomaly

if Ef > pi
    TAf = 360 - TAf;
end

COEai = [ai ei i RAANi ARGPi TAi];
COEaf = [af ef i RAANf ARGPf TAf];

%-displaying points-%

fprintf('Case 1: ISS \n')
fprintf('Initial Points \n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN = %f \n ARGP = %f \n TA = %f \n', COEai)
fprintf('\nFinal Points\n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN = %f \n ARGP = %f \n TA = %f \n', COEaf)

%-position and velocity-%

%tracking site at Schriever AFB

L = 38.806;             %Latitude (degrees)
theta = 22*15-104.54;    %Longitude (degrees east)
h = 1.915;          %altitude (km)

R = [ abs((Re/sqrt(1-e^2*sind(L)^2))+h)*cosd(L)*cosd(theta);
    abs((Re/sqrt(1-e^2*sind(L)^2))+h)*cosd(L)*sind(theta);
    (((Re*(1-e^2))/sqrt(1-(e^2)*(sind(L)^2)))+h)*sind(L)];

rfmag = (af*(1-ef^2))/(1+ef*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                %in PQW frame
vPQWf = sqrt(mu/(af*(1-ef^2))) * [ -sind(TAf) ; ef+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAANf)*cosd(ARGPf)-sind(RAANf)*sind(ARGPf)*cosd(i) -cosd(RAANf)*sind(ARGPf)-sind(RAANf)*cosd(ARGPf)*cosd(i) sind(RAANf)*sind(i);   
           sind(RAANf)*cosd(ARGPf)+cosd(RAANf)*sind(ARGPf)*cosd(i) -sind(RAANf)*sind(ARGPf)+cosd(RAANf)*cosd(ARGPf)*cosd(i) -cosd(RAANf)*sind(i);
           sind(ARGPf)*sind(i) cosd(ARGPf)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

%-Displaying position-%

fprintf('\n Tracking site: Schriver AFB ')
fprintf('\n Satellite and Site Points: \n\n')
CASE1f = [RIJKf VIJKf];
T1f = array2table(CASE1f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T1f)
SATEPOS1f = [rhomagf;Azf;Elf];
fprintf('\n Range = %fkm \n Azimuth = %f° \n Elevation = %f° \n\n', SATEPOS1f)
fprintf('================================================================\n\n')



%----------Part B----------%



ndot = 0.0000084*(2*pi)/(86400^2);%drag rate (rad/sec)
i = 97.8397;                    %inclination (degrees)
RAANi = 143.5394;               %right ascention of the ascending node (degrees)
ei = 0.0008529;                 %ecentricity
ARGPi = 285.0815;               %argument of perigee (degrees)
M = 74.9399*pi/180;             %mean anomaly (radians)
n = (14.643737148*(2*pi))/86400;%mean motion (radians/s)
TOF = 1.21e6;                   %time of flight (seconds)

ai = (mu/n^2)^(1/3);            %semi-major axis (km)
po = ai*(1-ei^2);               %semi-latus rectum (km)

%-True anomaly loop-%

Mi = M;
Ei = M;

for it = 1:10
    M = Ei - ei*sin(Ei);
    Ei = Ei + (Mi-M)/(1-ei*cos(Ei));
end

TAi = (180/pi)*acos((cos(Ei)-ei)/(1-(ei*cos(Ei))));      %initial true anomaly

%---Rate of change---%

adot = -(2*ndot*ai)/(3*n);                                          %km/s
edot = -(2/3)*(1-ei)*(ndot/n);                                      %rad/s
nbar = (1+(3/2)*J2*(Re/po)^2*sqrt(1-ei^2)*(1-(3/2)*sind(i)^2))*n;   %rad/s
RAANf = RAANi + (RAANdot*TOF)*(180/pi);    %Right ascention of ascending node 
ARGPf = ARGPi + (ARGPdot*TOF)*(180/pi);    %Argument of perigee

%-----FInal points-----%

af = ai + adot*TOF;             %semi-major axis (km)
ef = ei + edot*TOF;             %eccentricity
RAANf = RAANi + RAANdot*TOF;    %Right ascention of ascending node 
ARGPf = ARGPi + ARGPdot*TOF;    %Argument of perigee

Mf = Mi + n*TOF + ndot*(TOF^2);
Ef = Mf;
for it = 1:10
    M = Ef - ef*sin(Ef);
    Ef = Ef + (Mf-M)/(1-ei*cos(Ef));
end

TAf = (180/pi)*acos((cos(Ef)-ef)/(1-(ef*cos(Ef))));      %final true anomaly

if Ef > pi
    TAf = 360 - TAf;
end

COEbi = [ai ei i RAANi ARGPi TAi];
COEbf = [af ef i RAANf ARGPf TAf];

%-displaying points-%

fprintf('Case 2: EO-1\n')
fprintf('Initial Points \n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEbi)
fprintf('\nFinal Points\n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEbf)

%-position and velocity-%

%POGO greenland

L = 76.151;             %Latitude (degrees)
theta = 291.4;    %Longitude (degrees east)
h = .132;          %altitude (km)

R = [ ((Re/(sqrt(1-(ef^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((Re/sqrt(1-(ef^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((Re*(1-ef^2))/sqrt(1-(ef^2)*(sind(L)^2)))+h)*sind(L)];

rfmag = (af*(1-ef^2))/(1+ef*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                   %in PQW frame
vPQWf = sqrt(mu/(af*(1-ef^2))) * [ -sind(TAf) ; ef+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAANf)*cosd(ARGPf)-sind(RAANf)*sind(ARGPf)*cosd(i) -cosd(RAANf)*sind(ARGPf)-sind(RAANf)*cosd(ARGPf)*cosd(i) sind(RAANf)*sind(i);   
           sind(RAANf)*cosd(ARGPf)+cosd(RAANf)*sind(ARGPf)*cosd(i) -sind(RAANf)*sind(ARGPf)+cosd(RAANf)*cosd(ARGPf)*cosd(i) -cosd(RAANf)*sind(i);
           sind(ARGPf)*sind(i) cosd(ARGPf)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

%-Displaying position-%

fprintf('\n Tracking site: ')
fprintf('\n Satellite and Site Points: \n\n')
CASE2f = [RIJKf VIJKf];
T2f = array2table(CASE2f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T2f)
SATEPOS2f = [rhomagf;Azf;Elf];
fprintf('\n Range = %fkm \n Azimuth = %f° \n Elevation = %f° \n\n', SATEPOS2f)
fprintf('================================================================\n\n')


%----------Part C----------%


ndot = 0.00000136*(2*pi)/(86400^2);%drag rate (rad/sec)
i = 4.7447;                     %inclination (degrees)
RAANi = 58.5438;                %right ascention of the ascending node (degrees)
ei = 0.0000844;                 %ecentricity
ARGPi = 164.9477;               %argument of perigee (degrees)
M = 4.3783*pi/180;              %mean anomaly (radians)
n = (1.002381083*(2*pi))/86400; %mean motion (radians/s)
TOF = 1.21e6;                   %time of flight (seconds)

ai = (mu/n^2)^(1/3);            %semi-major axis (km)
po = ai*(1-ei^2);               %semi-latus rectum (km)

%-True anomaly loop-%

Mi = M;
Ei = M;

for it = 1:10
    M = Ei - ei*sin(Ei);
    Ei = Ei + (Mi-M)/(1-ei*cos(Ei));
end

TAi = (180/pi)*acos((cos(Ei)-ei)/(1-(ei*cos(Ei))));      %initial true anomaly

%---Rate of change---%

adot = -(2*ndot*ai)/(3*n);                                          %km/s
edot = -(2/3)*(1-ei)*(ndot/n);                                      %rad/s
nbar = (1+(3/2)*J2*(Re/po)^2*sqrt(1-ei^2)*(1-(3/2)*sind(i)^2))*n;   %rad/s
RAANf = RAANi + (RAANdot*TOF)*(180/pi);    %Right ascention of ascending node 
ARGPf = ARGPi + (ARGPdot*TOF)*(180/pi);    %Argument of perigee

%-----FInal points-----%

af = ai + adot*TOF;             %semi-major axis (km)
ef = ei + edot*TOF;             %eccentricity
RAANf = RAANi + RAANdot*TOF;    %Right ascention of ascending node 
ARGPf = ARGPi + ARGPdot*TOF;    %Argument of perigee

Mf = Mi + n*TOF + ndot*(TOF^2);
Ef = Mf;
for it = 1:10
    M = Ef - ef*sin(Ef);
    Ef = Ef + (Mf-M)/(1-ei*cos(Ef));
end

TAf = (180/pi)*acos((cos(Ef)-ef)/(1-(ef*cos(Ef))));      %final true anomaly

if Ef > pi
    TAf = 360 - TAf;
end

COEci = [ai ei i RAANi ARGPi TAi];
COEcf = [af ef i RAANf ARGPf TAf];

%-displaying points-%

fprintf('Case 3: METEOSAT-8\n')
fprintf('Initial Points \n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEci)
fprintf('\nFinal Points\n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEcf)

%-position and velocity-%

%GUAM

L = 13.615;             %Latitude (degrees)
theta = 144.8667;         %Longitude (degrees east)
h = .2184;             %altitude (km)

R = [ ((Re/(sqrt(1-(ef^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((Re/sqrt(1-(ef^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((Re*(1-ef^2))/sqrt(1-(ef^2)*(sind(L)^2)))+h)*sind(L)];

rfmag = (af*(1-ef^2))/(1+ef*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                   %in PQW frame
vPQWf = sqrt(mu/(af*(1-ef^2))) * [ -sind(TAf) ; ef+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAANf)*cosd(ARGPf)-sind(RAANf)*sind(ARGPf)*cosd(i) -cosd(RAANf)*sind(ARGPf)-sind(RAANf)*cosd(ARGPf)*cosd(i) sind(RAANf)*sind(i);   
           sind(RAANf)*cosd(ARGPf)+cosd(RAANf)*sind(ARGPf)*cosd(i) -sind(RAANf)*sind(ARGPf)+cosd(RAANf)*cosd(ARGPf)*cosd(i) -cosd(RAANf)*sind(i);
           sind(ARGPf)*sind(i) cosd(ARGPf)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

%-Displaying position-%

fprintf('\n Tracking site: ')
fprintf('\n Satellite and Site Points: \n\n')
CASE3f = [RIJKf VIJKf];
T3f = array2table(CASE3f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T3f)
SATEPOS3f = [rhomagf;Azf;Elf];
fprintf('\n Range = %fkm \n Azimuth = %f° \n Elevation = %f° \n\n', SATEPOS3f)
fprintf('================================================================\n\n')



%----------Case D----------%



ndot = 0.00000046*(2*pi)/(86400^2);%drag rate (rad/sec)
i = 82.9339;                    %inclination (degrees)
RAANi = 157.7326;               %right ascention of the ascending node (degrees)
ei = 0.0032539;                 %ecentricity
ARGPi = 43.7054;               %argument of perigee (degrees)
M = 95.6061*pi/180;             %mean anomaly (radians)
n = (13.729113879*(2*pi))/86400;%mean motion (radians/s)
TOF = 1.21e6;                   %time of flight (seconds)

ai = (mu/n^2)^(1/3);            %semi-major axis (km)
po = ai*(1-ei^2);               %semi-latus rectum (km)

%-True anomaly loop-%

Mi = M;
Ei = M;

for it = 1:10
    M = Ei - ei*sin(Ei);
    Ei = Ei + (Mi-M)/(1-ei*cos(Ei));
end

TAi = (180/pi)*acos((cos(Ei)-ei)/(1-(ei*cos(Ei))));      %initial true anomaly

%---Rate of change---%

adot = -(2*ndot*ai)/(3*n);                                          %km/s
edot = -(2/3)*(1-ei)*(ndot/n);                                      %rad/s
nbar = (1+(3/2)*J2*(Re/po)^2*sqrt(1-ei^2)*(1-(3/2)*sind(i)^2))*n;   %rad/s
RAANf = RAANi + (RAANdot*TOF)*(180/pi);    %Right ascention of ascending node 
ARGPf = ARGPi + (ARGPdot*TOF)*(180/pi);    %Argument of perigee

%-----Final points-----%

af = ai + adot*TOF;             %semi-major axis (km)
ef = ei + edot*TOF;             %eccentricity
RAANf = RAANi + RAANdot*TOF;    %Right ascention of ascending node 
ARGPf = ARGPi + ARGPdot*TOF;    %Argument of perigee

Mf = Mi + n*TOF + ndot*(TOF^2);
Ef = Mf;
for it = 1:10
    M = Ef - ef*sin(Ef);
    Ef = Ef + (Mf-M)/(1-ei*cos(Ef));
end

TAf = (180/pi)*acos((cos(Ef)-ef)/(1-(ef*cos(Ef))));      %final true anomaly

if Ef > pi
    TAf = 360 - TAf;
end

COEdi = [ai ei i RAANi ARGPi TAi];
COEdf = [af ef i RAANf ARGPf TAf];

%-displaying points-%

fprintf('Case 4: COSMOS 2361\n')
fprintf('Initial Points \n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEdi)
fprintf('\nFinal Points\n')
fprintf('COEs: \n a = %f \n e = %f \n i = %f \n RAAN =%f \n ARGP = %f \n TA = %f \n', COEdf)

%-position and velocity-%

%POGO

L = 76.515;             %Latitude (degrees)
theta = 291.4;         %Longitude (degrees east)
h = .132;             %altitude (km)

R = [ ((Re/(sqrt(1-(ef^2*(sind(L))^2))))+h)*cosd(L)*cosd(theta);
    ((Re/sqrt(1-(ef^2)*(sind(L)^2)))+h)*cosd(L)*sind(theta);
    (((Re*(1-ef^2))/sqrt(1-(ef^2)*(sind(L)^2)))+h)*sind(L)];

rfmag = (af*(1-ef^2))/(1+ef*cosd(TAf));           
rPQWf = [ (rfmag*cosd(TAf)) ; rfmag*sind(TAf) ; 0 ];                   %in PQW frame
vPQWf = sqrt(mu/(af*(1-ef^2))) * [ -sind(TAf) ; ef+cosd(TAf) ; 0 ]; %in PQW frame

Rtrans = [ cosd(RAANf)*cosd(ARGPf)-sind(RAANf)*sind(ARGPf)*cosd(i) -cosd(RAANf)*sind(ARGPf)-sind(RAANf)*cosd(ARGPf)*cosd(i) sind(RAANf)*sind(i);   
           sind(RAANf)*cosd(ARGPf)+cosd(RAANf)*sind(ARGPf)*cosd(i) -sind(RAANf)*sind(ARGPf)+cosd(RAANf)*cosd(ARGPf)*cosd(i) -cosd(RAANf)*sind(i);
           sind(ARGPf)*sind(i) cosd(ARGPf)*sind(i) cosd(i)];
       
RIJKf = Rtrans*rPQWf;       %turning PQW to IJK
VIJKf = Rtrans*vPQWf;       %turning PQW to IJK

rhof = RIJKf - R;           %finding final range vector
rhomagf = norm(rhof);       %magnitude of range vector

Elf = asind(rhof(3)/rhomagf);               %final elevation
Azf = acosd(rhof(1)/(-rhomagf*cosd(Elf)));  %final azimuth

%-Displaying position-%

fprintf('\n Tracking site: ')
fprintf('\n Satellite and Site Points: \n\n')
CASE4f = [RIJKf VIJKf];
T4f = array2table(CASE4f,'RowNames',{'I','J','K'},'VariableNames',{'Position','Velocity'});
disp(T4f)
SATEPOS4f = [rhomagf;Azf;Elf];
fprintf('\n Range = %fkm \n Azimuth = %f° \n Elevation = %f° \n\n', SATEPOS4f)
