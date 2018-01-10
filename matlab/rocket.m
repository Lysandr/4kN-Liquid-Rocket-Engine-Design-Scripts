clear all;
close all;

%% 8/29/18 the braenig optimal method

%inputss
theta = 0.261799;
mix_ratio = 1.25;
m_dot = 8;      % kg/s, mass flow rate
gam = 1.2;      % specific heat ratio
Tc = 3060;      % flame temperature, Kelvin
Pc = 2758000;   % chamber pressure, MPa
M = 22.9/1000;  % molecular weight km/mol
Pa = 101325;    % atmospheric pressure
R_un = 8.314;   % universal gas constant
Pe = Pa;        % assuming nozzle is designed properly..

% pressure and temperature at the throat
Pt = Pc * ((1+((gam-1)/2))^((-gam)/(gam-1)));
Tt = Tc * ((1+((gam-1)/2))^-1);
% throat area
At = (m_dot/Pt)*sqrt((R_un*Tt)/(M*gam));
Ar_radius = sqrt(At/pi)*100;
% mach value
Me = sqrt((2/(gam-1)) * (((Pc/Pa)^((gam-1)/gam))-1));
% exit area
Ae = (At/Me)*(((1+(((gam-1)/2)*Me*Me))/((gam+1)/2))^((gam+1)/(2*(gam-1))));
Ae_radius = sqrt(Ae/pi)*100;
% exit velocity
Ve = sqrt( (2*(gam/(gam-1))) * (R_un * (Tc/M)) * (1-((Pe/Pc)^((gam-1)/gam))));
% thrust
F = (m_dot*Ve) + ((Pe-Pa)*Ae)

% nozzle length
Ln = ((Ae_radius - Ar_radius)/sin(theta))*cos(theta);



figure;
hold on;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';  % setting y axis location to origin
title('Cone Nozzle Geometry in cm')
axis([-0.5*Ln 1.5*Ln -2*Ae_radius 2*Ae_radius ])
plot(zeros(100),linspace(-Ar_radius, Ar_radius))
plot(Ln+ zeros(100), linspace(-Ae_radius, Ae_radius))
slope = (Ae_radius - Ar_radius)/Ln;
plot(linspace(0,Ln) ,linspace(Ar_radius, Ae_radius))
plot(linspace(0,Ln) ,-linspace(Ar_radius, Ae_radius))
hold off;

% here I want to see how the thrust chamnges based on Pe shift
% figure;
% hold on;
% m_dot = 8;      % kg/s, mass flow rate
% gam = 1.2;      % specific heat ratio
% Tc = 3060;      % flame temperature, Kelvin
% Pc = 2758000;   % chamber pressure, MPa
% M = 22.9/1000;  % molecular weight km/mol
% Pa = 101325;    % atmospheric pressure
% Pe= linspace(0, 4*Pa);
% R_un = 8.314;   % universal gas constant
% Pa= linspace(101325, 0);
% Pt = Pc * ((1+((gam-1)/2))^((-gam)/(gam-1)));
% Tt = Tc * ((1+((gam-1)/2))^-1);
% At = (m_dot/Pt)*sqrt((R_un*Tt)/(M*gam));
% Me = sqrt((2/(gam-1)) .* (((Pc./Pa).^((gam-1)/gam))-1));
% Ae = (At./Me).*(((1+(((gam-1)/2).*Me.*Me))./((gam+1)/2)).^((gam+1)/(2*(gam-1))));
% Ve = sqrt((2*(gam/(gam-1))) * (R_un * (Tc/M)) .* (1-((Pe./Pc).^((gam-1)/gam))));
% F = (m_dot.*Ve) + ((Pe-Pa).*Ae);
% title('Thrust vs Exit Pressure');
% xlabel('Exit Pressure Pascals');
% ylabel('Total Thrust Newtons');
% plot(Pe, F);
% clear;
% 
% % here I want to see how the thrust chamnges based on Ambient pressure
% figure;
% hold on;
% m_dot = 8;      % kg/s, mass flow rate
% gam = 1.2;      % specific heat ratio
% Tc = 3060;      % flame temperature, Kelvin
% Pc = 2758000;   % chamber pressure, MPa
% M = 22.9/1000;  % molecular weight km/mol
% R_un = 8.314;   % universal gas constant
% Pe= 101325;
% Pa= linspace(0, 101325);
% Pt = Pc * ((1+((gam-1)/2))^((-gam)/(gam-1)));
% Tt = Tc * ((1+((gam-1)/2))^-1);
% At = (m_dot/Pt)*sqrt((R_un*Tt)/(M*gam));
% Me = sqrt((2/(gam-1)) .* (((Pc./Pa).^((gam-1)/gam))-1));
% Ae = (At./Me).*(((1+(((gam-1)/2).*Me.*Me))./((gam+1)/2)).^((gam+1)/(2*(gam-1))));
% Ve = sqrt( (2*(gam/(gam-1))) * (R_un * (Tc/M)) * (1-((Pe/Pc)^((gam-1)/gam))));
% F = (m_dot.*Ve) + ((Pe-Pa).*Ae);
% title('Thrust vs Ambient Pressure');
% xlabel('Ambient Pressure Pascals');
% ylabel('Total Thrust Newtons');
% plot(Pa, F);
% clear;
% 
% 
% % here I want to see how the thrust chamnges based on Ambient pressure
% figure;
% hold on;
% m_dot = linspace(0,16);      % kg/s, mass flow rate
% gam = 1.2;      % specific heat ratio
% Tc = 3060;      % flame temperature, Kelvin
% Pc = 2758000;   % chamber pressure, MPa
% M = 22.9/1000;  % molecular weight km/mol
% R_un = 8.314;   % universal gas constant
% Pe = 101325;
% Pa = Pe; 
% Pt = Pc * ((1+((gam-1)/2))^((-gam)/(gam-1)));
% Tt = Tc * ((1+((gam-1)/2))^-1);
% At = (m_dot/Pt)*sqrt((R_un*Tt)/(M*gam));
% Me = sqrt((2/(gam-1)) .* (((Pc./Pa).^((gam-1)/gam))-1));
% Ae = (At./Me).*(((1+(((gam-1)/2).*Me.*Me))./((gam+1)/2)).^((gam+1)/(2*(gam-1))));
% Ve = sqrt( (2*(gam/(gam-1))) * (R_un * (Tc/M)) * (1-((Pe/Pc)^((gam-1)/gam))));
% F = (m_dot.*Ve) + ((Pe-Pa).*Ae);
% title('Thrust vs M_dot');
% xlabel('M_dot, Kg/s');
% ylabel('Total Thrust Newtons');
% plot(m_dot, F);
% clear;





