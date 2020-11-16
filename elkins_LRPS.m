% LRPS rocket sizing

clc
clear all

% ---------delta-v calcs for orbit raise----------

r1 = 42164; % km
r2 = 554200; % km

% consts
mu_e = 0.399e+6; %in kg, whatever the units are for this
g0 = 9.81;

% hohmann transfer used, one-way

v_1c = sqrt(mu_e/r1); % v in circ orbit at r1

a_transfer = 0.5*(r1 + r2); %SMA of transfer ellipse
e_transfer = 1 - (r1/a_transfer); % ecc. of transfer ellipse

v_transf_p = sqrt(2*(mu_e*((1/r1)-(1/(2*a_transfer))))); % v after first maneuver

dv1 = v_transf_p - v_1c; % delta-v of first maneuver

v_transf_a = sqrt(2*(mu_e*((1/r2)-(1/(2*a_transfer))))); % v at apogee of transfer ellipse

v_2c = sqrt(mu_e/r2); % v in circ orbit at r2

dv2 = v_2c - v_transf_a; % delta-v of second maneuver

dv_tot = dv1 + dv2; %total dv of maneuver

fprintf('total delta-v of orbital maneuver: %.4f km/s \n', dv_tot);

%----------------------------
%-----begin LRPS stuff-------
%----------------------------

m_payload = 50000; % kg



% CHOOSE: initial max mass value, dv margin
m_initial = 127000; %kg
thrust_to_weight = 0.2;
dv_margin = 0.1; % edit this number around
dv_design = (1 + dv_margin)*dv_tot*1000; %in meters from here on out
%------------------------------

thrust_req = thrust_to_weight*g0*m_initial;

% rough estimates on thrust chamber sizing:
engine_thrust_to_weight = (0.0006098*thrust_req + 13.44);
m_engine = thrust_req/(g0*(engine_thrust_to_weight));
length_engine = (0.0054*thrust_req + 31.92)/100;
diam_engine = (0.00357*thrust_req +14.48)/100;

% chosen: O/F of 5. params are read from data at O/F = 5.
% assume combustion efficiency of 1.
ox_to_fuel = 5;
T_flame = 3250; %K
gamma_prop = 1.21;
molar_mass_prop = 11.8; % kg/kmol

P_chamber = 4000000; % Pa, picked from plot of historical data

nozzle_eff = 0.98;
comb_eff = 1.0;

% calculating spec impulse
R = 8314/molar_mass_prop; % exhaust gas const

% characteristic velocity of exhaust
char_velo = (comb_eff*sqrt(gamma_prop*R*T_flame))/(gamma_prop*(2/(gamma_prop+1))^((gamma_prop+1)/(2*gamma_prop - 2)));

% gear these vars up for plotting
plot_exp_ratios = 66:1:260;
plot_I_sp = [ ];

% iteratively solve for exit mach number
for exp_ratio = 66:1:260
    mach_exit = 0.01;
    arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
    while abs(exp_ratio - arg) > 2
        mach_exit = mach_exit + 0.01;
        %exp_ratio = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
        arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
    end
    
    P_exit = P_chamber *(1 + ((gamma_prop-1)/2)*mach_exit^2)^(gamma_prop/(1-gamma_prop));
    
    I_sp = nozzle_eff*((char_velo*gamma_prop/g0)*sqrt((2/(gamma_prop-1))*((2/(gamma_prop+1))^((gamma_prop+1)/(gamma_prop-1)))*(1 - ((P_exit/P_chamber)^((gamma_prop-1)/gamma_prop)))) + (char_velo*exp_ratio*P_exit/(g0*P_chamber)));
    
    plot_I_sp = [plot_I_sp I_sp];
end


figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');
plot(plot_exp_ratios, plot_I_sp);
xlabel('Nozzle Expansion Ratio');
ylabel('I_sp (s)');
box(axes1,'on');
grid on

% selected ratio: 70

exp_ratio = 70;

% now solve for exit mach num to get I_sp
mach_exit = 0.01;
arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
while abs(exp_ratio - arg) > 2
    mach_exit = mach_exit + 0.01;
    %exp_ratio = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
    arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
end

P_exit = P_chamber *(1 + ((gamma_prop-1)/2)*mach_exit^2)^(gamma_prop/(1-gamma_prop));

I_sp = nozzle_eff*((char_velo*gamma_prop/g0)*sqrt((2/(gamma_prop-1))*((2/(gamma_prop+1))^((gamma_prop+1)/(gamma_prop-1)))*(1 - ((P_exit/P_chamber)^((gamma_prop-1)/gamma_prop)))) + (char_velo*exp_ratio*P_exit/(g0*P_chamber)));


fprintf('for expansion ratio of %.1f (selected), I_sp = %.2f, exit pressure = %.2f , exit mach number = %.2f\n', exp_ratio, I_sp, P_exit, mach_exit);


% CHOOSE: f_inert (inert mass fraction)
f_inert = 0.44;

m_prop = (m_payload*(exp(dv_design/(I_sp*g0)) - 1)*(1 - f_inert))/(1 - f_inert*exp(dv_design/(I_sp*g0)));

m_inert = (f_inert/(1 - f_inert))*m_prop;

total_impulse = I_sp*g0*m_prop;

t_burn = total_impulse/thrust_req;

m_dot = thrust_req/(g0*I_sp);

m_dot_fuel = m_dot/(ox_to_fuel+1);

m_dot_ox = m_dot - m_dot_fuel;

% ---------thrust chamber round 2-----------
% CHOOSE: comb. chamber mach number, characteristic length
%ASSUME: mult factor of 3, pressure factor of 2, columbium for materiat
%(yield strength of 310 MPa), density 8500 kg/m3
% combustion chamber turn angle 45 deg.

mach_chamber = 0.2;
L_star = 0.9;
wall_density = 8500;
cc_angle = 45*pi/180;

A_throat = m_dot*char_velo/P_chamber;
throat_diam = sqrt(4*A_throat/pi);

A_exit = exp_ratio*A_throat;

chamb_to_throat_area_ratio = 8*((throat_diam*100)^(-0.6)) + 1.25;
A_chamber = chamb_to_throat_area_ratio*A_throat;

chamber_diam = sqrt(4*A_chamber/pi);

% nozzle contraction ratio
nozzle_contr_ratio = (1/mach_chamber)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_chamber^2))^((gamma_prop+1)/(2*gamma_prop - 2));

length_chamber = L_star*(1/nozzle_contr_ratio);

nozzle_exit_diam = sqrt(4*exp_ratio*A_throat/pi);

chamber_thickness = P_chamber*chamber_diam*3/(2*310000000);

%--comb. chamber mass--
m_cc = pi*wall_density*chamber_thickness*(2*(chamber_diam/2)*length_chamber + ((((chamber_diam/2)^2) - ((throat_diam/2)^2))/(tan(cc_angle))));



% --- nozzle length, angles---
% ASSUME: go from concial -> bell frome graphs
bell_nozzle_frac = 0.675;

conic_half_angle = acos(2*nozzle_eff - 1);

conic_nozzle_length = (nozzle_exit_diam - throat_diam)/(2*tan(conic_half_angle));

bell_nozzle_length = bell_nozzle_frac*conic_nozzle_length;

% ---- nozzle mass ----
%ASSUME: throat wall thickness ~ 0.5 chamber thickness, tapers to zero
%thickness
throat_thickness = 0.5*chamber_thickness;

f1 = -throat_thickness/bell_nozzle_length;
f2 = ((nozzle_exit_diam/2)-(throat_diam/2))/bell_nozzle_length;

m_nozzle = 2*pi*wall_density*bell_nozzle_length*((1/3)*f1*f2*(bell_nozzle_length^2) + 0.5*(f1*(throat_diam/2) + f2*throat_thickness)*bell_nozzle_length + (throat_diam/2)*throat_thickness);

% engine masses

m_engine = m_nozzle/(0.101);
m_feed = 0.062*m_engine;

%----------------tank sizing---------------
density_lox = 1142;
density_lh2 = 71;

m_ox = m_dot_ox*t_burn;
m_fuel = m_dot_fuel*t_burn;

V_fuel = m_fuel/density_lh2;
V_ox = m_fuel/density_lox;

%get tank pressures from curve fit to volume
P_fuel_tank = (10^(-0.1068*log10(V_fuel) - 0.2588))*10^6;
P_ox_tank = (10^(-0.1068*log10(V_ox) - 0.2588))*10^6;

%ASSUME: metallic tank so factor is 2500 below, burst factor of 2
m_fuel_tank = 2*P_fuel_tank*V_fuel/(2500*g0);
m_ox_tank = 2*P_ox_tank*V_ox/(2500*g0);


% -----pressurant tank: helium-----
gamma_He = 1.66;
molar_mass_He = 4.003;

%ASSUME: init temp 273K, init pressure 21 MPa, avg pressure 938,000
pressurant_Ti = 273;
pressurant_Pf = 938000;
pressurant_Pi = 21000000;

pressurant_Tf = pressurant_Ti*((pressurant_Pf/pressurant_Pi)^((gamma_He-1)/gamma_He));

% mass of actual helium mass used as pressurant (ideal gas w 5% margin)
m_pressurant = (1.05*pressurant_Pf*(V_ox + V_fuel)*molar_mass_He)/(8314*pressurant_Tf);

% now we need the tank mass (ideal gas), initial guess
V_pressurant_tank = m_pressurant*8314*pressurant_Ti/(pressurant_Pi*molar_mass_He);

% solution part?
m_press_tank = 1.05*pressurant_Pi*V_pressurant_tank/(g0*6350);


% strcture support sys: assume 10% inert increase
m_structure = 0.1*m_inert;

m_total = m_payload + m_fuel + m_ox + m_pressurant + m_feed+m_nozzle+m_cc + m_fuel_tank + m_ox_tank + m_press_tank + m_feed + m_structure;

% summarize results nice and pretty
fprintf('+------------------------------------+\n')
fprintf('|     initial mass: %.2f        |\n', m_initial)
fprintf('|     payload mass: %.2f         |\n', m_payload)
fprintf('|     thrust: %.2f              |\n', thrust_req)
fprintf('|     thrust duration: %.2f        |\n', t_burn)
fprintf('|     specific impulse: %.2f       |\n', I_sp)
fprintf('|     total impulse: %.2f    |\n', total_impulse)
fprintf('|     fuel flow rate: %.2f           |\n', m_dot_fuel)
fprintf('|     ox flow rate: %.2f            |\n', m_dot_ox)
fprintf('|     fuel mass: %.2f             |\n', m_fuel)
fprintf('|     ox mass: %.2f              |\n', m_ox)
fprintf('|     pressurant mass: %.2f        |\n', m_pressurant)
fprintf('|     thrust chamber mass: %.2f    |\n', m_feed+m_nozzle+m_cc)
fprintf('|     h2 tank mass: %.2f          |\n', m_fuel_tank)
fprintf('|     ox tank mass: %.2f           |\n', m_ox_tank)
fprintf('|     pressurant tank mass: %.2f  |\n', m_press_tank)
fprintf('|     feed system mass: %.2f        |\n', m_feed)
fprintf('|     support struct mass: %.2f   |\n', m_structure)
fprintf('|     total mass: %.2f          |\n', m_total)
fprintf('|     throat area: %.2f              |\n', A_throat)
fprintf('|     throat diam: %.2f              |\n', throat_diam)
fprintf('|     exit area: %.2f                |\n', A_exit)
fprintf('|     exit diam: %.2f                |\n', nozzle_exit_diam)
fprintf('|     chamber length: %.2f           |\n', length_chamber)
fprintf('|     chamber area: %.2f             |\n', A_chamber)
fprintf('|     chamber diam: %.2f             |\n', chamber_diam)
fprintf('|     wall thickness: %.2f           |\n', chamber_thickness)
fprintf('|     bell nozzle length: %.2f       |\n', bell_nozzle_length)
fprintf('+------------------------------------+\n')
