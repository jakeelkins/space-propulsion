% grad report code for AEM 529

clc
clear all

% take engine given parameters and calculate delta-v ability using ideal
% rocket

% ----this part for RD-180, at sea level, using assumptions about first
% stage ONLY-----
g0 = 9.81;
P_a = 101325;

m_dry = 21054; % kg
m_prop = 305143;

m_wet = m_prop + m_dry;

% Isp at sea level
I_sp_given = 311.3; % s

% thrust at SL
thrust = 3827e+3; % in Newtons

% calculate Mach exit, c_star, exit pressure, Isp ourselves and compare
% using fuel parameters and nozzle expansion value
exp_ratio = 36.87;

P_chamber = 26.7e+6; % Pa

%fuel params for RP-1, LOx at o/f 2.72
gamma_prop = 1.2175;
T_flame = 3600;
comb_eff = 1.05; % frozen flow assumed
nozzle_eff = 0.98;
molar_mass = 23.75; % kg/kmol

R = 8314/molar_mass;

%----------start----------
% iterate
mach_exit = 1.0;
arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));

while abs(exp_ratio - arg) > 0.001
    mach_exit = mach_exit + 0.00001;
    arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
end


% now get exit pressure for given chamber pressure
P_exit = P_chamber *(1 + ((gamma_prop-1)/2)*mach_exit^2)^(gamma_prop/(1-gamma_prop));

char_velo = (comb_eff*sqrt(gamma_prop*R*T_flame))/(gamma_prop*(2/(gamma_prop+1))^((gamma_prop+1)/(2*gamma_prop - 2)));

% and Isp using nozzle velo
I_sp = nozzle_eff*((char_velo*gamma_prop/g0)*sqrt((2/(gamma_prop-1))*((2/(gamma_prop+1))^((gamma_prop+1)/(gamma_prop-1)))*(1 - ((P_exit/P_chamber)^((gamma_prop-1)/gamma_prop)))) + (char_velo*exp_ratio*(P_exit-P_a)/(g0*P_chamber)));

fprintf('------for RD-180------\n')
fprintf('Isp given: %.2f s, Isp calc: %.2f s, delta: %.2f s\n', I_sp_given, I_sp, abs(I_sp - I_sp_given));

delta_v = I_sp*g0*log(m_wet/m_dry);

% equiv exit velo
v_exit_eq = I_sp*g0;

% mass flow rate of fuel and oxidizer assuming const thrust
m_dot = thrust/v_exit_eq;

fprintf('delta-v capability: %.2f m/s \n', delta_v);
fprintf('c star: %.2f m/s \n', char_velo);
fprintf('m_dot exit (const thrust): %.2f kg/s \n', m_dot);
fprintf('---------------------\n')







% ----this part for BE-4, at sea level, using assumptions about first
% stage ONLY-----
g0 = 9.81;
P_a = 101325;

m_dry = 21054; % kg
m_prop = 305143;

m_wet = m_prop + m_dry;

% thrust at SL
thrust = 2400e+3; % in Newtons

% calculate Mach exit, c_star, exit pressure, Isp ourselves and compare
% using fuel parameters and nozzle expansion value
exp_ratio = 31;

P_chamber = 13.4e+6; % Pa

%fuel params for LNG, LOx at o/f 3
gamma_prop = 1.23;
T_flame = 3000;
comb_eff = 1.05; % frozen flow assumed
nozzle_eff = 0.98;
molar_mass = 28; % kg/kmol

R = 8314/molar_mass;

%----------start----------
% iterate
mach_exit = 1.0;
arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));

while abs(exp_ratio - arg) > 0.001
    mach_exit = mach_exit + 0.00001;
    arg = (1/mach_exit)*((2/(gamma_prop+1))*(1 + ((gamma_prop-1)/2)*mach_exit^2))^((gamma_prop+1)/(2*gamma_prop - 2));
end


% now get exit pressure for given chamber pressure
P_exit = P_chamber *(1 + ((gamma_prop-1)/2)*mach_exit^2)^(gamma_prop/(1-gamma_prop));

char_velo = (comb_eff*sqrt(gamma_prop*R*T_flame))/(gamma_prop*(2/(gamma_prop+1))^((gamma_prop+1)/(2*gamma_prop - 2)));

% and Isp using nozzle velo
I_sp = nozzle_eff*((char_velo*gamma_prop/g0)*sqrt((2/(gamma_prop-1))*((2/(gamma_prop+1))^((gamma_prop+1)/(gamma_prop-1)))*(1 - ((P_exit/P_chamber)^((gamma_prop-1)/gamma_prop)))) + (char_velo*exp_ratio*(P_exit-P_a)/(g0*P_chamber)));

fprintf('------for BE-4------\n')
fprintf('Isp calc: %.2f s\n', I_sp);

delta_v = I_sp*g0*log(m_wet/m_dry);

% equiv exit velo
v_exit_eq = I_sp*g0;

% mass flow rate of fuel and oxidizer assuming const thrust
m_dot = thrust/v_exit_eq;

fprintf('delta-v capability: %.2f m/s \n', delta_v);
fprintf('c star: %.2f m/s \n', char_velo);
fprintf('m_dot exit (const thrust): %.2f kg/s \n', m_dot);
fprintf('---------------------\n')

