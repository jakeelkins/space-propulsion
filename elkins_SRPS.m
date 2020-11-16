% SRPS rocket design

clc
clear all

% ---------delta-v calcs for orbit raise----------

r1 = 6374 + 200; % km
r2 = 6374 + 500; % km

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
%-----begin SRPS stuff-------
%----------------------------

m_payload = 5000; % kg

% ASSUMPTIONS
I_sp = 290;   %s

%---performance---
P_c_avg = 5.17e+6;    %Pa, avg chamber pressure
t_burn = 80;    %s

%---motor case---
L_to_D_ratio = 2;
volumetric_loading = 0.92;
P_c_max_to_avg = 1.4;
safety_factor = 1.25;
F_tu = 1.34e+9; %Pa
density_case = 1550;

%---propellant---
density_prop = 1800;
char_velo = 1527;   %m/s

%---insulation---
density_insul = 1100;

%---nozzle---
% effective divergence angle
theta_conical = 14*pi/180;


%-------------begin calcs-------------
mass_ratio = exp((dv_tot*1000)/(g0*I_sp));

% ASSUME: mass fract of 0.9 from graphs
f_prop = 0.9;

% propellant mass
m_prop = m_payload*((f_prop*(mass_ratio - 1))/(1 - mass_ratio*(1 - f_prop)));

% some efficiency
eta_D = 0.92;

% case volume
V_case = m_prop/(eta_D*density_prop);

% case diameter
D_case = (V_case/((pi/6) + (pi/4)*(L_to_D_ratio-1)))^(1/3);

% entire case
L_case = L_to_D_ratio*D_case;

% length of cylinder part
L_cy = L_case - D_case;

% burst pressure
%ASSUME: 
P_b = P_c_avg*(P_c_max_to_avg)*safety_factor;

% case thickness
t_case = P_b*D_case/(2*F_tu);

% mass of pressure vessel:
m_pv = density_case*t_case*(D_case^2)*pi*(1 + (L_cy/D_case));

% mass of skirt
m_sk = density_case*t_case*pi*D_case^2;

% update motor case mass by 10% (not reported in proj)
m_case = 1.1*(m_pv + m_sk);

% area of exposed wall for insulation
A_w = pi*D_case^2 + pi*D_case*L_cy;

%--mass of insulation estimate--
% ASSUME: 0.05 m nozzle submergence length
X = 0.05;

% length of submergence
L_sub = 100*X/L_case;

% insulation mass
m_insul = (1.788e-9)*(m_prop^-1.33)*(t_burn^0.965)*(L_to_D_ratio^0.144)*(L_sub^0.058)*((A_w*100*100)^2.69);

% insulation thickness
t_insul = m_insul/(density_insul*A_w);



%-------nozzle, igniter--------
% throat diameter
D_throat = sqrt(4*char_velo*m_prop/(pi*t_burn*P_c_avg));

% need exit diam < case diam for packaging
% ASSUME: exit diam 80% case diam
exit_diam_ratio = 0.8;

% diam of nozzle exit
D_exit = exit_diam_ratio*D_case;

% area expansion ratio
exp_ratio = ((D_exit)^2)/(D_throat^2);

% length of nozzle from throat to exit
L_noz = (D_exit - D_throat)/(2*tan(theta_conical));

% mass of nozzle
m_noz = (0.256e-4)*(((((m_prop*char_velo)^1.2)*(exp_ratio^0.3))/(((P_c_avg*10^-6)^0.8)*(t_burn^0.6)*(tan(theta_conical))^0.4))^0.917);

% mass of nozzle system (igniter, TVC system)
% ASSUME: historical factor of 1.5
hist_factor = 1.5;
m_noz_sys = hist_factor*m_noz;


% ---ballistics---
%ASSUME: web fraction of 0.9
web_frac = 0.9;

% outer radius (radius of case minues thicknesses)
r_outer = (D_case/2) - t_case - t_insul;

% inner radius of bore (use web frac here)
r_inner = (D_case/2)*(1 - web_frac);

% burn area
%A_bo = 2*pi*((r_outer^2)-(r_inner^2));
A_bo = 2*pi*((r_outer^2)-(r_inner^2)) + L_cy*(r_outer-r_inner);% - 0.05*(r_outer-r_inner)+ 2*pi*r_outer - (r_outer*r_inner)

% do axis of rev to get volume of burn
V_bo = pi*(A_bo/2)^2;

% burning rate (in cm)
r_b = 0.316*(P_c_avg*10^-6)^0.3;

%---do updates---
m_inert = m_case + m_insul + m_noz_sys;

% update m_prop based on volume of burn
m_prop_new = V_bo*density_prop;

f_prop_new = m_prop_new/(m_prop_new + m_inert);

% new mass ratio R for rocket eq dv calc: (mi/mf)
new_mass_ratio = (m_payload + m_prop_new + m_inert)/(m_payload + m_inert);

% get actual dv for this design
new_dv = g0*I_sp*log(new_mass_ratio);

if new_dv >= dv_tot*1000
    fprintf('design success! dv margin of %.2f m/s \n', new_dv - dv_tot*1000);
end


%----------------------make prints nice and pretty-------------------
fprintf('+-----------------------------------+\n')
fprintf('|     calc delta-v: %.2f          |\n', dv_tot*1000)
fprintf('|     mass ratio: %.2f              |\n', mass_ratio)
fprintf('|     m_payload: %.2f            |\n', m_payload)
fprintf('|     case volume: %.2f             |\n', V_case)
fprintf('|     case length: %.2f             |\n', L_case)
fprintf('|     case diam: %.2f               |\n', D_case)
fprintf('|     cylinder length: %.2f         |\n', L_cy)
fprintf('|     burst pressure: %.2f    |\n', P_b)
fprintf('|     case thickness: %.5f       |\n', t_case)
fprintf('|     pressure vessel mass: %.2f    |\n', m_pv)
fprintf('|     skirt mass: %.2f              |\n', m_sk)
fprintf('|     wall area: %.2f               |\n', A_w)
fprintf('|     nozzle sub length (X): %.2f   |\n', X)
fprintf('|     insulation mass: %.2f        |\n', m_insul)
fprintf('|     insulation thickness: %.4f  |\n', t_insul)
fprintf('|     throat diam: %.2f             |\n', D_throat)
fprintf('|     exit diam: %.2f               |\n', D_exit)
fprintf('|     nozzle length: %.2f           |\n', L_noz)
fprintf('|     nozzle mass: %.2f             |\n', m_noz)
fprintf('|     nozzle system mass: %.2f     |\n', m_noz_sys)
fprintf('|     inner radius: %.2f            |\n', r_inner)
fprintf('|     outer radius: %.2f            |\n', r_outer)
fprintf('|     burn area: %.2f               |\n', A_bo)
fprintf('|     burn volume: %.2f             |\n', V_bo)
fprintf('|     new propellant mass: %.2f   |\n', m_prop_new)
fprintf('|     burn rate: %.2f               |\n', r_b)
fprintf('|     inert mass: %.2f             |\n', m_inert)
fprintf('|     new mass frac: %.2f           |\n', f_prop_new)
fprintf('|     new mass ratio: %.2f          |\n', new_mass_ratio)
fprintf('|     new delta-v: %.2f           |\n', new_dv)
fprintf('+-----------------------------------+\n')


