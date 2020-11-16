% cold gas thruster design, calculating time it will take to empty tank

%-----constants------
%propellant: nitrogen
R = 297;    %J/kgK
gamma = 1.4;

%--------initial conditions--------
mass_expelled = 150; %kg

P_0i = 40000000;    %Pa
T_0i = 298.15;    %K
thrust_eff = 0.98;
P_atm = 997.8;   %Pa
nozzle_throat_diam = 0.005;  %m
nozzle_exit_diam = 0.025;   %m
P_c = 1290000;  %MPa. total pressure after valve
tank_vol = 0.3735;   %cubic meters

%-----------------------------------
c_p = gamma*R/(gamma-1);
A_t = 0.25*pi*nozzle_throat_diam^2;

%-- init for iteration
T_0 = T_0i;
P_0 = P_0i;
m_thrusted = 0;
time_taken = 0;
%m_dot = ((A_t*P_c)/sqrt(gamma*R*T_0i))*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));

P0_list = [ ];
T0_list = [ ];
ue_list = [ ];
mdot_list = [ ];
thrust_list = [ ];

dt = 1;  %s

while m_thrusted < mass_expelled
%need variation in P_c, T_0 ,--------assume P_c const?

const = (A_t/tank_vol)*sqrt(R/gamma)*P_c*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))))*((T_0i^(gamma/(gamma-1)))/P_0i);

dT_0 = -const*(gamma-1)*dt*T_0^((gamma-3)/(2*(gamma-1)));
T_0 = T_0 + dT_0;

%which means.. (from 508 textbook)
dP_0 = (gamma/(gamma-1))*(dT_0/T_0)*P_0;
P_0 = P_0 + dP_0;

m_dot = ((A_t*P_c)/sqrt(gamma*R*T_0))*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));

dm = m_dot*dt;

m_thrusted = m_thrusted + dm;

time_taken = time_taken + dt;

%need exit pressure, which is assumed to be atmospheric pressure
P_e = P_atm;

exit_velo = sqrt(2*c_p*T_0*(1-(P_e/P_0)^((gamma-1)/gamma)));

thrust = thrust_eff*m_dot*exit_velo;

P0_list = [P0_list P_0];
T0_list = [T0_list T_0];
ue_list = [ue_list exit_velo];
mdot_list = [mdot_list m_dot];
thrust_list = [thrust_list thrust];

end

fprintf('time to expel 150 kg of propellant w/ varying T0: %.2f seconds\n',time_taken);

%-----find T_0 at end of burning----
%tspan = [0 5];
%y0 = 0;
%[t,y] = ode45(@(t,y) 2*t, tspan, y0);

%T_0_end = (((gamma+1)/(2*(gamma-1)))*(-const*(gamma-1))*time + T_0i^((gamma+1)/(2*(gamma+1))))^(2*(gamma-1)/(gamma+1));

x_list = 0.01:1:time_taken;

figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');
plot(x_list, P0_list);
xlabel('Time (s)');
ylabel('Tank Total Pressure (Pa)');
box(axes1,'on');
grid on

figure2 = figure;
axes2 = axes('Parent', figure2);
hold(axes2, 'on');
plot(x_list, T0_list);
xlabel('Time (s)');
ylabel('Tank Total Temperature (K)');
box(axes2,'on');
grid on


figure3 = figure;
axes3 = axes('Parent', figure3);
hold(axes3, 'on');
plot(x_list, ue_list);
xlabel('Time (s)');
ylabel('Nozzle Exit Velocity (m/s)');
box(axes3,'on');
grid on

figure4 = figure;
axes4 = axes('Parent', figure4);
hold(axes4, 'on');
plot(x_list, mdot_list);
xlabel('Time (s)');
ylabel('Mass Flow Rate Thru Nozzle (kg/s)');
box(axes4,'on');
grid on

figure5 = figure;
axes5 = axes('Parent', figure5);
hold(axes5, 'on');
plot(x_list, thrust_list);
xlabel('Time (s)');
ylabel('Thrust (N)');
ylim([0 50]);
box(axes5,'on');
grid on


%-----------------Part B-----------------
m_dot = ((A_t*P_c)/sqrt(gamma*R*T_0i))*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));
t_duration = mass_expelled/m_dot;

fprintf('time to expel 150 kg of propellant w/ const T0: %.2f seconds\n',t_duration);

%-- init for iteration
T_0 = T_0i;
P_0 = P_0i;
m_thrusted = 0;
time_taken = 0;
%m_dot = ((A_t*P_c)/sqrt(gamma*R*T_0i))*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));

P0_list = [ ];
T0_list = [ ];
ue_list = [ ];
mdot_list = [ ];
thrust_list = [ ];

dt = 1;  %s

while m_thrusted < mass_expelled
%need variation in P_c, T_0 ,--------assume P_c const?

const = (A_t/tank_vol)*sqrt(R/gamma)*P_c*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))))*((T_0i^(gamma/(gamma-1)))/P_0i);

dT_0 = 0;
T_0 = T_0 + dT_0;

%which means.. (from 508 textbook)
dP_0 = (gamma/(gamma-1))*(dT_0/T_0)*P_0;
P_0 = P_0 + dP_0;

m_dot = ((A_t*P_c)/sqrt(gamma*R*T_0))*(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));

dm = m_dot*dt;

m_thrusted = m_thrusted + dm;

time_taken = time_taken + dt;

%need exit pressure, which is assumed to be atmospheric pressure
P_e = P_atm;

exit_velo = sqrt(2*c_p*T_0*(1-(P_e/P_0)^((gamma-1)/gamma)));

thrust = thrust_eff*m_dot*exit_velo;

P0_list = [P0_list P_0];
T0_list = [T0_list T_0];
ue_list = [ue_list exit_velo];
mdot_list = [mdot_list m_dot];
thrust_list = [thrust_list thrust];

end


x_list = 0.01:1:time_taken;

figure6 = figure;
axes6 = axes('Parent', figure6);
hold(axes6, 'on');
plot(x_list, P0_list);
xlabel('Time (s)');
ylabel('Tank Total Pressure (Pa)');
box(axes6,'on');
grid on

figure7 = figure;
axes7 = axes('Parent', figure7);
hold(axes7, 'on');
plot(x_list, T0_list);
xlabel('Time (s)');
ylabel('Tank Total Temperature (K)');
box(axes7,'on');
grid on


figure8 = figure;
axes8 = axes('Parent', figure8);
hold(axes8, 'on');
plot(x_list, ue_list);
xlabel('Time (s)');
ylabel('Nozzle Exit Velocity (m/s)');
box(axes3,'on');
grid on

figure9 = figure;
axes9 = axes('Parent', figure9);
hold(axes9, 'on');
plot(x_list, mdot_list);
xlabel('Time (s)');
ylabel('Mass Flow Rate Thru Nozzle (kg/s)');
box(axes9,'on');
grid on

figure10 = figure;
axes10 = axes('Parent', figure10);
hold(axes10, 'on');
plot(x_list, thrust_list);
xlabel('Time (s)');
ylabel('Thrust (N)');
ylim([0 50]);
box(axes10,'on');
grid on




