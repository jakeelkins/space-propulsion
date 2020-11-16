% solving rocket diff eqs
clear all
clc

m_payload = 54;
m_prop = 316;
m_no_payload = 536;

thrust = 32068.5;

t_burnout = 30;

outer_diam = 0.352;

height = 4.95;

C_d = 0.2;

air_density = 1.225;

m_dot = m_prop/t_burnout;

A_f = pi*(outer_diam^2)/4;

g = 9.81;

%------------------------------------
% Part A, rocket launched veritcally
timestep = 1/240;
% ICs
time = 0;
down = 0;
u = 0;
y = 0;
mass = m_no_payload + m_payload;
dry_mass = mass - m_prop;

drag = 0.5*air_density*(u^2)*A_f*C_d;

du_dt = (thrust/mass) - (drag/mass) - g;

u_list = [u ];
y_list = [y ];
a_list = [du_dt ];
t_list = [time ];

while y >= 0

    mass = mass - m_dot*timestep;

    if abs(time-t_burnout) < 1e-6
        thrust = 0;
        m_dot = 0;
        fprintf('at burnout, mass is %.2f kg, %.2f kg expected.\n',mass, dry_mass);
    end
    
    if (abs(u) < 0.1) && (y > 0)
        down = 1;
        fprintf('time to peak: %.4f s\n', time);
    end
    
    drag = 0.5*air_density*(u^2)*A_f*C_d;
    
    if down == 1
        drag = -drag;
    end
    
    du_dt = (thrust/mass) - (drag/mass) - g;

    u = u + (du_dt)*timestep;
    y = y + (u*timestep);

    time = time + timestep;
    
    %append result
    u_list = [u_list u];
    y_list = [y_list y];
    a_list = [a_list du_dt];
    t_list = [t_list time];

end

fprintf('time to fall to Earth: %.4f s', time);


figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');
plot(t_list, u_list);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
box(axes1,'on');
grid on

figure2 = figure;
axes2 = axes('Parent', figure2);
hold(axes2, 'on');
plot(t_list, y_list);
xlabel('Time (s)');
ylabel('Altitude (m)');
box(axes2,'on');
grid on


figure3 = figure;
axes3 = axes('Parent', figure3);
hold(axes3, 'on');
plot(t_list, a_list);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
box(axes3,'on');
grid on
