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
% Part B, rocket launched at an angle of 45 deg
timestep = 1/240;
% ICs
time = 0;
const_phi = pi/4;
phi = const_phi;
down = 0;
burnout = 0;
ux = 0;
uy = 0;
x = 0;
y = 0;
mass = m_no_payload + m_payload;
dry_mass = mass - m_prop;

phi_list = [phi ];
ux_list = [ux ];
uy_list = [uy ];
x_list = [x ];
y_list = [y ];
t_list = [time ];

while y >= 0

    mass = mass - m_dot*timestep;
    
    u = sqrt(ux^2 + uy^2);

    if abs(time-t_burnout) < 1e-6
        thrust = 0;
        m_dot = 0;
        burnout = 1;
        fprintf('at burnout, mass is %.2f kg, %.2f kg expected.\n',mass, dry_mass);
    end
    
    if (abs(uy) < 0.1) && (y > 0)
        down = 1;
        phi = pi/2;
        fprintf('time to peak: %.4f s \n', time);
    end
    
    dragx = 0.5*air_density*(u^2)*A_f*C_d;
    dragy = 0.5*air_density*(u^2)*A_f*C_d;
    
    if down == 1
        dragy = -dragy;
    end
    
    duy_dt = ((thrust - dragy)/mass)*cos(phi) - g;
    dux_dt = ((thrust - dragx)/mass)*sin(phi);

    uy = uy + (duy_dt)*timestep;
    ux = ux + (dux_dt)*timestep;
    x = x + (ux*timestep);
    y = y + (uy*timestep);
    
    if burnout == 1
        phi = atan2(ux, uy);
    else
        phi = const_phi;
    end

    time = time + timestep;
    
    %append result
    phi_list = [phi_list phi];
    ux_list = [ux_list ux];
    uy_list = [uy_list uy];
    x_list = [x_list x];
    y_list = [y_list y];
    t_list = [t_list time];

end

x_to_mi = x/1609;
fprintf('range: %.4f m, approx. %.2f miles\n', x, x_to_mi);
fprintf('time of flight: %.4f s \n', time);


figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');
plot(t_list, ux_list);
xlabel('Time (s)');
ylabel('Velocity X (m/s)');
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
plot(t_list, phi_list);
xlabel('Time (s)');
ylabel('Flight Angle (rad.)');
box(axes3,'on');
grid on

figure4 = figure;
axes4 = axes('Parent', figure4);
hold(axes4, 'on');
plot(t_list, uy_list);
xlabel('Time (s)');
ylabel('Velocity Y (m/s)');
box(axes4,'on');
grid on


figure5 = figure;
axes5 = axes('Parent', figure5);
hold(axes5, 'on');
plot(t_list, x_list);
xlabel('Time (s)');
ylabel('Range (m)');
box(axes5,'on');
grid on



%---now we optimize to find initial flight angle that gives best---
ranges = [ ];
%------------------------------------
% go from 10 deg to 70 deg
for init_phi=0.087:0.001:(70*pi/180)
    timestep = 1/80;
    % ICs
    thrust = 32068.5;
    m_dot = m_prop/t_burnout;
    time = 0;
    phi = init_phi;
    down = 0;
    burnout = 0;
    ux = 0;
    uy = 0;
    duy_dt = 0;
    dux_dt = 0;
    x = 0;
    y = 0;
    mass = m_no_payload + m_payload;
    dry_mass = mass - m_prop;

    while y >= 0

        mass = mass - m_dot*timestep;

        u = sqrt(ux^2 + uy^2);

        if abs(time-t_burnout) < 1e-6
            thrust = 0;
            m_dot = 0;
            burnout = 1;
            %fprintf('at burnout, mass is %.2f kg, %.2f kg expected.\n',mass, dry_mass);
        end

        if (abs(uy) < 0.1) && (y > 0)
            down = 1;
            phi = pi/2;
            %fprintf('time to peak: %.4f s \n', time);
        end

        dragx = 0.5*air_density*(u^2)*A_f*C_d;
        dragy = 0.5*air_density*(u^2)*A_f*C_d;

        if down == 1
            dragy = -dragy;
        end

        duy_dt = ((thrust - dragy)/mass)*cos(phi) - g;
        dux_dt = ((thrust - dragx)/mass)*sin(phi);

        uy = uy + (duy_dt)*timestep;
        ux = ux + (dux_dt)*timestep;
        x = x + (ux*timestep);
        y = y + (uy*timestep);

        if burnout == 1
            phi = atan2(ux, uy);
        else
            phi = init_phi;
        end

        time = time + timestep;

    end
    ranges = [ranges x];
    %if abs(x) < 0.1
     %   fprintf('ux is %.4f, uy is %.4f\n', ux, uy);
    %end
end

phi_trials = 0.087:0.001:(70*pi/180);

[max_range, max_idx] = max(ranges);
max_phi = phi_trials(max_idx);
max_range_mi = max_range/1609;
max_phi_deg = max_phi*180/pi;
fprintf('Max range is %.4f km = %.2f mi, at an angle of %.4f rad. = %.2f deg.\n', max_range, max_range_mi, max_phi, max_phi_deg);

figure6 = figure;
axes6 = axes('Parent', figure6);
hold(axes6, 'on');
plot(phi_trials, ranges);
xlabel('Phi (rad.)');
ylabel('Range (m)');
box(axes6,'on');
grid on
