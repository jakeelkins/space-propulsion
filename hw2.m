clc
clear all
format long


% initial conditions
M = 2.016;    %kg/kmol
D = 0.01;    %m

T01 = 500;  %K
P1 = 300000;    %Pa

Pa = 200000;    %Pa

Qdot = 100000;  %W

%-------------------------
Apipe = 0.25*pi*(D^2);

Rbar = 8314.3; %J/kmolK

R = Rbar/M;

Q_guess = 1;

%guess T1 iteratively
T1 = 490; %K

while abs(Qdot - Q_guess) > 0.1

    T1 = T1 + 0.1;
    fprintf('curr T1 = %d\n\n\n',T1)

    tspan0 = [298.16 T01];
    [t, h01_int] = ode45(@calc_cp, tspan0, 0);

    h01_int = h01_int(length(h01_int)); 

    tspan = [298.16 T1];
    [t, h1_int] = ode45(@calc_cp, tspan, 0);

    h1_int = h1_int(length(h1_int));

    v1 = 2*sqrt(h01_int-h1_int);

    G = P1*v1/(R*T1);

    for T2 = 1500:0.1:2000
        rootA = G;
        rootB = -(G*v1 + P1);
        rootC = G*R*T2;

        v2_vec = roots([rootA rootB rootC]);

        v2 = v2_vec(isreal(v2_vec));
        v2 = v2(v2>0);

        if length(v2) < 1
            continue
        end

        if length(v2) > 1
            v2 = max(v2);
        end

        P2 = G*R*T2/v2;

        tspan2 = [298.16 T2];
        [t, h2_int] = ode45(@calc_cp, tspan2, 0);

        h2_int = h2_int(length(h2_int));

        Q_guess = (G*Apipe)*(h2_int + 0.5*v2^2 - h01_int);
        dq = (Q_guess - Qdot);
        abs_dq = abs(dq);
        fprintf('curr dQ = %d, abs: %d\n',dq, abs_dq);

        if abs(dq) < 0.5  
            flag = 1;
            break;
        end

    end
    
    if flag==1
        break
    end


end

fprintf('DONE. results:\n\nG = %d, T1 = %d, T2 = %d, v1 = %d, v2 = %d\n',G, T1, T2, v1, v2);







function cp = calc_cp(t, y)
    cp = (1000/2.016)*(56.505 - 702.79*(t/100)^(-0.75) + 1165*(t/100)^-1 - 560.7*(t/100)^-1.5);
end
