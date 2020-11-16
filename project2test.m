% Analysis of RL10 rocket engine bell nozzle with varying flow types

% assume: calorically perfect gas (constant ratio of specific heats)
% const molecular weight of gas
% no injected gases
% no internal drag

clear all

%polynomial fit for post-throat area:
%will be looped later
x = 0.0001;
radius = 0.1011*x^3 - 0.3837*x^2 + 0.6787*x + 0.0604;

%--------------------
gamma = 1.33;

%f = 0.005;
f = 0;
%assume dT_0 = some number(play with this number)
dT_0 = 0;

P_i0 = 50101325;
P_0 = P_i0;
T_i0 = 3000;
T_0 = T_i0;

%start analysis at choked flow at throat
M = 1;
m2 = M^2;

P = P_0*((1+((gamma-1)/2)*M^2)^((1-gamma)/gamma));
T = T_0*((1+((gamma-1)/2)*M^2)^(-1));

%---------------------
%throat
r_initial = 0.0604;
D_initial = 2*r_initial;
A_prev = pi*(D_initial^2)/4;
A_star = pi*(D_initial^2)/4;



T_list = [ ];
P_list = [ ];
M_list = [ ];
m2_list = [ ];

dx = 0.0001;

i = 0;
x = 0.0001;


i = i+1;

r = 0.1011*x^3 - 0.3837*x^2 + 0.6787*x + 0.0604;
D = 2*r;
A = pi*(D^2)/4;

T_0 = T_0 + dT_0;

%dA = A - A_prev;
dA = 2*pi*r*(0.3033*x^2 - 0.7674*x + 0.6787)*dx;

chi = (1 + ((gamma-1)/2)*M^2);

if M == 1
    Z = 1 - (1.01)^2;
end

dm2_over_m2 = ((-2*chi*dA)/(Z*A)) + (((1 + gamma*M^2)*chi*dT_0)/(Z*T_0)) + (((gamma*M^2)*chi*(4*f*dx))/(Z*D));

dm2 = dm2_over_m2*M^2;

%-----------now do dP/P-------------

dP_over_P = (((gamma*M^2)*dA)/(Z*A)) - ((gamma*M^2)*chi*dT_0/(Z*T_0)) - (((gamma*M^2)*(1 + (gamma - 1)*M^2)*4*f*dx)/(2*(Z)*D));

dP = dP_over_P*P;

%-----------dT/T------------------

dT_over_T = ((((gamma-1)*M^2)*dA)/(Z*A)) + (((1 - gamma*M^2)*chi*dT_0)/(Z*T_0)) - ((gamma*(gamma-1)*(M^4)*4*f*dx)/(2*(Z)*D));

dT = dT_over_T*T;