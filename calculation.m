
u = 0.5;  % east-west component in m/s
v = 0.2;  % north-south component in m/s
w = -0.1; % vertical component in m/s
sref = 7.25;  % Reference area m^2
rho = 1.8; %kg/m^3 
b=3.07; %m
c=1.02; %m
h=0.14*c ; %m
rho1= 5.43 * exp(-0.512 * 10^-3*h);  % Air density at given altitude on Titan
alpha_max = 45*pi/180;      %rad

% Calculate the magnitude of the velocity vector
V = sqrt(u^2 + v^2 + w^2);

% Calculate dynamic pressure
qbar = 0.5 * rho * V^2;

% Set control surface deflection
ds = 1; %(FA_max)             %  ds = 0(For FA_min);
da = 1; %(MA_max)             %  da = 0(For MA_min);
p = 0.1; % Roll rate in rad/s
q = 0.1; % Pitch rate in rad/s
r = 0.1; % Yaw rate in rad/s

% Define aerodynamic coefficients
CL0 = 0.091;
CLalpha = 0.9;
CLds1 = 0.21;
CD0 = 0.25;
CDalpha2 = 0.12;
CDds = 0.3;
CSbeta = -0.23;
% Calculate angles alpha and beta
alpha = inv(w/u);
%beta = atan(v / sqrt(u^2 + w^2));
beta=0;

Cl_beta=-0.036; Cl_da=-0.0035; Cl_p=-0.84;Cl_r=-0.082;
Cm_0=0.25; Cm_alpha=-0.72; Cm_q=-1.49;
Cn_beta=-0.0015; Cn_da=0.0155; Cn_p=-0.082; Cn_r=-0.27;

% Calculate coefficients of Lift, Drag, and Side force
CL = CL0 + 0*CLalpha * alpha_max + CLds1 * ds;
CD = CD0 + 0*CDalpha2 * alpha_max^2 + CDds * ds;  
CS = CSbeta * 0;

Cl=Cl_beta *beta +Cl_da *da +Cl_p * (b/(2*V)) * p +Cl_r *(b/(2*V))*r;
Cm=Cm_0 +Cm_alpha *alpha +Cm_q *(c/(2*V))*q;
Cn=Cn_beta *beta +Cn_da *da +Cn_p *(b/(2*V)) *p + Cn_r *(b/(2*V))*r;

% Calculate aerodynamic forces
FA = qbar * sref * [-CD; -CS; -CL];
MA = qbar * sref * [Cl*b; Cm*c ; Cn*b];

% Calculate the norm of the aerodynamic force vector
normFA = norm(FA);
normMA = norm(MA);

% Display the aerodynamic force vector and its norm
disp('Aerodynamic Force Vector FA:');
disp(FA);
disp('Norm of the Aerodynamic Force Vector FA:');
disp(normFA);

% Display the aerodynamic force vector and its norm
disp('Aerodynamic Force Vector FA:');
disp(MA);
disp('Norm of the Aerodynamic Force Vector FA:');
disp(normMA);
