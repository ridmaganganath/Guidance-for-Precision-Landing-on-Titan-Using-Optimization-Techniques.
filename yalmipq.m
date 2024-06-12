

% Ensure YALMIP is added to your MATLAB path
addpath(genpath('D:\Softwares\YALMIP\YALMIP-master')); % Update with the actual path to YALMIP
yalmip('clear');

 %Vehicle fixed parameters
m_dry = 1505;       % Vehicle mass without fuel [kg]
Isp = 225;          % Specific impulse [s]
g0 = 9.80665;       % Standard earth gravity [m/s^2]
g = [0 ; -3.7114];  % Mars gravity vector [m/s^2]
max_throttle = 0.8; % Max open throttle [.%]
min_throttle = 0.3; % Min open throttle [.%]
T_max = 6 * 3100;   % Max total thrust force at 1.0 throttle [N]
phi = 27;  

% Calculate alpha
alpha = 1 / (Isp * g0 * cosd(phi));
r1 = min_throttle * T_max * cosd(phi);
r2 = max_throttle * T_max * cosd(phi);

% Initial conditions
m_wet = 1905;           % Vehicle mass with fuel [kg]
r0 = [  2 ; 1.5] * 1e3; % Initial position [x;z] [m]
v0 = [100 ; -75];       % Initial velocity [x;z] [m/s]

% Target conditions
rf = [ 0 ; 0 ];
vf = [ 0 ; 0 ];

tf = 75;  % Target end time [s]
dt = 1.0; % Discrete node time interval [s]
N = (tf / dt) + 1;
tv = 0:dt:tf;


% Now define your variables using YALMIP
r = sdpvar(2, N, 'full');
v = sdpvar(2, N, 'full');
u = sdpvar(2, N, 'full');
z = sdpvar(1, N, 'full');
s = sdpvar(1, N, 'full');



% Objective: Maximize ln of final mass -> Minimize fuel used
Objective = z(N);

% Constraints
Constraints = [r(:,1) == r0, v(:,1) == v0, z(1) == log(m_wet), ...
               r(:,N) == rf, v(:,N) == vf];

% Dynamical constraints
for i = 1:N-1
    Constraints = [Constraints, ...
        v(:,i+1) == v(:,i) + dt*g + (dt/2)*(u(:,i) + u(:,i+1)), ...
        r(:,i+1) == r(:,i) + (dt/2)*(v(:,i) + v(:,i+1)) + (dt^2/12)*(u(:,i+1) - u(:,i)), ...
        z(i+1) == z(i) - (alpha*dt/2)*(s(i) + s(i+1))];
end

% Thrust limit, mass flow limit, and extremal mass limits
for i = 1:N
    z0_term = m_wet - alpha * r2 * (i-1) * dt;
    z1_term = m_wet - alpha * r1 * (i-1) * dt;
    z0 = log(z0_term);
    z1 = log(z1_term);
    mu_1 = r1 / z0_term;
    mu_2 = r2 / z0_term;
    Constraints = [Constraints, ...
        norm(u(:,i)) <= s(i), ...
        s(i) >= mu_1 * (1 - (z(i) - z0) + (1/2)*(z(i) - z0)^2), ...
        s(i) <= mu_2 * (1 - (z(i) - z0)), ...
        z(i) >= z0, z(i) <= z1, ...
        r(2,:) >= -1];  % No sub-surface flight constraint
end

% Set up options and solve using Gurobi
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
sol = optimize(Constraints, -Objective, options);

% Check solver output
disp(sol)

% Detailed diagnosis
checkset(Constraints)

% Solver status and errors
yalmiperror(sol.problem)


% Check if the solution was successful
if sol.problem == 0
    % Extract data
    r_value = value(r);
    v_value = value(v);
    u_value = value(u);
    z_value = value(z);
    m_vals = exp(z_value);
    plot_run2D(tv, r_value, v_value, u_value, m_vals);
else
    disp(['Problem solving the optimization: ' yalmiperror(sol.problem)]);
end
