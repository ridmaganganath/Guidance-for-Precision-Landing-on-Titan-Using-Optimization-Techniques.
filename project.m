%% Initialize YALMIP and Clear Workspace
yalmip('clear');
clear;
clc;

%% Define Given Constants and Parameters
g_I = [0; 0; -1.352]; % Gravitational acceleration on Titan in m/s^2
m = 201.4; % Mass of the spacecraft in kg
alpha_max = 45; % Maximum angle of attack in degrees
gamma_gs = 20; % Glide-slope angle in degrees
theta_max = 90; % Maximum tilt angle in degrees
omega_max = deg2rad(30); % Maximum angular rate in rad/s
W_vc = 1e5; % Virtual control weight factor
W_tr = 1e-3; % Trust region weight factor
tol_vc = 1e-9; % Virtual control tolerance
tol_tr = 1e-3; % Trust region tolerance
max_iter = 10; % Maximum number of iterations
N = 70; % Number of temporal nodes
tf = 76; % Time-of-flight in seconds
r_I0 = [9; -7; 183]; % Starting position in meters
v_I0 = [9; 0; -8]; % Starting velocity in m/s
F_A_min_sq = 1.6669^2;  % Squared minimum force magnitude
F_A_max_sq = 3.3919^2;  % Squared maximum force magnitude
M_A_min_sq = 0.1^2;  % Squared minimum moment magnitude
M_A_max_sq = 0.5^2;  % Squared maximum moment magnitude
I_B = diag([50, 50, 50]);  % Inertia matrix of the spacecraft

%% Define Variables
% Control variables
F_A = sdpvar(3, N-1); % Aerodynamic force vectors for each time step
M_A = sdpvar(3, N-1); % Moment vectors for each time step
delta_s = sdpvar(1, N-1); % Control input scaling factors

% State variables
r_I = sdpvar(3, N); % Position in inertial frame
v_I = sdpvar(3, N); % Velocity in inertial frame
q_BI = sdpvar(4, N); % Quaternion representing orientation
omega_B = sdpvar(3, N); % Angular velocity in body frame

% Initial conditions
constraints = [r_I(:,1) == r_I0, v_I(:,1) == v_I0, q_BI(:,1) == [1; 0; 0; 0], omega_B(:,1) == zeros(3,1)];

%% Define Dynamics and Constraints
objective = 0;
for k = 1:N-1
    dt = tf/N;  % Time step
    
    % Position and velocity dynamics
    r_next = r_I(:,k) + dt * v_I(:,k);
    v_next = v_I(:,k) + dt * (F_A(:,k)/m + g_I);
    
    % Ensure quaternion is a 1x4 vector for multiplication
    current_q = q_BI(:,k)';  % Transpose to make it a 1x4 vector

    % Create a quaternion from angular velocity, ensuring it's 1x4
    omega_quat = [0, omega_B(:,k)'];  % Include zero for real part, make it 1x4

    % Quaternion dynamics using quaternion multiplication
    % Calculate the quaternion derivative
    q_dot = 0.5 * quatmultiply(current_q, omega_quat);

    % Update quaternion by integrating the derivative
    q_next = current_q + q_dot * dt;

    % Normalize the quaternion to prevent drift
    q_next = q_next / norm(q_next);

    % Update quaternion back to 4x1 column vector for storage
    q_BI(:,k+1) = q_next';

    % Angular velocity dynamics
    omega_next = omega_B(:,k) + dt * (I_B \ (M_A(:,k) - cross(omega_B(:,k), I_B * omega_B(:,k))));

    % Update state variables
    constraints = [constraints, ...
                   r_I(:,k+1) == r_next, ...
                   v_I(:,k+1) == v_next, ...
                   omega_B(:,k+1) == omega_next];

    % Force and moment constraints
    constraints = [constraints, ...
                   F_A_min_sq <= sum(F_A(:,k).^2) <= F_A_max_sq, ...
                   M_A_min_sq <= sum(M_A(:,k).^2) <= M_A_max_sq, ...
                   0 <= delta_s(k) <= 1];

    % Objective (minimize control effort)
    objective = objective + sum(F_A(:,k).^2) + W_vc * delta_s(k)^2;
end


%% Boundary Conditions
final_conditions = [r_I(:,N) == zeros(3,1), v_I(:,N) == zeros(3,1), q_BI(:,N) == [1; 0; 0; 0], omega_B(:,N) == zeros(3,1)];
constraints = [constraints, final_conditions];

%% Setup Solver Options and Solve
options = sdpsettings('solver', 'gurobi', 'verbose', 1, 'gurobi.TimeLimit', 100);  %// Increased time limit

% Solving the problem
sol = optimize(constraints, objective, options)

% Check results
if sol.problem == 0
    disp('Solution Found');
    disp(value([r_I, v_I]));
else
    disp('No solution');
    disp(yalmiperror(sol.problem));
end
