%% Constants and Initialization
g_I = [0; 0; -1.352]; % Titan gravity in m/s^2
m = 201.4; % Spacecraft mass in kg
N = 70; % Number of temporal nodes
tf = 76; % Time-of-flight in seconds
dt = tf / (N - 1); % Time step

% Initial conditions
r_I0 = [9; -7; 183]; % Starting position in meters
v_I0 = [9; 0; -8]; % Starting velocity in m/s

% Control limits
F_A_min = 251.7137; % Minimum force magnitude
F_A_max = 1.1215e+03; % Maximum force magnitude
M_A_min = 0.1; % Minimum moment magnitude
M_A_max = F_A_max; % Maximum moment magnitude
omega_max = deg2rad(30); % Maximum angular rate in rad/s
gamma_gs = deg2rad(20); % Glide-slope angle in radians
alpha_max = deg2rad(45); % Maximum angle of attack in radians
theta_max = deg2rad(90); % Maximum tilt angle in radians

% Spacecraft inertia matrix
I_B = diag([50, 50, 50]); % Inertia matrix of the spacecraft
n1 = [1; 0; 0]; % Unit vector along the body X-axis
n2 = [0; 1; 0]; % Unit vector along the body Y-axis
n3 = [0; 0; 1]; % Unit vector along the body Z-axis

%% Variables
F_A = sdpvar(3, N-1, 'full'); % Aerodynamic forces
M_A = sdpvar(3, N-1, 'full'); % Aerodynamic moments
delta_s = sdpvar(1, N-1, 'full'); % Control input scaling factors
r_I = sdpvar(3, N, 'full'); % Position
v_I = sdpvar(3, N, 'full'); % Velocity
q_BI = sdpvar(4, N, 'full'); % Orientation quaternions
omega_B = sdpvar(3, N, 'full'); % Angular velocity

%% Dynamics and Constraints
Constraints = [r_I(:,1) == r_I0, v_I(:,1) == v_I0, q_BI(:,1) == [1; 0; 0; 0], omega_B(:,1) == zeros(3,1)];

for i = 1:N-1
    % Kinematic equations
    Constraints = [Constraints,
        r_I(:,i+1) == r_I(:,i) + dt * v_I(:,i),
        v_I(:,i+1) == v_I(:,i) + dt * (F_A(:,i) / m + g_I)];
    
    % Force and moment constraints
    Constraints = [Constraints,
        F_A_min <= norm(F_A(:,i)) <= F_A_max,
        M_A_min <= norm(M_A(:,i)) <= M_A_max,
        0 <= delta_s(i) <= 1];
    
    % Glide-slope constraint
    Constraints = [Constraints,
        tan(gamma_gs) * norm([n2, n3]' * (r_I(:,i+1) - r_I0)) <= (n1' * (r_I(:,i+1) - r_I0))];
    
    % Tilt angle constraint
    Constraints = [Constraints,
        cos(theta_max) <= 1 - 2*(q_BI(2,i+1)^2 + q_BI(3,i+1)^2)];
    
    % Angle of attack constraint
    Constraints = [Constraints,
        v_I(:,i+1)' * n2 <= tan(alpha_max) * (v_I(:,i+1)' * n1)];
    
    % Angular rate constraint
    Constraints = [Constraints,
        norm(omega_B(:,i+1)) <= omega_max];


% Terminal constraints for soft landing
Constraints = [Constraints, r_I(:,N) == zeros(3,1), v_I(:,N) == zeros(3,1),
               omega_B(:,N) == zeros(3,1), q_BI(:,N) == [1; 0; 0; 0]];

%% Cost Function
Objective = sum(arrayfun(@(i) norm(F_A(:,i))^2 + norm(M_A(:,i))^2, 1:N-1));

%% Solver Setup and Execution
options = sdpsettings('solver', 'gurobi', 'verbose', 2, 'debug', 1);
sol = optimize(Constraints, Objective, options);

%% Post-Processing
if sol.problem == 0
    % Successful optimization
    r_val = value(r_I);
    v_val = value(v_I);
    F_val = value(F_A);
    
    % Plotting trajectory
    figure;
    plot3(r_val(1,:), r_val(2,:), r_val(3,:));
    grid on;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('Spacecraft Trajectory');
else
    % Handle errors
    disp(['Problem solving the optimization: ', yalmiperror(sol.problem)]);
end
end
