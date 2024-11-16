%% The following code on the LIF model is based on the Neuromatch Academy Tutorials that I have followed for reference 

% Define Basic Variables
t_max = 0.15;   % second
dt = 0.001;     % second
tau = 0.02;     % second
el = -60e-3;    % millivolt
vr = -0.07;     % millivolt
vth = -0.05;    % millivolt
r = 100000000;  % ohm
i_mean = 2.5e-10; % ampere

% Display values with specific formatting
disp('The fundamental Parameters that have been predefined are as follows:')
fprintf('t_max:%.2f, time step %.3f, Tau:%.2f, El:%.2f, Vr:%.2f, Vth:%.2f r:%.0f, I_mean %.1e\n', t_max, dt, tau, el, vr, vth, r, i_mean);
fprintf('\n'); 
% Display type of i_mean
%disp(class(i_mean)); Used for testing and validating code

%% Caluclate V for I approx as a sigmoid for 10 iterations 
% I is approximated as a sigmoid and V is calculated as per the LIF linear
% dif. equation 



% Initialize step_end and v0
step_end = 10;
v = el;   % Initial voltage

% Pre-allocate arrays for storing time, current, and voltage values
t_vals = zeros(step_end, 1);
i_vals = zeros(step_end, 1);
v_vals = zeros(step_end, 1);

% Loop for step_end steps
for step = 1:step_end
    % Compute value of t
    t = (step - 1) * dt;
    
    % Compute value of i at this time step
    i = i_mean * (1 + sin((t * 2 * pi) / 0.01));
    
    % Compute v using LIF equation
    v = v + (dt / tau) * (el - v + r * i); % el - v(t) - R*i

    % Store computed values
    t_vals(step) = t;
    i_vals(step) = i;
    v_vals(step) = v;
end

% Display the results in a table format (time, current, voltage)
disp('Given below are the values of Voltage for input current approx. as a sigmoid at different times:');fprintf('\n');
fprintf('%6s %12s %16s\n', 'Time(s)', 'Current(A)', 'Voltage(V)');
for step = 1:step_end
    fprintf('%.3f     %.4e      %.4e\n', t_vals(step), i_vals(step), v_vals(step));
end

%% Plotting Current and Voltage 
% Initialize step_end
step_end = floor(t_max / dt);

% Initialize v0 (membrane potential)
v = el;

% Pre-allocate arrays for time, current, and voltage values
t_vals = zeros(step_end, 1);
i_vals = zeros(step_end, 1);
v_vals = zeros(step_end, 1);

% Loop for step_end steps to compute both current and voltage
for step = 1:step_end
    % Compute value of t
    t = (step - 1) * dt;
    
    % Compute value of i (current) at this time step
    i = i_mean * (1 + sin((t * 2 * pi) / 0.01));
    
    % Compute v (membrane potential) using LIF equation
    v = v + (dt / tau) * (el - v + r * i);

    % Store computed values
    t_vals(step) = t;
    i_vals(step) = i;
    v_vals(step) = v;
end
% Plot Current vs Time 
figure;
plot(t_vals, i_vals, 'b-'); % 'b-' for blue solid line instead of dots to visulaise sigmoidal shape
title('Current vs Time');
xlabel('Time (s)');
ylabel('Current (A)');
legend('Current');
ylim([min(i_vals)-1e-10, max(i_vals)+1e-10]); % Adjust y-axis limits
grid on; % Optional, to add grid lines

% Plot Membrane Potential vs Time
figure;
plot(t_vals, v_vals, 'r.'); % 'r.' for red  so kinks are clearly visible 
title('Membrane Potential (V_m) vs Time');
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
legend('V_m');
grid on; % Adding grid lines 

%% Plotting voltage and Current for random synaptic input
% Set random number generator
rng(2024);

% Initialize step_end and v
step_end = t_max / dt;
v = el;

% Initialize the figure for V_m
figure;

% Prepare vectors for time and current for smoother plotting
time_vector = (0:dt:t_max-dt); % Create a time vector
current_values = zeros(1, step_end); % Preallocate current values array

% Loop for step_end steps
for step = 1:step_end
    % Compute value of t
    t = (step - 1) * dt;

    % Get random number in the range of -1 to 1
    random_num = 2 * rand() - 1; % random number between -1 and 1

    % Compute the value of i at this time step with random synaptic input
    i = i_mean * (1 + 0.1 * sqrt(t_max/dt) * random_num);
    current_values(step) = i; % Store the current value

    % Compute v (membrane potential) with random current input
    v = v + (dt / tau) * (el - v + r * i);

    % Plot v (membrane potential) using 'k.' to get small black dots
    plot(t, v, 'k.'); 
    hold on; 
end

title('V_m with Random I(t)', 'FontSize', 14); % Increase title font size
xlabel('Time (s)', 'FontSize', 12); % Increase label font size
ylabel('V_m (V)', 'FontSize', 12); % Increase label font size

% Add legend
legend('V_m');
grid on; % grid for better readability

% Display plot
hold off;

% Initialize the figure for Current
figure;

% Plot current values as a continuous line
plot(time_vector, current_values, 'k-', 'LineWidth', 1.5); % Use line for continuous effect
hold on;
% Plot the individual points for clarity
plot(time_vector, current_values, 'k.', 'MarkerSize', 6); % Overlay dots
% Add legend and customize
legend('I(t)');
grid on; % Add a grid for better readability
title('Current I(t) with Random Input', 'FontSize', 14); %title font size
xlabel('Time (s)', 'FontSize', 12); % Increase label font size
ylabel('Current I(t) (A)', 'FontSize', 12); % Increase label font size

% Display plot
hold off;


%% Simulating n= 50 different realisations of membrane potential (V_m), along with its sample mean and standard deviation
% Set random seed for reproducibility
rng(2024);

% Calculate total number of steps
step_end = round(t_max / dt); % Total number of steps, rounded for integer steps
n = 50; % Number of simulations
v_n = el * ones(n, 1); % Initialize with membrane potential

% Initialize figure
figure;
hold on;

% Loop through time steps
for step = 1:step_end
    % Compute current time step
    t = (step - 1) * dt;

    % Loop through each realization
    for j = 1:n
        % Compute current with random noise component
        i = i_mean * (1 + 0.1 * sqrt(t_max/dt) * (2 * rand - 1));

        % Update membrane potential for each simulation
        v_n(j) = v_n(j) + (dt / tau) * (el - v_n(j) + r * i);
    end

    % Compute sample mean and standard deviation after updating all realizations
    v_mean = mean(v_n);
    v_std = std(v_n);

    % Plot the realizations for the current time step
    plot(t * ones(n, 1), v_n, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); % Filled black circles

    % Plot sample mean for the current time step
    plot(t, v_mean, 'bo', 'MarkerSize', 8); % Mean as a blue point

    % Plot mean + standard deviation for the current time step
    plot(t, v_mean + v_std, 'go', 'MarkerSize', 6); % Mean + std as a green point

    % Plot mean - standard deviation for the current time step
    plot(t, v_mean - v_std, 'go', 'MarkerSize', 6); % Mean - std as a green point
end

% Finalise the plot
title('Multiple realizations of V_m');
xlabel('Time (s)');
ylabel('V_m (V)');
grid on;

% Add a legend to the plot
legend('Realizations', 'Mean', 'Mean + Std', 'Mean - Std');

% Display the plot
hold off;
