% Set random seed for reproducibility
rng(2024);

% Define Basic Variables
t_max = 0.15;   % seconds
dt = 0.001;     % seconds
tau = 0.02;     % seconds
el = -60e-3;    % millivolts
vr = -0.07;     % millivolts (reset potential)
vth = -0.05;    % millivolts (threshold potential)
r = 100000000;  % ohms
input_currents = [2.5e-10, 3e-10]; % Different input currents (in amperes)

% Initialize figure
figure;
hold on;

% Loop through different input currents
for i = 1:length(input_currents)
    i_mean = input_currents(i); % Current for this simulation
    step_end = round(t_max / dt); % Total number of steps

    % Initialize membrane potential
    v_n = el; % Start at membrane leak potential
    v_record = []; % Record potential for plotting
    time_record = []; % Record time for plotting
    threshold_crossings = 0; % Count how many times threshold is crossed

    % Loop through time steps
    for step = 1:step_end
        t = (step - 1) * dt; % Compute current time step

        % Update membrane potential with current
        v_n = v_n + (dt / tau) * (el - v_n + r * i_mean);
        
        % Record the potential and time for plotting
        v_record = [v_record; v_n];
        time_record = [time_record; t];

        % Check for threshold crossing
        if v_n >= vth
            threshold_crossings = threshold_crossings + 1; % Increment crossing count
            v_n = vr; % Reset potential after threshold crossing
            if threshold_crossings >= 2 % Break after crossing threshold twice
                break;
            end
        end
    end
    
    % Continue to rise again after threshold resets
    for step = length(v_record)+1:step_end
        t = (step - 1) * dt; % Compute current time step

        % Keep updating with the same current after reset
        v_n = v_n + (dt / tau) * (el - v_n + r * i_mean);
        
        % Record the potential and time for plotting
        v_record = [v_record; v_n];
        time_record = [time_record; t];

        % Check for another threshold crossing
        if v_n >= vth
            threshold_crossings = threshold_crossings + 1; % Increment crossing count
            v_n = vr; % Reset potential after crossing
            if threshold_crossings >= 3 % Break after three crossings
                break;
            end
        end
    end
    
    % Plot results for this current
    plot(time_record, v_record * 1000, 'DisplayName', sprintf('I_{mean} = %.1e A', i_mean)); % Convert to millivolts for plotting
end

% Finalize the plot
title('Membrane Potential Rise and Fall with Different Input Currents');
xlabel('Time (s)');
ylabel('Membrane Potential V_m (mV)');
legend('show'); % Show legend
grid on;

% Set axis limits for clarity
ylim([-70 -40]); % Adjusted Y-axis limits to clearly show potential behavior
xlim([0 t_max]); % X-axis limits to match maximum time
hold off;
