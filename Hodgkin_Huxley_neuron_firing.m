%% Source: I have followed the instructions in andysbrainblog.blogspot.com which was a very helpful resource 
% The code below is reproduced from the blog with certain parameters changed
%% Setting Constants and Simulation Parameters 
% simulation time
simulationTime = 100; %in milliseconds
deltaT=.01;
t=0:deltaT:simulationTime; % to store timepoints in simulation 


% specify the external current I that changes over time (50 (500 ms), 0 (1500 ms) then 50 again)
changeTimes = [0]; %in milliseconds
currentLevels = [50]; %Change this to see effect of different currents on voltage (values chosen: 3, 20, 50, 1000)

%Set externally applied current across time
%Here, first 500 timesteps are at current of 50, next 1500 timesteps at
%current of zero (resets resting potential of neuron), and the rest of
%timesteps are at constant current
I(1:500) = currentLevels; I(501:2000) = 0; I(2001:numel(t)) = currentLevels;
%Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
%I(1:numel(t)) = currentLevels;


%constant parameters from can be found in Table 3 of paper
gbar_K=36; gbar_Na=120; g_L=.3;
E_K = -12; E_Na=115; E_L=10.6;
C=1;


%% set the initial states- based on stead state of HH model (all equation numberings wrt paper)

V=0; %Baseline voltage
alpha_n = 0.01 * ( (10-V) / (exp((10-V)/10)-1) ); %Equation 12
beta_n = 0.125*exp(-V/80); %Equation 13
alpha_m = 0.1*( (25-V) / (exp((25-V)/10)-1) ); %Equation 20
beta_m = 4*exp(-V/18); %Equation 21
alpha_h = 0.07*exp(-V/20); %Equation 23
beta_h = 1/(exp((30-V)/10)+1); %Equation 24

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18

%% Main Loop: Compute coefficients, currents, and derivates at each time step
for i=1:numel(t)-1
   
    %calculate the coefficients
    % Same equations as above but calculating at each time step
    alpha_n(i) = 0.01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = 0.125*exp(-V(i)/80);
    alpha_m(i) = 0.1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = 0.07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
   
   
    %calculate the currents (Ionic current is external current - individual currents)
    I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L = g_L *(V(i)-E_L); %Equation 5
    I_ion = I(i) - I_K - I_Na - I_L;
   
   
    % calculate the derivatives using Euler's method (first order) of numerical
    % integration 
    V(i+1) = V(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16

end


V = V-70; %Set resting potential to -70mv

%% Plotting 

%===plot Voltage===%
plot(t,V,'LineWidth',3)
hold on
legend({'voltage'})
ylabel('Voltage (mv)')
xlabel('time (ms)')
title('Voltage over Time in HH Neuron')


%===plot Conductance===%
figure
p1 = plot(t,gbar_K*n.^4,'LineWidth',2);
hold on
p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
ylabel('Conductance (S)')
xlabel('time (ms)')
title('Conductance for Potassium and Sodium Ions in HH Neuron')
