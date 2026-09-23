%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 2, Problem 5c
% September 22, 2026

%--------------------------------------------------------------------------
%% Simulating a noisy IF model
% Simulate the model by Euler-Maruyama with reset, taking µ = θ = 1 and σ^2 = 0.25, 1, 4.
% Histogram the subthreshold voltage and overlay p(u) from part (a), and report the measured
% firing rate and the measured fraction of time below the reset against the predictions. 
%  Discard an initial transient before collecting statistics, and say why that matters here.

mu = 1;
sigma = 4; % 1, 4
T = 100;
T_burn = 10;
dt = 0.001;
urest = 0;
uth = 1; % theta

% using du/dt = mu + sqrt(2*sigma)*randn

% collecting code from lif_sim.m and noisy_if.m
nt = round(T/dt) + 1;
t  = linspace(0, T, nt)';
u  = zeros(nt,1);
u(1) = urest;
spikes = [];

for j = 1:nt-1
    du = dt*mu + sqrt(2*dt*sigma)*randn; % replace this ode with the mu one
    u(j+1) = u(j) + du;
    if u(j+1) > uth
        u(j+1) = urest;
        if t(j+1) >= T_burn
            spikes(end+1,1) = t(j+1);   %#ok<SAGROW>
        end
    end
end

if numel(spikes) >= 2
    rate = 1000*(numel(spikes)-1)/(spikes(end)-spikes(1));   % Hz
else
    rate = 0;
end


% implement burn-in
u = u(find(t==T_burn):end);

% next: define the probability func
p_u = @(u,theta,sigma) ((1-exp(-1/sigma))*exp((u./sigma))).*(u < 0) + (1-exp((u-1)./sigma)).*(u >= 0 & u <= theta);

p_for_hist = p_u(u, uth, sigma);

% collect some statistics to report compared to the expected stats:
% pred: firing rate
nu = mu/uth;
obs_fire = rate/1000;

pred_time_below = sigma*(1-exp(-1/sigma));
time_below = 0;
for i=1:length(u)
    if u(i) < 0
        time_below = time_below + dt;
    end
end

obs_time_below = time_below/(T-T_burn);

fprintf("pred rate: %8.3f\nobs rate: %.8f\n pred time: %.8f\n obs time: %.8f", nu, obs_fire, pred_time_below, obs_time_below)

figure (1)
histogram(u,'Normalization','pdf')
hold on
plot(u, p_for_hist, 'LineWidth',1)
hold off
xlabel("Voltage (mV)", 'FontWeight','bold','FontSize',14)
ylabel("pdf", 'FontWeight','bold','FontSize',14)
legend("u(t)", "p(u)")
title("Subthreshold voltage: $\sigma^2=4$",'interpreter','latex','FontWeight','bold','FontSize', 16)


%==========================================================================







