%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 3, Problem 2c
% October 4, 2026

%--------------------------------------------------------------------------
%% Measuring PRCs numerically
% at phase φ kick the voltage by a small δ (say 10−3), record the time T0 of the next spike, 
% and compute R(φ) ≈ (∆ − T0)/(∆δ), repeating on a grid of φ. 
% Firstdo this for the LIF of part (a) with τ = 1, vr = 0, vth = 1, I = 1.5, and overlay your formula. 
% Then do it for the FitzHugh–Nagumo oscillator
% kicking v, with phase zero at the upward crossing of v = 0. Here there is no reset, so measure the shift in
% a crossing a few cycles later, after the perturbation has relaxed back to the limit cycle. Report the period,
% the extreme values of R, and where R < 0. Which PRC is Type I and which is Type II, and what does a
% negative lobe mean physically?

tau = 1;
vr = 0;
vth = 1;
I = 1.5;
kicksize = 10^(-3);
dt = 0.001;

T = lif_period(I, dt, vth); 

phi_range_unscaled = 0:dt:T;
phi_range_scaled = phi_range_unscaled/T;
Z_LIF = zeros(length(phi_range_scaled),1);

% code adapted from prc.m
for i = 1:length(phi_range_scaled) % going through each point on the phi grid
    u = 0; t = 0; kicked = false; % simulate the voltage
    while u <= vth
        u = u + dt*(I - u); % numerically solve the LIF model
        t = t + dt;
        if ~kicked && t >= phi_range_unscaled(i) % once we get to a time point past the unscaled phi point
            u = u + kicksize; kicked = true; % kick the voltage, then resume solving the model until we hit a spike
        end
    end
    Z_LIF(i) = (T - t)/(T*kicksize); % t is the spike time, so this simulates the PRC according to the provided formula
end

comp_LIF_PRC = @(phi, delt, vr) (tau*exp((phi*delt)/tau))/(delt*(I*tau-vr)); % derived PRC from HW/text
comp_LIF_PRC_scaled = comp_LIF_PRC(phi_range_unscaled, T, vr);

% FH-N
opts   = odeset('RelTol',1e-9,'AbsTol',1e-11);
[t, Y] = ode45(@(tt,y) fn_rhs(y, 0.5, 0.08, 0.8, 0.7), 0:dt:20, [0.1, 0], opts);
u = Y(:,1); w = Y(:,2); % this gives us a vector to use for w

[minval, index] = min(abs(u));
first_cross = t(index);

[t, Y] = ode45(@(tt,y) fn_rhs(y, 0.5, 0.08, 0.8, 0.7), 0:dt:60, [0.1, 0], opts);
u = Y(:,1); w = Y(:,2); 
[minval, index] = min(abs(u));
next_cross = t(index);

fhn_phi_range_unscaled = 0:500*dt:(next_cross-first_cross);
T1 = fhn_phi_range_unscaled(end);
fhn_phi_range_scaled = fhn_phi_range_unscaled/(next_cross-first_cross);


Z_FHN = zeros(length(fhn_phi_range_scaled),1);
FHN_kicksize = 0.5;
for i = 1:length(fhn_phi_range_scaled) % going through each point on the phi grid
    u_old = 0.1; t = 0; kicked = false; % simulate the voltage
    u_vec = [];
    t_vec = [];
    u_vec(1) = u_old;
    t_vec(1) = t;
    w_old = 0; % initial val for w
    crossings = [];
    j=1;
    while (t <= 250)
        u = u_old + dt*(u_old - u_old^3/3 - w_old + 0.5); % numerically solve the LIF model
        w = w_old + dt*(0.08*(u_old + 0.7 - 0.8*w_old)) ;
        t = t + dt;
        if ~kicked && t >= fhn_phi_range_unscaled(i) % once we get to a time point past the unscaled phi point
            % fprintf('Voltage pre-kick: %8.3f\n', u)
            u = u + FHN_kicksize; kicked = true; % kick the voltage, then resume solving the model until we hit a spike
            % fprintf('Voltage was kicked! At t=%8.3f\n', t)
            % fprintf('Voltage post-kick: %8.3f\n', u)
        end
        u_old = u;
        w_old = w;
        u_vec(end+1) = u_old;
        t_vec(end+1) = t;
        j= j+1;
        if abs(u_old) < 3*10^(-4) && t >= fhn_phi_range_unscaled(i) && u_vec(j-1) > 0
            crossings(end+1) = t; 
        end
    end
    T_new = crossings(2) - crossings(1); % sample the period in which the kicker kicked
    Z_FHN(i) = (T1 - T_new)/(T1);
end

fprintf('max est. FH-N PRC: %8.3f\n, min est. FH-N PRC: %8.3f\n', max(Z_FHN), min(Z_FHN))


% LIF plot:
figure(1)
plot(phi_range_scaled, Z_LIF)
y = ylim;
x = xlim;
hold on
plot(phi_range_scaled, comp_LIF_PRC_scaled)
%plot([left_bdry left_bdry], [y(1) y(2)], 'r--', 'LineWidth',1)
%plot([right_bdry right_bdry], [y(1) y(2)], 'b--', 'LineWidth',1)
hold off
xlabel("Phase $\phi$", 'interpreter','latex','FontWeight','bold','FontSize',14)
ylabel("R(\phi)", 'FontWeight','bold','FontSize',14)
legend("Measured PRC", "Expected PRC")
title("Measured vs. Expected PRC for the LIF Model",'interpreter','latex','FontWeight','bold','FontSize', 16)

% FH-N Plot
figure(2)
%plot(t, u)
plot(fhn_phi_range_scaled, Z_FHN)
y = ylim;
x = xlim;
hold on
%plot(phi_range_scaled, comp_LIF_PRC_scaled)
%plot([left_bdry left_bdry], [y(1) y(2)], 'r--', 'LineWidth',1)
%plot([right_bdry right_bdry], [y(1) y(2)], 'b--', 'LineWidth',1)
hold off
xlabel("Phase $\phi$", 'interpreter','latex','FontWeight','bold','FontSize',14)
ylabel("Z(\phi)", 'FontWeight','bold','FontSize',14)
%legend("Measured Firing Numbers", "Expected Firing Numbers", 'Left tongue bdry, $\sim 1.5034', 'Right tongue bdry, $\sim 1.6606')
title("Measured PRC for the FH-N Model",'interpreter','latex','FontWeight','bold','FontSize', 16)


%==========================================================================
function T = lif_period(I, dt, thresh) % adapted from prc.m
u = 0; T = 0; % note we use u= 0 as a starting point, not u= -thresh from prc.m
while u <= thresh
    u = u + dt*(I - u);
    T = T + dt;
end
end

% FitzHugh-Nagumo model: code from fn_sim.m
function dy = fn_rhs(y, I, eps, gamma, b0) % encodes the 2d system
dy = [y(1) - y(1)^3/3 - y(2) + I; eps*(y(1) + b0 - gamma*y(2)) ];
end

