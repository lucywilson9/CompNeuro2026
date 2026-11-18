%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 6, Problem 2c
% November 6, 2026

%--------------------------------------------------------------------------
%% Foraging over patches
% Simulate 4 × 105 possessions with random outcomes, η = 0.02, a 
% running-average baseline, from θ = (0.3, −0.3). Add an entropy bonus by 
% putting −τθ p(1 − p) inside the bracket that multiplies η, and run τ = 0 
% and 0.05. (i) Plot the mean-flow orbit of (b) from (0.2, 0.3); report the
% range of x, y and the variation of K. (ii) Plot (x, y) for both runs; 
% report the last-quarter average and mean distance from (x^∗, y^∗) early 
% and late. (iii) Why does the τ = 0 learner drift outward, and what does 
% the bonus change?

% params/mixed eq
x_star = 8/17;
y_star = 9/17;
eta = 0.02;

% first: stuff from b
ode = @(t,z) [
    eta * (z(1)*(1-z(1)))^2*0.85*(z(2)-y_star); % x dot
   -eta * (z(2)*(1-z(2)))^2*0.85*(z(1)-x_star) % y dot
];

[~,Z] = ode45(ode,[0 10000],[0.2;0.3]);
x_traj = Z(:,1); % get x evolution
y_traj = Z(:,2); % get y evolution

K = @(s,p) (2*p-1).*log(s./(1 - s)) + (p - s)./(s.*(1 - s)); % integral from K def, need to add two to get K value
K_vals = (K(x_traj, x_star) - K(x_star, x_star)) + (K(y_traj, y_star) - K(y_star, y_star));

fprintf('Mean-Flow Orbit\n');
fprintf('x range: [%.4f, %.4f]\n', min(x_traj), max(x_traj));
fprintf('y range: [%.4f, %.4f]\n', min(y_traj), max(y_traj));
fprintf('K diff (max - min): %.3e\n', max(K_vals) - min(K_vals));

num_possessions = 4*10^5;
tau_vec = [0, 0.05];
rng(9); % seed

figure(1)
plot(x_traj, y_traj, 'b', 'LineWidth', 1.5) 
hold on
plot(x_star, y_star, 'ro', 'MarkerFaceColor', 'r','MarkerSize', 10, 'LineWidth', 2)
plot(0.2, 0.3, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10)
hold off
xlim([0 1])
ylim([0 1])
xlabel('x') 
ylabel('y')
title('Mean-Flow Orbit from (0.2, 0.3)')
legend('Orbit', 'Mixed Eq (x^*, y^*)', 'Initial Conditions')

for i = 1:length(tau_vec)
    tau = tau_vec(i);
    theta1 = 0.3;
    theta2 = -0.3;
    b1 = 0.0; % starting running baseline val
    b2 = 0.0; % starting running baseline val
    
    x_hist = zeros(num_possessions, 1); % initialize for storage
    y_hist = zeros(num_possessions, 1);
    
    for t = 1:num_possessions
        x = 1/(1 + exp(-theta1));
        y = 1/(1 + exp(-theta2));
        x_hist(t) = x;
        y_hist(t) = y;
        
        a1 = double(rand<x); % actions for a given possession
        a2 = double(rand<y); 
        
        % outcomes based on possible action combos
        if a1==0 && a2 == 1
            pts = 2*double(rand < 0.40); 
        elseif a1 == 0 && a2 == 0
            pts = 2*double(rand < 0.60);
        elseif a1 == 1 && a2 == 1
            pts = 3*double(rand < 0.40);
        else
            pts = 3*double(rand < 0.25);
        end
        
        r1 = pts;
        r2 = -pts; % defender rewared is minus points
        
        % "running-average baseline," update at each step
        b1 = b1+(r1-b1)/t;
        b2 = b2+(r2-b2)/t;
        
        theta1 = theta1+eta *((r1-b1)*(a1-x)-tau*theta1*x*(1-x)); % new term for entropy bonus
        theta2 = theta2 + eta * ((r2 - b2)*(a2 - y) - tau * theta2 * y * (1 - y));
    end
    
    q1 = 1:(N/4); % first quarter of "game"
    q4 = (3*N/4 + 1):N; % last quarter of "game"
    dist = sqrt((x_hist-x_star).^2+(y_hist-y_star).^2); % dist between the two points
    
    fprintf('Sim Run (tau = %.2f)\n', tau);
    fprintf('Last quarter avg (x, y): (%.4f, %.4f)\n', mean(x_hist(q4)), mean(y_hist(q4)));
    fprintf('Mean distance "early": %.4f\n', mean(dist(q1)));
    fprintf('Mean distance "late":  %.4f\n\n', mean(dist(q4)));
    
    % Plot trajectory (subsampled for clarity)
    figure(i + 1)
    plot(x_hist(1:20:end), y_hist(1:20:end), 'b', 'LineWidth', 1)
    hold on
    plot(x_star, y_star, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2)
    plot(mean(x_hist(q4)), mean(y_hist(q4)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7)
    hold off
    xlim([0 1])
    ylim([0 1])
    xlabel('x')
    ylabel('y')
    title(sprintf('Trajectory (\\tau = %.2f)', tau))
    legend('(x,y)', 'Mixed Eq (x^*, y^*)', 'Q4 Avg')
end
