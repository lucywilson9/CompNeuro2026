%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 4, Problem 2d
% October 19, 2026

%--------------------------------------------------------------------------
%% Competitive decision making 
% Simulate the two-population model with µ = 0.2, σ = 0.5, I0 = 5, ws = wc = (1 + λ)/2 for λ = 0 and λ = 0.2,
% starting at r1 = r2 = I0. Use Euler–Maruyama with dt = 0.002 and 10^4 trials for each θ ∈ {0.5, 1, . . . , 4}. 
% Plot accuracy and mean decision time against θ together with the predictions from (b) and (c). At λ = 0.2 and θ = 4,
% compare the accuracy and decision time with those of the perfect integrator at the threshold giving the same accuracy

% parameter stuff:
lambda_vec = [0 0.2];

[mean_time_vec_1, prob_vec_1] = popsim(lambda_vec(1));
[mean_time_vec_2, prob_vec_2] = popsim(lambda_vec(2));

theta_vec = 0.5:0.5:4;
% predicted vals:

p_correct_b = @(theta) 1./(1+exp((-0.2.*theta)/0.5^2));
mean_t_b = @(theta) (theta./0.2).*tanh((0.2.*theta)/(2*(0.5^2)));

gfun = @(y) exp((-0.2*y.^2 - 2*0.2*y)/(2*0.5^2));
p_correct_c = @(theta) (integral(gfun, -theta,0))/(integral(gfun,-theta,theta));

b_prob_pred = p_correct_b(theta_vec);
c_prob_pred = [];
for i = 1:length(theta_vec)
    theta = theta_vec(i);
    c_prob_pred(i) = p_correct_c(theta); %#ok<SAGROW>
end
mean_time_pred = mean_t_b(theta_vec);

% perfect integrator, "accuracy" for lambda = 0.2, theta = 4;
perf_acc = cdf('Normal', 0.2/(0.5*sqrt(0.2)), 0, 1);

figure (1)
plot(theta_vec, prob_vec_1, '*-', 'LineWidth',1)
hold on
plot(theta_vec, b_prob_pred, '*-', 'LineWidth',1)
hold off
xlabel("Threshold value", 'FontWeight','bold','FontSize',14)
ylabel("Measured Val", 'FontWeight','bold','FontSize',14)
legend("Observed Accuracy", "Predicted Accuracy")
title("Two population simulation: $\lambda=0$",'interpreter','latex','FontWeight','bold','FontSize', 16)

figure (2)
plot(theta_vec, prob_vec_2, '*-', 'LineWidth',1)
hold on
plot(theta_vec, c_prob_pred, '*-', 'LineWidth',1)
hold off
xlabel("Threshold value", 'FontWeight','bold','FontSize',14)
ylabel("Measured Val", 'FontWeight','bold','FontSize',14)
legend("Observed Accuracy", "Predicted Accuracy")
title("Two population simulation: $\lambda=0.2$",'interpreter','latex','FontWeight','bold','FontSize', 16)

figure (3)
plot(theta_vec, mean_time_vec_1, 'o-', 'LineWidth',1)
hold on
plot(theta_vec, mean_time_pred, 'o-', 'LineWidth',1)
hold off
xlabel("Threshold value", 'FontWeight','bold','FontSize',14)
ylabel("Measured Val", 'FontWeight','bold','FontSize',14)
legend("Observed Mean Time", "Predicted Mean Time")
title("Two population simulation: $\lambda=0$",'interpreter','latex','FontWeight','bold','FontSize', 16)

figure (4)
plot(theta_vec, mean_time_vec_2, 'o-', 'LineWidth',1)
hold on
plot(theta_vec, mean_time_pred, 'o-', 'LineWidth',1)
hold off
xlabel("Threshold value", 'FontWeight','bold','FontSize',14)
ylabel("Measured Val", 'FontWeight','bold','FontSize',14)
legend("Observed Mean Time", "Predicted Mean Time")
title("Two population simulation: $\lambda=0.2$",'interpreter','latex','FontWeight','bold','FontSize', 16)



%==========================================================================
function [mean_time_vec, prob_vec] = popsim(lambda)
mu = 0.2;
sigma = 0.5;
I0 = 5;
ws = (1+lambda)/2;
wc = (1+lambda)/2;
r10 = I0;
r20 = I0;
num_trials = 10^4;
theta_vec = 0.5:0.5:4;
dt = 0.002;

prob_vec = [];
mean_time_vec = [];
for i = 1:length(theta_vec)
    theta = theta_vec(i);

    results_vec = [];
    time_vec = [];

    for j = 1:num_trials
    r1 = r10;
    r2 = r20;
    d = 0;
    t = 0;
    while abs(d) < theta
        dr1 = dt*(-r1 + ws*r1 - wc*r2 + I0 + mu/2) + sqrt(2*dt*sigma)*randn; %sqrt(2*dt*sigma)*randn
        dr2 = dt*(-r2 + ws*r2 - wc*r1 + I0 - mu/2) + sqrt(2*dt*sigma)*randn;
        r1 = r1 + dr1;
        r2 = r2 + dr2;
        d = r1 - r2;
        t = t + dt;
    end

    if d >= theta
        results_vec(j) = 1;  %#ok<AGROW>
    elseif d <= -theta
        results_vec(j) = -1; %#ok<AGROW>
    end

    time_vec(j) = t; %#ok<AGROW>
    end

    prob_vec(i) = sum(results_vec == 1)/num_trials; %#ok<AGROW>
    mean_time_vec(i) = mean(time_vec); %#ok<AGROW>
end
end



