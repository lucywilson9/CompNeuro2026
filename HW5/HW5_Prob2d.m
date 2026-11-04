%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 5, Problem 2d
% November 2, 2026

%--------------------------------------------------------------------------
%% Binocular rivalry with adaptation
% Replace the heaviside by f = 1/(1+e^{-nx}) with n = 1000, integrate with
% a "stiff" method to t=30*ta, starting from (uL, uR, aL, aR) = (1,0,0,0).
% Call a percept dominant when its rate exceeds the other's more by than
% 1/2, and measure dominance durations between hand-overs, discarding the
% first three. 

% Plot uL, uR, aL, aR against time for IL = IR = 1.3, mark B and 1-B on the
% adaptation panel. 

% For IL = 1.3, 1.45, 1.6, and IR = 1.3 report the mean TL and TR and
% compare with part (c). 

% For w = 0.3 (below the window), report Tdom, compare with
% Ta(ln(w+I)/(g-w-I)). Report durations to 3 sig digits, comment on
% discrepancies.


% stiff method: BDF via ode15s
% params
w = 0.8; % sub 0.3 for other analysis
Be = 1;
g = 2;
ta = 120;
IL = 1.3; % sub 1.45, 1.6
IR = 1.3;
tspan = [0 30*ta];
y0 = [1 0 0 0];
[t,y] = ode15s(@(t,y) odefcn(t,y,w,Be,g,ta,IL,IR), tspan, y0);
uL = y(:,1);
uR = y(:,2);
aL = y(:,3);
aR = y(:,4);

B = (IL - Be)/g;

dom_vec = zeros(size(t));
dom_vec(uL - uR > 0.5) = 1; % TL timepoints
dom_vec(uR - uL > 0.5) = -1; % TR timepoints
for i = 2:length(dom_vec)
    if dom_vec(i) == 0
        dom_vec(i) = dom_vec(i-1); % assign to be the previous state
    end
end

handovers = find(diff(dom_vec) ~= 0); % find where the difference between two adjacent elements is non-zero

for i = 1:length(handovers)
    id = handovers(i);
    t1 = t(id); t2 = t(id+1);
    diff1 = uL(id) - uR(id);
    diff2 = uL(id+1) - uR(id+1);
    handovers_t(i) = t1 - diff1 * (t2 - t1) / (diff2 - diff1); % approx time where a switch happened
end

if length(handovers_t) > 3
    handovers = handovers_t(4:end);
    dom_ints = diff(handovers);
end

% sort dom_ints into TL, TR vectors
TL_vec = dom_ints(1:2:end);
avg_TL = mean(TL_vec);
TR_vec = dom_ints(2:2:end);
avg_TR = mean(TR_vec);

%Tdom = mean(dom_ints); % if w = 0.3
%pred_Tdom = ta*log((w+IL)/(g-w-IL));

cutoff_time_vec = t - handovers(1); % for getting rid of the first three ints
index = find(cutoff_time_vec == min(abs(cutoff_time_vec)));
threshold = t(index);

figure (1)
plot(t(index:end),uL(index:end), 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
plot(t(index:end), uR(index:end), 'LineWidth',1)
hold off
title('$u_L, u_R, I_L = I_R = 1.3$', 'interpreter','latex','FontSize', 14)
legend('u_L', 'u_R')
xlabel('time','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('rates','interpreter','latex', 'FontWeight','bold','FontSize',16)

figure (2)
plot(t(index:end),aL(index:end), 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
plot(t(index:end), aR(index:end), 'LineWidth',1)
plot([x(1) x(2)], [B B], 'k--', 'LineWidth',1)
plot([x(1) x(2)], [1-B 1-B], 'k*-', 'LineWidth',1)
hold off
title('$a_L, a_R, I_L = I_R = 1.3$', 'interpreter','latex','FontSize', 14)
legend('a_L', 'a_R', 'B', '1-B')
xlabel('time','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$a_L, a_R$','interpreter','latex', 'FontWeight','bold','FontSize',16)



% =========================================================================
function dydt = odefcn(~,y,w,Be,g,ta,IL,IR)
dydt = zeros(4,1);
dydt(1) = -y(1) + 1./(1 + exp(-1000.*(w*y(1)-Be*y(2)-g*y(3)+IL))); % uL
dydt(2) = -y(2) + 1./(1 + exp(-1000.*(w*y(2)-Be*y(1)-g*y(4)+IR))); % uR
dydt(3) = -y(3)/ta + y(1)/ta; % aL
dydt(4) = -y(4)/ta + y(2)/ta; % aR
end