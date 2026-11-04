%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 5, Problem 3d
% November 2, 2026

%--------------------------------------------------------------------------
%% Bumps for a Mexican Hat
% Simulate the field on [-20, 20] with deltax = 0.02, an explicit euler
% step deltat = 0.02, FFT convolution with zero padding, and H replaced by
% a steep sigmoid 1/(1+e^{-100(u-k)}). For k = 0.3, start from u(x,0) =
% W(x+a0) - W(x-a0) with a0 = a_ \pm 0.02 and a0 = a+ \pm 0.05, plot the
% half-width of the active region against time. Which return, which run
% away, and which die? 
% For k = 0.1, 0.3, 0.45, start the wide bump 0.05 too wide and fit the
% relaxation rate of |a(t) - a+| on a log scale; compare with lambda_e(a+).


% params
T = 20; 
dx = 0.02;
dt = 0.02;
x = (-T:dx:T)'; % spatial grid
N = length(x);
tmax = 15;
nsteps = round(tmax/dt); % temporal grid
t_vec = (1:nsteps)*dt;
    
%k = 0.3;
k_vec = [0.1 0.3 0.45];
%S = sqrt(1 - 2*k);
%a_minus = -log(0.5*(1 + S)); % logging the exp formulae from prob statement
%a_plus  = -log(0.5*(1 - S));

% initial conditions
%a0_vec = [a_minus - 0.02, a_minus + 0.02]; % initial conditions close to unstable width
%a0_vec = [a_plus - 0.05, a_plus + 0.05]; % initial conditionsn close to stable width

%lambda_e_a_plus = -S*(1-S)/(1+((1-S)*S)/2);
          

% FFT with zero padding????

% case-sensitive 
w_func = @(x) 2*exp(-abs(x))-exp(-abs(x)/2);
W_func = @(x) sign(x).*(2*exp(-abs(x)/2)-2*exp(-abs(x)));
    
x_vec = [0:dx:T, -T:dx:-dx]'; % shift x grid to start at 0
w_vec = w_func(x_vec);
W_fft = fft(w_vec) * dx; 

relax_vec = zeros(3*nsteps, 1);
for i = 1:length(k_vec)
    figure(1)
% for i = 1:length(a0_vec)
    k = k_vec(i);
    S = sqrt(1 - 2*k);
    a_plus  = -log(0.5*(1 - S));
    a0 = a_plus + 0.05;
    lambda_e_a_plus = -S*(1-S)/(1+((1-S)*S)/2)
    %a0 = a0_vec(i);
    u = W_func(x+a0)-W_func(x-a0); % u(x,0), this is a vector    
    half_widths = zeros(nsteps, 1); % will compute the half width at each time step
    relax_vec = zeros(nsteps, 1);
        
    for j = 1:nsteps
        H = 1./(1+exp(-100*(u-k)));            
        H_fft = fft(H);
        conv = real(ifft(W_fft.*H_fft)); % FFT of product is product of FFTs  
        u = u + dt * (-u + conv); % explicit Euler
            
        active_indices = find(u>= k); % heaviside-esque
        if isempty(active_indices)
           half_widths(j) = 0;
        else
           half_widths(j) = (x(active_indices(end)) - x(active_indices(1))) / 2; % finds a, where range of x points where active is|x| <= a
        end
        half_widths(j);
        %relax_vec(j+((i-1)*nsteps)) = abs(half_widths(j) - a_plus);
        relax_vec(j) = abs(half_widths(j) - a_plus);
   end
        
     %plot(t_vec, half_widths, 'LineWidth', 2)
    log_relax = log10(relax_vec);
    p = polyfit(t_vec, log10(relax_vec), 1);
    logrelax_fit = polyval(p, t_vec);
    slope = p(1);
    intercept = 10^p(2);
    relax_fit = 10.^logrelax_fit;

        semilogy(t_vec, relax_vec, '-', 'LineWidth', 2)
        hold on
        semilogy(t_vec, relax_fit, '-', 'LineWidth',2)
end
    
% first_k = relax_vec(1:nsteps);
% second_k = relax_vec((nsteps+1):2*nsteps);
% third_k = relax_vec((2*nsteps+1):end);


    %hold off
    % need to add line of expected rate of conv.
    %figure (2)
    %loglog(t_vec, second_k, 'LineWidth',1.5)
    %hold on
    %loglog(t_vec, second_k, 'LineWidth',1.5)
    %loglog(t_vec, third_k, 'LineWidth',1.5)
   % hold off
    hold off
    xlabel('time (t)')
    ylabel('$|a(t) - a_+|$', 'Interpreter','latex')
    title('Relaxation Rate vs. time, varying \kappa')
    legend('k = 0.1', 'k=0.1 fit','k=0.3', 'k=0.3 fit', 'k=0.45', 'k=0.45 fit')
