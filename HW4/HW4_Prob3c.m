%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 4, Problem 3c
% October 19, 2026

%--------------------------------------------------------------------------
%% Simulating a Hopfield network
% Simulate a Hebbian network with N = 500 for α = 0.02, 0.06, . . . , 0.62. 
% For each α draw a pattern set, and measure the fraction of unstable bits, ξµihµi < 0 
% over all i and µ, and compare it with Φ(−1/√α). Then test recall: start 
% from a stored pattern with 10% of its bits flipped, update asynchronously 
% in random order until a full sweep changes nothing, and record the overlap m =1/Nξµ· s 
% with the original. Average over 10 cues per α, plot hmi against α, and say where recall breaks down.

alpha_vec = 0.02:0.04:0.62;
P_vec = round(500*alpha_vec);

% task one: "draw a pattern set": need to randomly make (and save?)
frac_unstable_vec = zeros(length(alpha_vec),1);
cdf_vec = zeros(length(alpha_vec),1);
overlap_vec = zeros(length(alpha_vec), 1);

for i = 1:length(P_vec)
    P = P_vec(i);
    cdf_vec(i) = normcdf(-1 / sqrt(alpha_vec(i)));
    xi = sign(rand(N,P)-0.5); % combine everything into the NxP matrix
    xi(xi == 0) = 1; % correct for any accidental 0 entries
    W = (xi * xi') / N;
    W = W - diag(diag(W));
    H = W * xi;
    unstable_bits = (xi .* H) < 0;
    frac_unstable_vec(i) = mean(unstable_bits(:)); % looks at unstable bits as a vector, computes fraction

    % recall stuff
    m = 0;
    for j = 1:10 % number of cues
        s = xi(:, 9); % pick ninth pattern as the one to examine.
        flipped = randi([1, 500], 1, 50);
        s(flipped) = -s(flipped); % flip 10% of the entries, 

        fixed = false;
        while ~fixed
            fixed = true;
            update_order = randperm(N); % randomize the order
            for k = update_order
                h_i = W(k, :) * s; 
                s_new = sign(h_i); 
                if s_new == 0 
                   s_new = 1; % correct for improbable zero again
                end
                
                if s_new ~= s(k)
                    s(k) = s_new;
                    fixed = false; % need to iterate again if a bit doesn't flip back
                end
            end
        end
        
       m = m + (xi(:, 9)' * s) / N;% Calculate final overlap
    end
    overlap_vec(i) = m / 10;
end

figure (1)
plot(alpha_vec, frac_unstable_vec, 'o-', 'LineWidth', 1)
hold on
plot(alpha_vec, cdf_vec, '*-', 'LineWidth', 1)
hold off
xlabel('\alpha') 
ylabel('Fraction of Unstable Bits')
legend('Observed', 'Theoretical')
title('Unstable Bits in Pattern Simulation')

figure(2)
plot(alpha_vec, overlap_vec, '*-', 'LineWidth', 2)
xlabel('\alpha'); 
ylabel('\langle m \rangle');
title('Average Recall Overlap')