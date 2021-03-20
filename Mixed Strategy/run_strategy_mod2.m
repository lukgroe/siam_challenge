%% Model and Simulator Initialization 

% Initialize Model Parameters
T = 500;
d = 20;
eta = 0.0002;

Mrank = floor(0.25*d);
[U,S,V] = svd( randn(d,d) );
diagM = diag( [ normrnd(0,1,Mrank,1) ; zeros(d-Mrank,1) ] );
M = 5e-3 * U * diagM * V'; % Randomly generated matrix of rank Mrank

mu = 2*10^(-5) * normrnd(0,1,d,1).^2;
c = 1e-8 * normrnd(0,1,d,1).^2;
%c=zeros(d,1);
s0 = 100*ones(d,1);

% Initialize Simulation Environment
model_params = struct('mu',mu,'M',M,'c',c,'eta',eta);
sim_obj = MarketSimulator(T,s0,model_params);

%% Visualization of a Single Simulation for a Strategy

% Run strategy on environment
%sim_obj = example_strategy_2(sim_obj); % proportional weight strategy
% sim_obj = example_strategy_1(sim_obj); % constant weight strategy
%sim_obj = simulation_strategy(sim_obj, mu, c, M);
%sim_obj = parameter_estimation_strategy(sim_obj, mu);
%sim_obj = mcts(sim_obj, mu, c);
[sim_obj, mu_t, c_t, cov, cor] = mixed_strategy(sim_obj, mu, c);

% Plot simulated price history
figure(1);
clf();
plot(1:(T+1),sim_obj.s_hist);
title('Stock Price Evolution')

% Plot portfolio weights
figure(2);
clf();
plot(1:T,sim_obj.w_hist);
title('Portfolio Weight Evolution')

% Plot portfolio 1-period returns + mean
figure(3);
clf();
hold on;
plot(1:T,sim_obj.r_hist);
plot(1:T,ones(1,T) * mean(sim_obj.r_hist))
hold off;
title('Portfolio 1-Period-Return Evolution')

% Plot portfolio cumulative growth
figure(4);
clf();
plot(1:T,sim_obj.R_hist-1);
title('Portfolio Cumulative Growth')


%% Computing the Target Objective for a Strategy
return; %stops the skript

nsims = 500;
lambda = 0.25;
cumret_array = zeros(nsims,1);

for k=1:nsims
    % Store each simulation's result in array
    sim_obj = example_strategy_2(sim_obj,lambda);
    cumret_array(k) = sim_obj.R_hist(end);
end


loss_value = mean(cumret_array) - 0.5*lambda*var(cumret_array);


