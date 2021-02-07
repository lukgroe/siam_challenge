% Output file
fout = fopen('test.dat', 'w');

T = 200; % time steps
iterations = 10; %simulation iterations


% Vector containing number of shares for each iteration
%it_d = randi([2,20], 1, iterations);
it_d = 2.*ones(1,iterations);


% output variables, net_x is one share's w_hist and s_hist i each column
% net_mu, net_c is vector of each shares mu and c
columns = sum(it_d);
net_x=zeros(2*T+1, columns);
net_mu=zeros(1, columns);
net_c=zeros(1, columns);
writing_index = 1;

for i = 1:iterations
    % Model params initialization
    d = it_d(i);
    eta = rand / 20000;

    Mrank = floor(0.25*d);
    [U,S,V] = svd( randn(d,d) );
    diagM = diag( [ normrnd(0,1,Mrank,1) ; zeros(d-Mrank,1) ] );
    M = 5e-3 * U * diagM * V'; % Randomly generated matrix of rank Mrank

    %mu = 2e-5 * normrnd(0,1,d,1).^2;
    %c = 1e-8 * normrnd(0,1,d,1).^2;
    mu = [1;-1];
    c = [.5;.5];
    
    % starting prices
    s0 = randi(100)*ones(d,1);

    % Initialize Simulation Environment
    model_params = struct('mu',mu,'M',M,'c',c,'eta',eta);
    sim_obj = MarketSimulator(T,s0,model_params);

    % Run strategy on environment
    random_strategy(sim_obj);
    %example_strategy_1(sim_obj);
    %example_strategy_2(sim_obj);
    
    % writing to variables net_x, net_mu, net_c
    for j = 1:d % iterate over shares
        x_temp = [sim_obj.s_hist(j,:), sim_obj.w_hist(j,:)];
        for k = 1:2*T+1 % iterate over x_temp containing s_hist and w_hist
            net_x(k, writing_index) = x_temp(k);
        end
        net_mu(writing_index) = mu(j);
        net_c(writing_index) = c(j);
        writing_index = writing_index + 1;
    end
    
    % Output to test.dat
    fprintf(fout, '%d %d %d \n \n', [T, d, eta]);
    fprintf(fout, '%d %d \n', [mu, c]');
    fprintf(fout, '\n');
    fprintf(fout,[repmat(' %d ', 1, T+1) '\n'], sim_obj.s_hist');
    fprintf(fout, '\n');
    fprintf(fout,[repmat(' %d ', 1, T) '\n'], sim_obj.w_hist');
    fprintf(fout, '\n');
end
fclose(fout);

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