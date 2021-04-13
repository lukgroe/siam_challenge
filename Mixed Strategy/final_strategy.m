function simObj = final_strategy(simObj, lambda)
    d = simObj.d;
    T = simObj.T;
    eta = simObj.eta;
    
    w_next = zeros(d,1);            % Weights for the next step
    
    % Initialization of first weight vector.
    K = min(floor(1+3*lambda),d); 
    w = randperm(d);                % Get random permutation and choose the first K for assets.
    w = w(1:K);
    for i = 1:K
        w_next(w(i)) = 1/K;
    end
    
    simObj = simObj.reset();        % reset simulation environment
    mu_t = zeros(d,2);              % [mu approximation, number of data points used]
    cov = zeros(d,d);               % Estimator for covariance
    cov_counter = zeros(d,d);       % Counter of number of data points considered in cov
    cor = zeros(d,d);               % Estimator for correlation matrix
    
    % Simulate.
    for i=1:T
        
        % Differentiate between initial strategy and later strategy.
        if i<min(T/4, lambda*200) || mod(i,max(floor(50000*eta), 30)) ~= 0
            simObj = simObj.step(w_next);
        else
            
            for j=1:d
                for k=1:d % Calculation of correlation from covariance
                    cor(j,k) = cov(j,k)/sqrt(cov(j,j)*cov(k,k));
                end
            end
            
            goodPairs = cor <= zeros(simObj.d, simObj.d); % Returns logical matrix where (i,j)-th entry is 1 iff Cor(i,j) <= 0. 
            w_next = zeros(d,1);
            % Search for floor(1+3*lambda) pairs.
            for j=1:floor(1+3*lambda)
                maxInd = zeros(2,1);
                % Iterate over all pairs for a good pair.
                for k = 1:simObj.d
                    for m = k:simObj.d
                        if goodPairs(k,m)
                            if maxInd(1) == 0
                                maxInd = [k,m];
                            end
                            if mu_t(maxInd(1)) + mu_t(maxInd(2)) < mu_t(k) + mu_t(m)
                                maxInd = [k,m];
                            end
                        end
                    end
                end
                if maxInd(1) == 0 
                    % Check if all shares are positively correlated.
                    if j == 1 
                        % If so, leave strategy unchanged.
                        simObj = simObj.step(simObj.w_cur);
                    else 
                        % Otherwise just use the already picked pairs.
                        simObj = simObj.step(simObj.w_next);
                    end 
                    break
                end
                % Neutralize pair from being picked again.
                goodPairs(:, maxInd(1)) = zeros(d, 1);
                goodPairs(maxInd(1), :) = zeros(1, d);
                goodPairs(:, maxInd(2)) = zeros(d, 1);
                goodPairs(maxInd(2), :) = zeros(1, d);
                w_next(maxInd(1)) = 1;
                w_next(maxInd(2)) = 1;
            end
            w_next = w_next/sum(w_next);
            simObj = simObj.step(w_next);
        end
        
        % Compute difference of weights.
        if simObj.t == 1
            w_delta = zeros(d, 1);
        else
            w_delta = simObj.w_cur - simObj.w_hist(:,simObj.t-1);
        end
        
        % Prepare data by looking at the differences.
        diff = log(simObj.s_cur) - log(simObj.s_hist(:,simObj.t));
        
        % Calculate current mu approximation from this round.
        for j = 1:(d)
            if w_delta(j) == 0
                mu_t(j,1) = mu_t(j,2) * mu_t(j,1) + diff(j);
                mu_t(j,2) = mu_t(j,2) + 1;
                mu_t(j,1) = mu_t(j,1) / mu_t(j,2);
            end
        end
        
        % Calculate cov from current round        
        for j=1:d
            for k=1:d
                if w_delta(j) == 0 && w_delta(k) == 0
                    cov(j,k) = (cov(j,k)*max(1, cov_counter(j,k)-1) + (diff(j)-mu_t(j,1))*(diff(k)-mu_t(k,1)))/max(1,cov_counter(j,k));
                    cov_counter(j,k) = cov_counter(j,k)+1;
                end
            end
        end
    end
end
