function [simObj, mu_t, c_t, cov, cor] = mixed_strategy2(simObj, lambda)
    d = simObj.d;
    T = simObj.T;
    eta = simObj.eta;
    w_next = zeros(d,1);
    w = min(floor(1+3*lambda),1/d);
    while sum(w_next)<1-w % Random initial strategy; equal weights over floor(1+3*lambda) random shares
        w_next(randi([1,d])) = w;
    end
    simObj = simObj.reset(); % reset simulation environment
    s_w_delta_nonzero = zeros(d, T); % Delta of s for delta w not equal zero
    w_w_delta_nonzero = zeros(d, T); % Delta of w for delta w not equal zero
    w_nonzero_counter = ones(d,1); % Counter for number of data points
    mu_t = zeros(d, 2); % [mu approximation, number of data points used]
    c_t = zeros(d, 1); % Estimator for c
    cov = zeros(d, d); % Estimator for covariance
    cov_counter = zeros(d, d); % Counter of number of data points considered in cov
    cor = zeros(d, d); % Estimator for correlation matrix
    
    % Simulate.
    for i=1:T
        
        if i<min(T/4, lambda*200) || mod(i,max(50000*eta, 30)) ~= 0
            simObj = simObj.step(w_next);
        else
            
            c_t = zeros(d, 1); %currently unnecessary calculation of c
            for j=1:d
                for k=1:d % Correlation is being calculated
                    cor(j,k) = cov(j,k)/sqrt(cov(j,j)*cov(k,k));
                end
                if w_nonzero_counter(j)-1 ~= 0
                    for k=1:(w_nonzero_counter(j)-1)
                        c_t(j) = c_t(j) + (s_w_delta_nonzero(j,k)-mu_t(j,1))/(sign(w_w_delta_nonzero(j,k)*sqrt(abs(w_w_delta_nonzero(j,k))))); % Adding up the differences
                    end
                    c_t(j) = c_t(j)/w_nonzero_counter(j)-1; % Dividing by number of data points
                end
            end
            
            
            w_next = zeros(d,1);
            mu_sort = sortrows([mu_t(:,1), (1:d)'],1);
            for j=1:floor(1+3*lambda)
                mu_high_ind = mu_sort(d-j+1,2);
                w_next(mu_high_ind) = 1;
                cor_smaller = -1*ones(d, 1)*lambda > cor(:, mu_high_ind);
                [maximum, min_cor_ind] = max(mu_t.*cor_smaller);
                if maximum == 0
                    [~, min_cor_ind] = min(abs(cor(:, mu_high_ind)-max(lambda*ones(d, 1),ones(d,1))));
                end
                w_next(min_cor_ind) = 1;
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
        
        % Put deltas in different matrices for c and cov
        for j=1:d
            if w_delta(j) ~= 0
                s_w_delta_nonzero(j,w_nonzero_counter(j)) = diff(j);
                w_w_delta_nonzero(j,w_nonzero_counter(j)) = w_delta(j);
                w_nonzero_counter(j) = w_nonzero_counter(j) + 1;
            end
        end

        
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
