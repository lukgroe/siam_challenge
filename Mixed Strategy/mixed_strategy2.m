function [simObj, mu_t, c_t, cov, cor] = mixed_strategy2(simObj, lambda) % mu only given to compare approximation to exact value.
    %w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    w_next = zeros(simObj.d,1);
    w = max(1/floor(1/lambda),1/simObj.d);
    while sum(w_next)<1-w
        w_next(randi([1,simObj.d])) = w;
    end
    simObj = simObj.reset(); % reset simulation environment
    s_w_delta_nonzero = zeros(simObj.d, simObj.T); % Delta of s for delta w not equal zero
    w_w_delta_nonzero = zeros(simObj.d, simObj.T); % Delta of w for delta w not equal zero
    w_nonzero_counter = ones(simObj.d,1); % Counter for number of data points
    mu_t = zeros(simObj.d, 2); % [mu approximation, number of data points used]
    c_t = zeros(simObj.d, 1);
    cov = zeros(simObj.d, simObj.d);
    cov_counter = zeros(simObj.d, simObj.d);
    diff_hist = zeros(simObj.d, simObj.T);
    cor = zeros(simObj.d, simObj.d);
    
    % Simulate.
    for i=1:simObj.T
        
        if i<simObj.T/3 || mod(i,10) ~= 0
            simObj = simObj.step(w_next); 
        else
            
            c_t = zeros(simObj.d, 1); %currently unnecessary calculation of c
            for j=1:simObj.d
                for k=1:simObj.d % Correlation is being calculated
                    cor(j,k) = cov(j,k)/sqrt(cov(j,j)*cov(k,k));
                end
                if w_nonzero_counter(j)-1 ~= 0
                    for k=1:(w_nonzero_counter(j)-1)
                        c_t(j) = c_t(j) + (s_w_delta_nonzero(j,k)-mu_t(j,1))/(sign(w_w_delta_nonzero(j,k)*sqrt(abs(w_w_delta_nonzero(j,k))))); % Adding up the differences
                    end
                    c_t(j) = c_t(j)/w_nonzero_counter(j)-1; % Dividing by number of data points
                end
            end
            
            
            w_next = zeros(simObj.d,1);
            mu_sort = sortrows([mu_t(:,1), (1:simObj.d)'],1);
            for j=1:floor(1/lambda)
                mu_high_ind = mu_sort(simObj.d-j+1,2);
                w_next(mu_high_ind) = 1;
                [~, min_cor_ind] = min(cor(:, mu_high_ind));
                w_next(min_cor_ind) = 1;
            end
            w_next = w_next/sum(w_next);
            simObj = simObj.step(w_next);
        end
        
        % Compute difference of weights.
        if simObj.t == 1
            w_delta = zeros(simObj.d, 1);
        else
            w_delta = simObj.w_cur - simObj.w_hist(:,simObj.t-1);
        end
        
            
        % Prepare data by looking at the differences.
        diff = log(simObj.s_cur) - log(simObj.s_hist(:,simObj.t));
        diff_hist(:,i) = diff;
        
        % Put deltas in different matrices for c and cov
        for j=1:simObj.d
            if w_delta(j) ~= 0
                s_w_delta_nonzero(j,w_nonzero_counter(j)) = diff(j);
                w_w_delta_nonzero(j,w_nonzero_counter(j)) = w_delta(j);
                w_nonzero_counter(j) = w_nonzero_counter(j) + 1;
            end
        end

        
        % Calculate current mu approximation from this round.
        for j = 1:(simObj.d)
            if w_delta(j) == 0
                mu_t(j,1) = mu_t(j,2) * mu_t(j,1) + diff(j); 
                mu_t(j,2) = mu_t(j,2) + 1;
                mu_t(j,1) = mu_t(j,1) / mu_t(j,2);
            end
        end
        
        for j=1:simObj.d
            for k=1:simObj.d
                if w_delta(j) == 0 && w_delta(k) == 0
                    cov(j,k) = (cov(j,k)*max(1, cov_counter(j,k)-1) + (diff(j)-mu_t(j,1))*(diff(k)-mu_t(k,1)))/max(1,cov_counter(j,k));
                    cov_counter(j,k) = cov_counter(j,k)+1;
                end
            end
        end
    end
end
