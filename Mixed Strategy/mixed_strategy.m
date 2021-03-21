function [simObj, mu_t, c_t, cov, cor] = mixed_strategy(simObj, mu, c) % mu only given to compare approximation to exact value.
    w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    simObj = simObj.reset(); % reset simulation environment
    %s_w_delta_zero = zeros(simObj.d, simObj.T); %delta s for delta w equals zero
    s_w_delta_nonzero = zeros(simObj.d, simObj.T); %delta s for delta w not equal zero
    w_w_delta_nonzero = zeros(simObj.d, simObj.T); %delta w for delta w not equal zero
    w_zero_counter = ones(simObj.d,1);
    w_nonzero_counter = ones(simObj.d,1);
    mu_t = zeros(simObj.d, 2); % [mu approximation, number of data used]
    c_t = zeros(simObj.d, 1);
    cov = zeros(simObj.d, simObj.d);
    cov_counter = zeros(simObj.d, simObj.d);
    diff_hist = zeros(simObj.d, simObj.T);
    cor = zeros(simObj.d, simObj.d);
    
    % Simulate.
    for i=1:simObj.T
        
        if i<simObj.T/2
            simObj = simObj.step(w_const); % Strategy only temporary.
        else
            w_random = poissrnd(10,simObj.d, 1) + ones(simObj.d,1);
            w_random = w_random./sum(w_random);
            simObj = simObj.step(w_random); %random temporary strategy
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
        
        % Put deltas in different matrices
        for j=1:simObj.d
            %if w_delta(j) == 0
            %    s_w_delta_zero(j,w_zero_counter(j)) = diff(j);
            %    w_zero_counter(j) = w_zero_counter(j) + 1;
            %else
            if w_delta(j) ~= 0
                s_w_delta_nonzero(j,w_nonzero_counter(j)) = diff(j);
                w_w_delta_nonzero(j,w_nonzero_counter(j)) = w_delta(j);
                w_nonzero_counter(j) = w_nonzero_counter(j) + 1;
            end
        end

        
        % Calculate current mu and c approximation from this round.
        for j = 1:(simObj.d)
            if w_delta(j) == 0
                mu_t(j,1) = mu_t(j,2) * mu_t(j,1) + diff(j); 
                mu_t(j,2) = mu_t(j,2) + 1;
                mu_t(j,1) = mu_t(j,1) / mu_t(j,2);
            end
            if w_nonzero_counter(j)-1 ~= 0
                for k=1:(w_nonzero_counter(j)-1)
                    c_t(j) = c_t(j) + (s_w_delta_nonzero(j,k)-mu_t(j,1))/(sign(w_w_delta_nonzero(j,k)*sqrt(abs(w_w_delta_nonzero(j,k)))));
                    %c_t(j) = c_t(j) + (s_w_delta_zero(j,mod(k, w_zero_counter(j)-1)+1)-s_w_delta_nonzero(j,mod(k, w_nonzero_counter(j)-1)+1))/(sign(w_w_delta_nonzero(j,mod(k, w_nonzero_counter(j)-1)+1))*sqrt(abs(w_w_delta_nonzero(j,mod(k, w_nonzero_counter(j)-1)+1))));
                end
                c_t(j) = c_t(j)/max(w_zero_counter(j)-1, w_nonzero_counter(j)-1);
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
        
        for j=1:simObj.d
            for k=1:simObj.d
                cor(j,k) = cov(j,k)/sqrt(cov(j,j)*cov(k,k));
            end
        end
        
        % Compare to exact mu.
        norm(mu - mu_t(:,1))
        norm(c-c_t)
    end
end
