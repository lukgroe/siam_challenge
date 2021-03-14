function [simObj, mu_t, c_t] = mixed_strategy(simObj, mu, c) % mu only given to compare approximation to exact value.
    w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    simObj = simObj.reset(); % reset simulation environment
    s_w_delta_zero = zeros(simObj.d, simObj.T); %delta s for delta w equals zero
    s_w_delta_nonzero = zeros(simObj.d, simObj.T); %delta s for delta w not equal zero
    w_w_delta_nonzero = zeros(simObj.d, simObj.T); %delta w for delta w not equal zero
    w_zero_counter = ones(simObj.d,1);
    w_nonzero_counter = ones(simObj.d,1);
    mu_t = zeros(simObj.d, 2); % [mu approximation, number of data used]
    c_t = zeros(simObj.d, 1);
    
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
        
        % Put deltas in different matrices
        for j=1:simObj.d
            if w_delta(j) == 0
                s_w_delta_zero(j,w_zero_counter(j)) = diff(j);
                w_zero_counter(j) = w_zero_counter(j) + 1;
            else
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
            if w_zero_counter(j)-1 ~= 0 && w_nonzero_counter(j)-1 ~= 0
                for k=1:max(w_zero_counter(j)-1, w_nonzero_counter(j)-1)
                    c_t(j) = c_t(j) + (s_w_delta_zero(j,mod(k, w_zero_counter(j)-1)+1)-s_w_delta_nonzero(j,mod(k, w_nonzero_counter(j)-1)+1))/w_w_delta_nonzero(j,mod(k, w_nonzero_counter(j)-1)+1);
                end
            end
        end
        
        % Compare to exact mu.
        norm(mu - mu_t(:,1))
        norm(c-c_t)
    end
end
