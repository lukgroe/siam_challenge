function simObj = parameter_estimation_strategy(simObj, mu) % mu only given to compare approximation to exact value.
    w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    simObj = simObj.reset(); % reset simulation environment
    
    mu_t = zeros(simObj.d, 2); % [mu approximation, number of data used]
    
    % Simulate.
    for i=1:simObj.T
                
        % Compute difference of weights.
        if simObj.t == 0
            w_delta = zeros(simObj.d, 1);
        else
            w_delta = simObj.w_cur - simObj.w_hist(:,simObj.t);
        end
        
        simObj = simObj.step(w_const); % Strategy only temporary.
        
        % Prepare data by looking at the differences.
        diff = log(simObj.s_cur) - log(simObj.s_hist(:,simObj.t));
        
        % Calculate current mu approximation from this round.
        for j = 1:(simObj.d)
            if w_delta(j) == 0 
                mu_t(j,1) = mu_t(j,2) * mu_t(j,1) + diff(j); 
                mu_t(j,2) = mu_t(j,2) + 1;
                mu_t(j,1) = mu_t(j,1) / mu_t(j,2);
            end
        end
        
        % Compare to exact mu.
        norm(mu - mu_t(:,1))
    end
end
