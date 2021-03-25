function [simObj, mu_t, c_t, cov, cor] = real_strategy(simObj) % mu only given to compare approximation to exact value.
    w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    simObj = simObj.reset(); % reset simulation environment
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
    for i = 1:simObj.T
        
        if i < min(simObj.T/2,40)
            simObj = simObj.step(w_const); % Strategy only temporary.
        else
            %maxS = simObj.s_cur == max(max(mu_t(:,1))); % Find maxS s.
            %w_t = zeros(simObj.d,1);
            %w_t(maxS) = 1;
            goodPairs = cor <= zeros(simObj.d, simObj.d); % returns logical matrix where (i,j)-th entry is 1 iff Cor(i,j) <= 0. 
            maxInd = ones(2,1);
            s_normed = simObj.s_cur - simObj.s0; % To get the relative growth.
            for j = 1:simObj.d
                for k = j:simObj.d
                    if goodPairs(j,k) 
                        if s_normed(maxInd(1)) + s_normed(maxInd(2)) < s_normed(j) + s_normed(k)
                            maxInd = [j,k];
                        end
                    end
                end
            end
            w_t = zeros(simObj.d,1);
            if min(s_normed(maxInd(1)), s_normed(maxInd(2))) < 0
                if s_normed(maxInd(1)) < 0 && s_normed(maxInd(2)) < 0
                    w_t(maxInd(1)) = 1/2;
                    w_t(maxInd(2)) = 1/2;
                end 
                
                if s_normed(maxInd(1)) < 0 
                    w_t(maxInd(2)) = 1;
                else 
                    w_t(maxInd(1)) = 1;
                end
            else
                w_t(maxInd(1)) = s_normed(maxInd(1));
                w_t(maxInd(2)) = s_normed(maxInd(2));
                w_t = w_t / norm(w_t,1);
                % Just a bunch of different ways how the previous w's
                % should be weight in for the new w.
                %simObj = simObj.step((1/3)*(w_t + simObj.w_cur + simObj.w_hist(:,simObj.t-1)));
                %simObj = simObj.step((1/6)*(3*w_t + 2*simObj.w_cur + simObj.w_hist(:,simObj.t-1)));
                %simObj = simObj.step((1/2)*(w_t + simObj.w_cur));
                %simObj = simObj.step(w_t);
                %simObj = simObj.step((w_t + simObj.w_cur + transpose(sum(transpose(simObj.w_hist(:,1:(i-1))))))/(i+1));
            end
            %simObj = simObj.step((1/2)*(w_t + simObj.w_cur));
            %simObj = simObj.step(((1/2)*w_t + (1/4)*simObj.w_cur + (1/4)*simObj.w_hist(:,simObj.t-1)));
            simObj = simObj.step(((4/7)*w_t + (2/7)*simObj.w_cur + (1/7)*simObj.w_hist(:,simObj.t-1)));
        end
        
        % Compute difference of weights.
        if simObj.t == 1
            w_delta = zeros(simObj.d, 1);
        else
            w_delta = simObj.w_cur - simObj.w_hist();
        end
        
            
        % Prepare data by looking at the differences.
        diff = log(simObj.s_cur) - log(simObj.s_hist(:,simObj.t));
        diff_hist(:,i) = diff;
        
        % Put deltas in different matrices
        for j=1:simObj.d
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
    end
end
