function simObj =  simulation_strategy(simObj, mu_in, c_in, M_in)
    simObj.reset();
    d = simObj.d;
    T = simObj.T;
    lambda = 0.2;
    cur_strat = ones(d,1)/d;
    for i=1:T
        if mod(i,floor(T/15))==0 %happens 15 times over simulation
            strategies = [get_strategies(d, mu_in),cur_strat]; %set of strategies
            mu = [mu_in, mu_in]; %confidence intervals, can/will be changed to estimator
            c = [c_in, c_in];
            res = simulate(strategies, 200, simObj.s_cur, simObj.w_cur, mu, c, simObj.eta, 20, M_in);
            max = -1;
            max_ind = -1;
            for j=1:size(strategies,2)
                if res(j,1) - lambda*res(j,2) > max
                    max = res(j,1) - lambda*res(j,2);
                    max_ind = j;
                end
            end
            cur_strat = strategies(:,max_ind)/sum(strategies(:,max_ind));
        end
        simObj.step(cur_strat);
    end
    
    
end



function result = simulate(strategies, T_remain, s_curr, w_prev, mu_conf, c_conf, eta, n_sim, M)
    num_strat = size(strategies,2);
    d=size(strategies,1);
    result = zeros(2, num_strat);
    parfor i=1:num_strat %simulation for each strategy
        strategy = strategies(:,i);
        strategy = strategy / sum(strategy);
        strategy_result = zeros(1, n_sim);
        for j=1:n_sim
            mu = rand(d,1).*abs(mu_conf(:,2) - mu_conf(:,1)) + min(mu_conf')';
            c = rand(d,1).*abs(c_conf(:,2) - c_conf(:,1)) + min(c_conf')';
            model_params = struct('mu', mu,'M',M,'c',c , 'eta', eta);
            sim_obj = MarketSimulator_copy(T_remain,s_curr,model_params, w_prev); %parameters are passed to modded Simulator
            sim_obj.step(strategy);
            for k=1:T_remain-1
                sim_obj.step(strategy);
            end
            strategy_result(j) = sim_obj.r_cur;
        end
        result(:,i) = [mean(strategy_result), var(strategy_result)]; %evaluation parameters
    end
    result = result';
end


function strat = get_strategies(d, mu) %returns all combinations of the shares with the highest mu
    num_shares = floor(log(d)+1);
    %s_list = sim_obj.s_cur;
    s_list = mu;
    strat = zeros(d,2^num_shares);
    max_indx = zeros(1, num_shares);
    for i=1:num_shares
        [~, argmax] = max(s_list);
        max_indx(i) = argmax;
        s_list(argmax) = -1;
    end
    for i=1:2^num_shares-1
        bin_str = dec2bin(i);
        while strlength(bin_str)<num_shares
            bin_str = "0" + bin_str;
        end
        for j=1:num_shares
            strat(max_indx(j), i) = str2double(extract(bin_str,j));
        end
    end
    strat(:,2^num_shares) = ones(d,1)/d;
end
