function simObj = random_strategy(simObj)
    simObj = simObj.reset();
    sub_strat = rand;
    if (sub_strat > .5)
        for i=1:simObj.T
            weights = zeros(simObj.d, 1);
            for j=1:simObj.d
                weights(j) = rand*(1-sum(weights));
            end
            weights = weights(randperm(simObj.d));
           simObj = simObj.step(weights);
        end
    else
        for i=1:simObj.T
            weights = zeros(simObj.d, 1);
            for j=1:simObj.d
                weights(j) = rand;
            end
            weights = weights ./sum(weights);
           simObj = simObj.step(weights);
        end
    end
end