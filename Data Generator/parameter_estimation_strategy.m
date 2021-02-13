function simObj = parameter_estimation_strategy(simObj)
    w_const = ones(simObj.d,1)/simObj.d; % Only really important thing is that it remains constant.
    simObj = simObj.reset(); % reset simulation environment
    for i=1:simObj.T
       simObj = simObj.step(w_const);
    end
end
