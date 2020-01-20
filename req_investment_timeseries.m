function req_investment_timeseries = req_investment_timeseries(plastic_flow, initial_rate, learning_rate)
    
    plastic_flow = reshape(plastic_flow, numel(plastic_flow),1);
    req_investment_timeseries = zeros(numel(plastic_flow),1);
    cost_rate = zeros(length(plastic_flow), 1);
    cost_rate(1) = initial_rate(1);
    req_investment_timeseries(1) = cost_rate(1) * plastic_flow(1);
    
    for i = 2:numel(plastic_flow)
        
        cost_rate(i) = cost_rate(i-1) * (plastic_flow(i) / plastic_flow(i-1))^-learning_rate;

        if isnan(cost_rate(i)) || plastic_flow(i-1) == 0
            cost_rate(i) = cost_rate(i-1);
        end 
        req_investment_timeseries(i) = cost_rate(i) * plastic_flow(i);
        
    end   

end