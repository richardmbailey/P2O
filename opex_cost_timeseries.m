function opex_cost_timeseries = opex_cost_timeseries(plastic_flow, initial_rate, learning_rate)
    
    plastic_flow = reshape(plastic_flow, numel(plastic_flow),1);
    opex_cost_timeseries = zeros(numel(plastic_flow),1);
    cost_rate = zeros(length(plastic_flow), 1);
    cost_rate(1) = initial_rate(1);
    opex_cost_timeseries(1) = cost_rate(1) * plastic_flow(1);
    
    for i = 2:numel(plastic_flow)
        
        cost_rate(i) = cost_rate(i-1)...
                          * (   (plastic_flow(1) * 10 + sum(plastic_flow(1:i)))...
                              / (plastic_flow(1) * 10 + sum(plastic_flow(1:i-1)))...
                            )^-log2(1+learning_rate);

        if isnan(cost_rate(i)) || plastic_flow(i-1) == 0
            cost_rate(i) = cost_rate(i-1);
        end 
        opex_cost_timeseries(i) = cost_rate(i) * plastic_flow(i);
        
    end

end




        
        

        
