function CAGR_timeseries = CAGR_timeseries(plastic_flow, initial_rate, CAGR_rate)
    
    plastic_flow = reshape(plastic_flow, numel(plastic_flow),1);
    rate = zeros(numel(plastic_flow), 1);
    rate(1) = initial_rate;
    
    for i=2:numel(plastic_flow)
        rate(i) = rate(i-1) * (1 + CAGR_rate);
    end
    
    rate(isnan(rate)) = 0;
    CAGR_timeseries = rate .* plastic_flow;

end