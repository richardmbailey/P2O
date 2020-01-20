function abs_flow = calc_absolute_flow(...
    equation_or_timeseries, interaction_number, flow_proportions,...
    plug, a, b, c, d, t, i, interpolated_timeseries_value,...
    function_type, plastic_type,...
    test_box, call_origin, call, zone_proportions, z, k)

if equation_or_timeseries(interaction_number,1) == 1  % absolute, equation
    
    if plug(interaction_number,1) == 1
        abs_flow = flow_proportions(interaction_number);
    elseif plug(interaction_number,1) == 0
        
        abs_flow = ...
            calculate_flow_rate(a(interaction_number,plastic_type),...
            b(interaction_number,plastic_type),c(interaction_number,plastic_type),...
            d(interaction_number,plastic_type),...
            t, function_type(interaction_number,plastic_type));
    else
    end
    
elseif equation_or_timeseries(interaction_number,1) == 2  % absolute, timeseries
    
    if plug(interaction_number,1) == 1
        abs_flow = flow_proportions(interaction_number);
    elseif plug(interaction_number,1) == 0
        
        abs_flow = interpolated_timeseries_value(interaction_number);
        
    else
    end
       
    if i==test_box && call_origin==call
        
        disp('------------------------');
        disp(['Box ', num2str(i)]);
        disp(['t: ', num2str(t)]);
        disp(['Inflow (1), Outflow (2)', num2str(k)]);
        disp(['equation_or_timeseries(1,2)', num2str(equation_or_timeseries(interaction_number,1))]);
        disp(['interaction number: ', num2str(interaction_number)]);
        disp(['abs_flow: ', num2str(abs_flow)]);
     
    end
  
else
end

abs_flow = abs_flow * zone_proportions(z);

end
