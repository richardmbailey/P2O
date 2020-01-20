function [capex_cost_timeseries, cost_rate] = capex_cost_timeseries(specific_plastic_flow, assest_cost, asset_capacity, asset_duration,...
    learning_rate, initial_rate, timeseries_type, formal_sorting_price, CAPEX_rate_timeseries, plastic_flow_box_total)

% ------------------------------------------
% Options preserved below, but all calculations performed in current model
% are of type 4 here, hence this is hard-coded
timeseries_type = 4;
% ------------------------------------------

assert(timeseries_type==1 || timeseries_type==2 || timeseries_type==3 || timeseries_type==4, 'Error in CAPEX time series type');

% -------------------------------------------------------------------------
if timeseries_type == 1
    
    plastic_flow_box_total = reshape(plastic_flow_box_total, numel(plastic_flow_box_total),1);
    capex_cost_timeseries = zeros(numel(plastic_flow_box_total),1);
    cost_rate = zeros(length(plastic_flow_box_total), 1);
    cost_rate(1) = initial_rate;
    capex_cost_timeseries(1) = cost_rate(1) * plastic_flow_box_total(1);
    
    for i = 2:numel(plastic_flow_box_total)
        
        cost_rate(i) = cost_rate(i-1) * (plastic_flow_box_total(i) / plastic_flow_box_total(i-1))^-learning_rate;
        
        if isnan(cost_rate(i)) || plastic_flow_box_total(i-1) == 0
            cost_rate(i) = cost_rate(i-1);
        end
        
        capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);
        
    end
    
    % ---------------------------------------------------------------------
elseif timeseries_type == 2
    
    plastic_flow_box_total = reshape(plastic_flow_box_total, numel(plastic_flow_box_total),1);
    capex_cost_timeseries = zeros(numel(plastic_flow_box_total),1);
    cost_rate = zeros(length(plastic_flow_box_total), 1);
    
    for i = 1:numel(plastic_flow_box_total)
        
        cost_rate(i) = CAPEX_rate_timeseries(i);
        capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);
        
    end
    % ---------------------------------------------------------------------
elseif timeseries_type == 3
    
    plastic_flow_box_total = reshape(plastic_flow_box_total, numel(plastic_flow_box_total),1);
    capex_cost_timeseries = zeros(numel(plastic_flow_box_total),1);
    cost_rate = zeros(length(plastic_flow_box_total), 1);
    cost_rate(1) = initial_rate(1);
    capex_cost_timeseries(1) = cost_rate(1) * specific_plastic_flow(1);
    
    for i = 2:numel(plastic_flow_box_total)
        
        cost_rate(i) = ( (cost_rate(i-1)-formal_sorting_price(i-1)) * (plastic_flow_box_total(i)/plastic_flow_box_total(i-1))^-learning_rate )    + formal_sorting_price(i);
        
        if isnan(cost_rate(i))             
            
            cost_rate(i) = cost_rate(i-1) + formal_sorting_price(i);
            capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);
            
        else
            capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);    
        end  
        

    end
       
    % ---------------------------------------------------------------------
elseif timeseries_type == 4
    
    plastic_flow_box_total = reshape(plastic_flow_box_total, numel(plastic_flow_box_total),1);
    capex_cost_timeseries = zeros(numel(plastic_flow_box_total),1);
    cost_rate = zeros(length(plastic_flow_box_total), 1);
    cost_rate(1) = initial_rate(1);
    capex_cost_timeseries(1) = cost_rate(1) * specific_plastic_flow(1);
    
    for i = 2:numel(plastic_flow_box_total)
        
        cost_rate(i) = cost_rate(i-1)...
                          * (   (plastic_flow_box_total(1) * 10 + sum(plastic_flow_box_total(1:i)))...
                              / (plastic_flow_box_total(1) * 10 + sum(plastic_flow_box_total(1:i-1)))...
                            )^-log2(1+learning_rate);
       
        if isnan(cost_rate(i)) || cost_rate(i)==0  
            
            cost_rate(i) = cost_rate(i-1);
            capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);
            
        else
            capex_cost_timeseries(i) = cost_rate(i) * specific_plastic_flow(i);    
        end  
        

    end
    
    
else
end

end

