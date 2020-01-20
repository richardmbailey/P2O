

[n_rows_demand_table_original, n_cols_demand_table_original] = size(demand_table_original);
demand_table = zeros(n_rows_demand_table_original, n_cols_demand_table_original-1);

for i=1: n_rows_demand_table
    
    demand_table(i,:) = demand_table_original(i, 2:end) .* demand_table_MC_multipliers(i, mc);
    
end

n_years = n_cols_demand_table_original - 1;

population = demand_table(1, 1:n_years);
waste_per_capita = demand_table(2, 1:n_years);
waste_proportion = demand_table(3:5, 1:n_years);

plastic_demand = zeros(n_plastic_types, n_years);

for i=1:n_plastic_types
    plastic_demand(i,:)         = population .* zone_proportions(z) .* waste_per_capita .* waste_proportion(i,:);
    gros_plastic_demand(i,:, z) = population .* zone_proportions(z) .* waste_per_capita .* waste_proportion(i,:);
end

% Reduce and substitute

% (re-)initialize to zero
reduce_eliminate = zeros(n_plastic_types, n_years);
reduce_reuse = zeros(n_plastic_types, n_years);
reduce_new_delivery = zeros(n_plastic_types, n_years);
substitute_paper = zeros(n_plastic_types, n_years);
substitute_coated_paper = zeros(n_plastic_types, n_years);
substitute_compostables = zeros(n_plastic_types, n_years);

proportion_shift_multi_to_rigid =    demand_table(24, 1:end);
proportion_shift_multi_to_flexible = demand_table(25, 1:end);
proportion_shift_flexible_to_rigid = demand_table(26, 1:end);

for i=1:n_plastic_types
    reduce_eliminate(i,:) =        plastic_demand(i,:) .* demand_table(6 + (i - 1) * 6, :);
    reduce_reuse(i,:) =            plastic_demand(i,:) .* demand_table(7 + (i - 1) * 6, :);
    reduce_new_delivery(i,:) =     plastic_demand(i,:) .* demand_table(8 + (i - 1) * 6, :);
    substitute_paper(i,:) =        plastic_demand(i,:) .* demand_table(9 + (i - 1) * 6, :);
    substitute_coated_paper(i,:) = plastic_demand(i,:) .* demand_table(10 + (i - 1) * 6, :);
    substitute_compostables(i,:) = plastic_demand(i,:) .* demand_table(11 + (i - 1) * 6, :);
end

sum_R_and_S = reduce_eliminate + reduce_reuse + reduce_new_delivery + substitute_paper...
    + substitute_coated_paper + substitute_compostables;

% Shift from Multi to Rigid
shift_multi_to_rigid = (plastic_demand(2,:) - sum_R_and_S(2, :)) .* proportion_shift_multi_to_rigid;
% Shift from Multi to Flexible
shift_multi_to_flexible = (plastic_demand(2,:) - sum_R_and_S(2, :)) .* proportion_shift_multi_to_flexible;
% Shift from Flexible to Rigid
shift_flexible_to_rigid = (plastic_demand(3,:) - sum_R_and_S(3, :)) .* proportion_shift_flexible_to_rigid;

% Update waste calculations with R&S and shift calculations
plastic_demand(1,:) = plastic_demand(1,:) - sum_R_and_S(1, :) + shift_multi_to_rigid + shift_flexible_to_rigid;
plastic_demand(2,:) = plastic_demand(2,:) - sum_R_and_S(2, :) + shift_multi_to_flexible - shift_flexible_to_rigid;
plastic_demand(3,:) = plastic_demand(3,:) - sum_R_and_S(3, :) - shift_multi_to_rigid - shift_multi_to_flexible;

if baseline_run == 1
    
    production_rates = zeros(1, n_plastic_types);
    
    
    for j = plastic_types
        production_rates(1,j) = plastic_demand(j,1);
    end
    
elseif baseline_run == 0
    
    production_rates = plastic_demand';
    production_rates_table(:,:,z) = production_rates;
 
    for j = plastic_types
        
        del = interp1(1:duration,...
            plastic_demand(j,1:duration),...
            [0.5:24.5], 'linear', 'extrap');
        
        production_rates(:,j) = del';
        production_rates_table(:,j,z) = del';
        
    end
    
    
else
end

