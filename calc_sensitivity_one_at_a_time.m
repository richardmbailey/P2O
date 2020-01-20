function [sensitivity_x, sensitivity_y] = calc_sensitivity_one_at_a_time(...
            analysis_year, x_analysis_flow, y_analysis_flows, plastic_types, n_plastic_types,...
            MC_interaction_flows, n_MC_iterations, dt, x_flow_sensitivity_table, y_flow_sensitivity_table, output_mass)

sensitivity_x = zeros(n_plastic_types, n_MC_iterations);
sensitivity_y = zeros(n_plastic_types, n_MC_iterations);

n_y_analysis_flows = sum(y_analysis_flows > 0);

for k = plastic_types
    
        for s = 1: n_y_analysis_flows
            sensitivity_y(k,:) = sensitivity_y(k,:)...
                + reshape( y_flow_sensitivity_table(analysis_year, y_analysis_flows(s), :),1, n_MC_iterations) ;

        end
        
        n_flows = numel(x_analysis_flow);
        
        for i = 1:n_flows
            sensitivity_x(k,:) = sensitivity_x(k,:)...
            + reshape( x_flow_sensitivity_table(analysis_year, x_analysis_flow, :),1, n_MC_iterations) ;
        end
        
 
end

end