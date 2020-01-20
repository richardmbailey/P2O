function[output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
    compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
    Mass, MC_interaction_flows, MC_stdev_flows, MC_stdev_masses)


if baseline_run == 0
      
    output_interaction_flow = zeros(1+(duration * output_resolution), n_flows, n_plastic_types);
    output_mass = zeros(1+(duration * output_resolution), n_stocks, n_plastic_types);
    
    output_interaction_flow_stdev = zeros(1+(duration * output_resolution), n_flows, n_plastic_types);
    output_mass_stdev = zeros(1+(duration * output_resolution), n_stocks, n_plastic_types);
    
    % Create summary flow results version at specified resolution
    for p = 1:n_plastic_types
        
        for f = 1 : n_flows
            
            output_interaction_flow(1,:,p) = MC_interaction_flows(1,:,p); %add first point
        
            for i = 1:duration
                
                for j = 1 : output_resolution
                    
                    if i == 1
                        start_sum_index = 1;
                        end_sum_index = 1 + round((datapoints_per_year / output_resolution), 0) - 1;
                        
                        output_interaction_flow(1 + ((i-1) * output_resolution) + j, f, p) =...
                        sum(MC_interaction_flows(start_sum_index : end_sum_index, f, p,1)) * dt * output_resolution; %(1/(1-0.1));
                        
                        output_interaction_flow_stdev(1 + ((i-1) * output_resolution) + j, f, p) =...
                        mean(MC_stdev_flows(start_sum_index : end_sum_index, f, p));
                    
                    else
                        start_sum_index = 1 + ((i - 1) * datapoints_per_year) + (j - 1)...
                            * round((datapoints_per_year / output_resolution), 0);
                        end_sum_index = start_sum_index + round((datapoints_per_year / output_resolution), 0) - 1;
                        
                        output_interaction_flow(1 + ((i-1) * output_resolution) + j, f, p) =...
                        sum(MC_interaction_flows(start_sum_index : end_sum_index, f, p,1)) * dt * output_resolution;
                                          
                        output_interaction_flow_stdev(1 + ((i-1) * output_resolution) + j, f, p) =...
                        mean(MC_stdev_flows(start_sum_index : end_sum_index, f, p));
                    end                                     
                    % in units of Mt/yr

                end
            end
        end
        
        % Create Mass results version at specified resolution
        for f = 1 : n_stocks
            
            output_mass(1,:,p) = Mass(1,:,p,1); %add first point
            
            for i = 1:duration
                
                for j = 1 : output_resolution
                    
                    start_sum_index = 1 + ((i - 1) * datapoints_per_year) + (j - 1)...
                        * round((datapoints_per_year / output_resolution), 0);
                    end_sum_index = start_sum_index + round((datapoints_per_year / output_resolution), 0) - 1;
                    
                    output_mass(1 + ((i-1) * output_resolution) + j, f,p) =...
                        mean(Mass(start_sum_index : end_sum_index, f,p,1));
                    
                    output_mass_stdev(1 + ((i-1) * output_resolution) + j, f,p) =...
                        mean(MC_stdev_masses(start_sum_index : end_sum_index, f,p));
                    % in units of Mt
                    
                end
            end
        end
    end
    
end

end