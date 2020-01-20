function [flow_proportions_list] = calc_flow_proportions3(n_flows, n_stocks, plug,...
    enforced_proportion, plastic_type, n_outflows, out_interaction_number,...
    intervention_t_start, intervention_t_end,...
    relative_absolute,a,b,c,d,t_real, function_type, production_interaction, mc,...
    interpolated_timeseries_value, local_M, processing_rate, equation_or_timeseries,...
    proportion_to_zone)


global box_active;

flow_proportions_list = zeros(n_flows,1);

for i = 1: n_stocks
   
    n_flows = n_outflows;
    flows_out = out_interaction_number(i, 1:n_outflows(i));
    
    absolute_values = find(relative_absolute(flows_out) == 2);
    a(flows_out(absolute_values), plastic_type) = 0;
    b(flows_out(absolute_values), plastic_type) = 0;
    c(flows_out(absolute_values), plastic_type) = 0;
    d(flows_out(absolute_values), plastic_type) = 0;
    
    %check if any of these are plug numbers
    n_plugs = sum(plug(flows_out, plastic_type));
    plug_flows = [];
    
    %if plugs present...
    if n_plugs > 0
        
        %identify the plug flows
        p_indices = plug(flows_out, plastic_type) == 1;
        plug_flows = flows_out(p_indices);
        
        %find and count the non-plug flows
        non_p_indices = plug(flows_out, plastic_type) == 0;
        n_non_plugs = sum(non_p_indices);
        non_plug_flows = flows_out(non_p_indices);
    end
    
    %if there are no plug values for this box...
    if n_plugs == 0
        non_plug_flows = flows_out;
        n_non_plugs = length(flows_out);
    end
    
    sum_plugs = 0;
    sum_non_plugs = 0;
    total_proportions = 0;
    non_plug_flow = zeros(n_flows(i), 1);
    
    % Check if mixed rel/abs/plugs
    n_relative_values = sum(relative_absolute(flows_out) == 1);
    n_absolute_values = sum(relative_absolute(flows_out) == 2);
    
    if n_relative_values>0 && n_absolute_values>0
        mixed = true();
    else
        mixed = false();
    end
    
    if ~mixed %not mixed, all same type of flow (rel or abs)
        
        %calculate the non-plug flows
        for k=1:n_non_plugs
            if t_real < intervention_t_start(non_plug_flows(k),plastic_type)
                t = 1;
            elseif t_real > intervention_t_end(non_plug_flows(k),plastic_type)
                t = intervention_t_end(non_plug_flows(k),plastic_type)...
                    - intervention_t_start(non_plug_flows(k),plastic_type) + 1;
            else
                t = t_real - intervention_t_start(non_plug_flows(k),plastic_type) + 1;
            end
            
            if equation_or_timeseries(non_plug_flows(k)) == 1 % eqn
                non_plug_flow(k) = ...
                    calculate_flow_rate(a(non_plug_flows(k),plastic_type),...
                    b(non_plug_flows(k),plastic_type),c(non_plug_flows(k),plastic_type),...
                    d(non_plug_flows(k),plastic_type),...
                    t,function_type(non_plug_flows(k),plastic_type));
            elseif equation_or_timeseries(non_plug_flows(k)) == 2 % timeseries
                non_plug_flow(k) = interpolated_timeseries_value( non_plug_flows(k) );
            else
            end
            
            flow_proportions_list(non_plug_flows(k)) = non_plug_flow(k);
            sum_non_plugs = sum_non_plugs + non_plug_flow(k);
        end
        
        %calculate the plug flows
        for k=1:n_plugs
            
            if relative_absolute( plug_flows(1), plastic_type  ) == 1
                reference = 1;
            elseif relative_absolute( plug_flows(1), plastic_type  ) == 2
                reference = local_M(i) * processing_rate(plug_flows(k), plastic_type);
            else
            end
            
            plug_flow = ...
                max(0, (reference - sum_non_plugs) * enforced_proportion(plug_flows(k),plastic_type));
            flow_proportions_list(plug_flows(k)) = plug_flow; %store here as either rel or abs
            sum_plugs = sum_plugs + plug_flow; %sum here in rel format
        end
        
    elseif mixed % flows are a mixture of rel/abs/plug
        
        sum_non_plugs_absolute = 0;
        sum_non_plugs_relative = 0;
        
        %calculate the non-plug absolute flows
        for k=1:n_non_plugs
            if t_real < intervention_t_start(non_plug_flows(k),plastic_type)
                t = 1;
            elseif t_real > intervention_t_end(non_plug_flows(k),plastic_type)
                t = intervention_t_end(non_plug_flows(k),plastic_type)...
                    - intervention_t_start(non_plug_flows(k),plastic_type) + 1;
            else
                t = t_real - intervention_t_start(non_plug_flows(k),plastic_type) + 1;
            end
            
            if relative_absolute( non_plug_flows(k), plastic_type  ) == 2  && plug(non_plug_flows(k)) == 0
                if equation_or_timeseries(non_plug_flows(k)) == 1 % eqn
                    non_plug_flow(k) = ...
                        calculate_flow_rate(a(non_plug_flows(k),plastic_type),...
                        b(non_plug_flows(k),plastic_type),c(non_plug_flows(k),plastic_type),...
                        d(non_plug_flows(k),plastic_type),...
                        t,function_type(non_plug_flows(k),plastic_type));
                elseif equation_or_timeseries(non_plug_flows(k)) == 2 % timeseries
                    non_plug_flow(k) = interpolated_timeseries_value( non_plug_flows(k) ) * proportion_to_zone;
                else
                end
                
                sum_non_plugs_absolute = sum_non_plugs_absolute + non_plug_flow(k);
                flow_proportions_list(non_plug_flows(k)) = non_plug_flow(k);
                
            end
        end
        
        residual_mass = local_M(i); 

        
        %calculate sum of non-plug relative flows
        for k=1:n_flows(i)
            if t_real < intervention_t_start(flows_out(k),plastic_type)
                t = 1;
            elseif t_real > intervention_t_end(flows_out(k),plastic_type)
                t = intervention_t_end(flows_out(k),plastic_type)...
                    - intervention_t_start(flows_out(k),plastic_type) + 1;
            else
                t = t_real - intervention_t_start(flows_out(k),plastic_type) + 1;
            end

            if relative_absolute( flows_out(k), plastic_type  ) == 1 && ...
                    plug( flows_out(k), plastic_type  ) == 0 %rel & not plug
                if equation_or_timeseries(flows_out(k)) == 1 % eqn
                    non_plug_flow(k) = ...
                        calculate_flow_rate(a(flows_out(k),plastic_type),...
                        b(flows_out(k),plastic_type),c(flows_out(k),plastic_type),...
                        d(flows_out(k),plastic_type),...
                        t,function_type(flows_out(k),plastic_type));
                elseif equation_or_timeseries(flows_out(k)) == 2 % timeseries
                    % save rel and absolute seperately then combine at the end
                    %                disp('Timeseries')
                    non_plug_flow(k) = interpolated_timeseries_value( flows_out(k) );
                else
                end
                
                % Sum relative flows
                sum_non_plugs_relative = sum_non_plugs_relative + non_plug_flow(k);
            end
        end
        
        % Calculate all relative flows (plug and non-plug) - all apart from the abs flow
        for k=1:n_flows(i)
            
            %if non-plug absolute flow
            if plug(flows_out(k), plastic_type ) == 0 ...
                    && relative_absolute( flows_out(k), plastic_type  ) == 2
                flow_proportions_list(flows_out(k)) = sum_non_plugs_absolute; %assumes just one abs flow       
            end
            
            %if non-plug relative
            if plug(flows_out(k), plastic_type ) == 0 ...
                    && relative_absolute( flows_out(k), plastic_type  ) == 1
                if local_M(i)>0
                    flow_proportions_list(flows_out(k)) = (non_plug_flow(k) * residual_mass) / local_M(i);
                else
                    flow_proportions_list(flows_out(k)) = 0;
                end
                
            end
            
            %if plug relative
            if plug(flows_out(k), plastic_type ) == 1 ...
                    && relative_absolute( flows_out(k), plastic_type  ) == 1
                
                if local_M(i)>0
                    flow_proportions_list(flows_out(k)) = ...
                        ((1 - sum_non_plugs_relative) * residual_mass) / local_M(i);
                else
                    flow_proportions_list(flows_out(k)) = 0;
                end
                
            end
        end
      
    else
    end %if ~mixed
    
    
end

end
