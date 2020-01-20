function [in_flow, out_flow, interaction_flow_at_t] = calc_flows(M,t,tspan,n_datapoints,dt,scaling,...
    n_stocks,n_flows,plastic_type,stock_stock_interactions,...
    intervention_t_start, intervention_t_end, equation_or_timeseries,arrow_timeseries,...
    relative_absolute,a,b,c,d,box_total_capacity,box_flow_capacity,box_flow_capacity_t_start, box_flow_capacity_t_end,...
    box_capacity_t_start, box_capacity_t_end,...
    production_interaction, production_rates, imports_interaction, exports_interaction, imports_rates, waste_generated_box_number,...
    plug, function_type, enforced_proportion, processing_rate,max_annual_flow_rate,...
    flow_flow_interactions, in_interaction_number, out_interaction_number,...
    n_inflows, n_outflows, in_index, in_source, out_source, out_index, mc, call_origin,...
    zone_proportions, z, finite_sinks, micro, macro, baseline_run, baseline_individual_flow)


% -------- Diagnostics ------------
test_box = -13; % set to box number to get diagnostic readouts as this function runs
test_box2 = -13;
call = 2; %for call_origin number
% ---------------------------------

%Note: M here is a vector of masses for each box for one time point

global box_active;
global switched_on;
global switched_off;
global switched_on_times;
global switched_off_times;
global year_counter;

in_flow_component = zeros(n_stocks,n_stocks);
out_flow_component = zeros(n_stocks,n_stocks);
production_rate = 0;


if t >= max(tspan)
    production_rate = production_rates(end, plastic_type);
elseif t~= max(tspan)
    
    % retained in case production_rate input parameter needed
    %{
    tx = find(tspan > t,1);
    first_t_index = max(tx - 1, 1);
    second_t_index = min(tx + 0, n_datapoints);
    
    production_rate = interp1([tspan(first_t_index), tspan(second_t_index)],...
        [production_rates(first_t_index, plastic_type), production_rates(second_t_index, plastic_type)],...
        t, 'linear');

    production_rate = interp1([tspan(first_t_index), tspan(second_t_index)],...
        [production_rates(first_t_index, plastic_type), production_rates(second_t_index, plastic_type)],...
        t, 'linear');
    %}
    
    [temp1, ~] = size(production_rates);
    
    if temp1 == 1
        production_rate = production_rates(1, plastic_type);
    elseif temp1 ~= 1
        %{
    int_time = [0.5:24.5]';
    production_rate = interp1(int_time,...
         production_rates(:, plastic_type),...
         t, 'linear', 'extrap');
        %}
        rowref = fix(t) + 1;
        production_rate = production_rates(rowref,plastic_type);
        
    else
    end
    
else
end


current_year_number = fix(t) + 1;

in_flow = zeros(n_stocks,1);
out_flow = zeros(n_stocks,1);
interaction_flow_at_t = zeros(n_flows,1); %(number of interactions can't exceed n_stocks for each stock/box)

%interpolate time series inputs
interpolated_timeseries_value = interpolate_arrow_at_t(equation_or_timeseries, plastic_type, n_flows,...
    tspan, arrow_timeseries, t, imports_interaction, exports_interaction);


flow_proportions = calc_flow_proportions3(n_flows, n_stocks, plug, enforced_proportion, plastic_type,...
    n_outflows, out_interaction_number, intervention_t_start, intervention_t_end,...
    relative_absolute, a,b,c,d,t, function_type, production_interaction, mc,...
    interpolated_timeseries_value, M, processing_rate, equation_or_timeseries, zone_proportions(z));

if macro; first_box = 2; end
if micro; first_box = 1; end

for i = first_box: n_stocks
    
    for k=1:2 %k=1 for in_flows, k=2 for out_flows
        
        if k==1; n_flows = n_inflows; end
        if k==2; n_flows = n_outflows; end
        
        for j = 1:n_flows(i)
            
            if k==1
                interaction_source = in_index(i,j);
                interaction_number = in_interaction_number(i,j);
            elseif k==2
                outflow_target = out_index(i,j);
                interaction_number = out_interaction_number(i,j);
                
            else
            end
            
            recycling_inflows = [];
            
            if interaction_number == production_interaction 
                %special case of first input to system, if production (demand) calculated internally is being used
                
                recycling_inputs = 0;
                
                % Find all recycled flows going in to waste production box
                [x,~] = ind2sub(size(out_index), find(stock_stock_interactions(:,waste_generated_box_number)>1));
                
                recycling_inflows = stock_stock_interactions(x,waste_generated_box_number);
                del = recycling_inflows == production_interaction;
                recycling_inflows(del) = []; %remove virgin production from list
                x(del) = [];
                
                %recycling_inflows
                for r = 1:length(recycling_inflows)
                    
                    recycling_inputs = recycling_inputs + ...
                        M(x(r)) * flow_proportions(recycling_inflows(r))...
                        * processing_rate(x(r),plastic_type) ;
                    
                end
                
                if k==1 % inflow
                    in_flow_component(i,j) = production_rate;% - recycling_inputs;
                    
                elseif k==2 % outflow
                    out_flow_component(i,j) = 0;
                    interaction_flow_at_t(interaction_number) = production_rate;% - recycling_inputs;
                else
                end
                
            elseif interaction_number ~= production_interaction %a normal flow between boxes
                
                if k==1 % inflow
                    
                    if box_active(i) == 1
                        
                        if relative_absolute(interaction_number) == 1 %relative
                            
                            if M(i) < box_total_capacity(current_year_number,i, plastic_type)
                                
                                if equation_or_timeseries(interaction_number,1) == 1    % rel equation
                                    in_flow_component(i,j) = min(...
                                        M(interaction_source) * flow_proportions(interaction_number)...
                                        * processing_rate(interaction_number,plastic_type),...
                                        max_annual_flow_rate(interaction_number) );
                                    
                                    if i==test_box && call_origin==call
                                        i
                                        disp('Inflow - equation - relative');
                                        fprintf('Plastic type: %d\n', plastic_type);
                                        interaction_number
                                        interaction_source
                                        flow_proportions(interaction_number)
                                        t
                                        M(i)
                                        box_total_capacity(current_year_number,i, plastic_type)
                                    end
                                    
                                elseif equation_or_timeseries(interaction_number,1) == 2    %  rel timeseries
                                    
                                    in_flow_component(i,j) = min(...
                                        M(interaction_source) * interpolated_timeseries_value(interaction_number)...
                                        * processing_rate(interaction_number,plastic_type),...
                                        max_annual_flow_rate(interaction_number) );
                                    
                                    
                                    if i==test_box && call_origin==call
                                        i
                                        disp('Inflow - timeseries - relative');
                                        fprintf('Plastic type: %d\n', plastic_type);
                                        interaction_number
                                        interaction_source
                                        interpolated_timeseries_value(interaction_number)
                                        t
                                    end
                                    
                                else
                                end
                                
                            elseif M(i) > box_total_capacity(current_year_number,i, plastic_type)
                                
                                in_flow_component(i,j) = 0;
                                
                            else
                            end
                            
                        elseif relative_absolute(interaction_number) == 2 %absolute
                            
                            %if capacity exists in the current box for inflow
                            if M(i) < box_total_capacity(current_year_number,i, plastic_type)
                                
                                testA = calc_absolute_flow(...
                                    equation_or_timeseries, interaction_number, flow_proportions,...
                                    plug, a, b, c, d, t, i, interpolated_timeseries_value,...
                                    function_type, plastic_type, test_box, call_origin, call, zone_proportions, z, k);
                                testB = max_annual_flow_rate(interaction_number);
                                testC = (M(interaction_source) * 1 * processing_rate(interaction_source,plastic_type));
                                
                                in_flow_component(i,j) = max(0, min([testA,testB,testC]) );
                                % defined absolute flow can't exceed max possible given mass and flow rate
                                
                            else
                                in_flow_component(i,j) = 0;
                            end
                            
                            if interaction_number == imports_interaction % over-write if imports
                                
                                if equation_or_timeseries(interaction_number,1) == 1
                                    interaction_flow_at_t(interaction_number,1) = calculate_flow_rate(a(interaction_number,plastic_type),...
                                        b(interaction_number,plastic_type),c(interaction_number,plastic_type),...
                                        d(interaction_number,plastic_type),...
                                        t,function_type(interaction_number,plastic_type));
                                elseif equation_or_timeseries(interaction_number,1) == 2
                                    interaction_flow_at_t(interaction_number,1) = calc_absolute_flow(...
                                        equation_or_timeseries, interaction_number, flow_proportions,...
                                        plug, a, b, c, d, t, i, interpolated_timeseries_value,...
                                        function_type, plastic_type, test_box, call_origin, call, zone_proportions, z, k);
                                else
                                end
                                in_flow_component(i,j) = interaction_flow_at_t(interaction_number,1);
                                
                            end
                            
                        else
                        end
                        
                    elseif box_active(i) == 0
                        in_flow_component(i,j) = 0;
                    else
                    end
                    
                    
                elseif k==2 % outflow
                    
                    if box_active(outflow_target) == 1 %true()
                        
                        if relative_absolute(interaction_number) == 1 %relative
                            
                            %if there is still room in the target box for inflow
                            if M(outflow_target) < box_total_capacity(current_year_number, outflow_target, plastic_type)
                                if equation_or_timeseries(interaction_number,1) == 1,...
                                        out_flow_component(i,j) = min(...
                                        M(i) * flow_proportions(interaction_number)...
                                        * processing_rate(interaction_number,plastic_type),...
                                        max_annual_flow_rate(interaction_number) );
                                    if i==test_box && call_origin==call
                                        i
                                        disp('Outflow - equation - relative');
                                        fprintf('Plastic type: %d\n', plastic_type);
                                        interaction_number
                                        M(i)
                                        flow_proportions(interaction_number)
                                        out_flow_component(i,j)
                                        t
                                    end
                                elseif equation_or_timeseries(interaction_number,1) == 2
                                    
                                    out_flow_component(i,j) = min(...
                                        M(i) * interpolated_timeseries_value(interaction_number)...
                                        * processing_rate(interaction_number,plastic_type),...
                                        max_annual_flow_rate(interaction_number) );
                                    if i==test_box && call_origin==call
                                        
                                        disp('Outflow - timeseries - relative');
                                        fprintf('Plastic type: %d\n', plastic_type);
                                        i
                                        box_active(i)
                                        call_origin
                                        interaction_number
                                        M(i)
                                        interpolated_timeseries_value(interaction_number)
                                        out_flow_component(i,j)
                                        t
                                    end
                                else
                                end
                                
                                
                            elseif M(outflow_target) > box_total_capacity(current_year_number, outflow_target, plastic_type)
                                
                                out_flow_component(i,j) = 0;
                                
                            else
                            end
                            interaction_flow_at_t(interaction_number,1) = out_flow_component(i,j);
                            
                        elseif relative_absolute(interaction_number) == 2 %absolute
                            
                            %if capacity in the target box for inflow
                            if M(outflow_target) < box_total_capacity(current_year_number, outflow_target, plastic_type)
                                
                                testA = calc_absolute_flow(...
                                    equation_or_timeseries, interaction_number, flow_proportions,...
                                    plug, a, b, c, d, t, i, interpolated_timeseries_value,...
                                    function_type, plastic_type, test_box, call_origin, call, zone_proportions, z, k);
                                testB = max_annual_flow_rate(interaction_number);
                                testC = (M(i) * 1 * processing_rate(i,plastic_type));
                                
                                out_flow_component(i,j) = max( 0, min([testA,testB,testC]) );
                                % defined absolute flow can't exceed max possible given mass and flow rate
                                
                                
                            elseif M(outflow_target) > box_total_capacity(current_year_number, outflow_target, plastic_type)
                                out_flow_component(i,j) = 0;
                                
                            else
                            end
                            interaction_flow_at_t(interaction_number,1) = out_flow_component(i,j);
                            if i==test_box && call_origin==call
                                
                                disp('------Outflow------');
                                disp(i);
                                fprintf('Plastic type: %d\n', plastic_type);
                                disp(testA);
                                disp(testB);
                                disp(testC);
                                disp(out_flow_component(i,j));
                                disp(t);
                            end
                            
                            if interaction_number == exports_interaction || interaction_number == imports_interaction
                                % over-write value if this is either the IMPORTS or EXPORTS flow
                                
                                if equation_or_timeseries(interaction_number,1) == 1
                                    interaction_flow_at_t(interaction_number,1) = calculate_flow_rate(a(interaction_number,plastic_type),...
                                        b(interaction_number,plastic_type),c(interaction_number,plastic_type),...
                                        d(interaction_number,plastic_type),...
                                        t,function_type(interaction_number,plastic_type));
                                elseif equation_or_timeseries(interaction_number,1) == 2
                                    
                                    % Set flow equal to the annual inputted time series value
                                    interaction_flow_at_t(interaction_number,1) = interpolated_timeseries_value(interaction_number) * zone_proportions(z);
                                else
                                end
                                
                                if interaction_number == imports_interaction
                                    out_flow_component(i,j) = 0;
                                end
                            end
                            
                        else
                        end
                        
                    elseif box_active(outflow_target)==0
                        out_flow_component(i,j) = 0;
                        interaction_flow_at_t(interaction_number,1) = 0;
                    else
                    end
                    
                else
                end
                
            end
        end
        
        
        
    end
    
    
end



if macro
    %% Special case: force flows A1, A2, B1, B2 equal to inputs, independent of above conditions
    
    % So that MC sampling doesn't result in incorrect total flows from Box B
    prescribed_1 = interpolated_timeseries_value(1) * zone_proportions(z);
    prescribed_2 = interpolated_timeseries_value(2) * zone_proportions(z);
    prescribed_3 = interpolated_timeseries_value(3) * zone_proportions(z);
    prescribed_4 = interpolated_timeseries_value(4) * zone_proportions(z);
    
    interaction_flow_at_t(1) = prescribed_1;
    interaction_flow_at_t(2) = prescribed_2;
    proportion_3 = prescribed_3 / (prescribed_3 + prescribed_4);
    proportion_4 = 1 - proportion_3;
    interaction_flow_at_t(3) = prescribed_1 * proportion_3;
    interaction_flow_at_t(4) = prescribed_1 * proportion_4;
    
    out_flow_component(1,1) = interaction_flow_at_t(1);
    out_flow_component(1,2) = interaction_flow_at_t(2);
    out_flow_component(2,1) = interaction_flow_at_t(3);
    out_flow_component(2,2) = interaction_flow_at_t(4);
    
    in_flow_component(2,1) = interaction_flow_at_t(1);
    in_flow_component(17,1) = interaction_flow_at_t(2);
    in_flow_component(3,1) = interaction_flow_at_t(3);
    in_flow_component(4,1) = interaction_flow_at_t(4);
    
    
    %% Special case: flows R2 (33), V1 (36), V2 (37), V3 (38)
    % define variables
    R2_original = interaction_flow_at_t(33,1);
    V1_original = interaction_flow_at_t(36,1);
    V2_original = interaction_flow_at_t(37,1);
    V3_original = interaction_flow_at_t(38,1);
    Q2_original = interaction_flow_at_t(30,1);
    T1_original = interaction_flow_at_t(34,1);
    
    proportion_to_V2 = interpolated_timeseries_value(37);
    proportion_to_V3 = interpolated_timeseries_value(38);
    proportion_to_T1 = interpolated_timeseries_value(34);
    
    % re-calculate flows
    if V1_original >= R2_original
        V1_new = R2_original;
        V2_new = 0;
        V3_new = 0;
        R2_new = 0;
    else
        V1_new = V1_original;
        V2_new = (R2_original - V1_original ) * proportion_to_V2;
        V3_new = (R2_original - V1_original ) * proportion_to_V3;
        R2_new = R2_original - V1_original - V2_new - V3_new;
    end
    
    T1_new = Q2_original * proportion_to_T1;
    
    % update existing interaction_flow arrays
    interaction_flow_at_t(33) = R2_new; % R2
    interaction_flow_at_t(34) = T1_new; % T1
    interaction_flow_at_t(36) = V1_new; % V1
    interaction_flow_at_t(37) = V2_new; % V2
    interaction_flow_at_t(38) = V3_new; % V3
    
    
    % Flow R2 (33)
    out_flow_component(18,2) = R2_original; % Box R (18), flow R2 (flow 33) is going to Box R (22), which is 2nd on list in out_index for Box 18
    in_flow_component(22,1) = R2_new; % Box V (22), flow R2 (flow 33) is coming from Box R (18), which is 1st on list in in_index for Box 22 OK
    
    % Flow V1 (36)
    out_flow_component(22,1) = 0; % Box V (22), flow V1 (flow 36) is going to Box D (4), which is 1st on list in out_index for Box 22  OK
    in_flow_component(4,2) = V1_new; % Box D (4), flow V1 (flow 36) is coming from Box V (22), which is 2nd on list in in_index for Box 4  OK
    
    % Flow V2 (37)
    out_flow_component(22,2) = 0; % Box V (22), flow V2 (flow 37 is going to Box 19, which is 2nd on list in out_index for Box 22 ) OK
    in_flow_component(19,2) = V2_new; % Box S (19), flow V2 (flow 37) is coming from Box 22, which is 2nd on list in in_index for Box 19  OK
    
    % Flow V3 (38)
    out_flow_component(22,3) = 0; % Box V (22), flow V3 (flow 38 is going to Box 23, which is 3rd on list in out_index for Box 22 ) OK
    in_flow_component(23,3) = V3_new; % Box W (23), flow V3 (flow 38) is coming from Box 22, which is 3rd on list in in_index for Box 23 OK
    
    % Flow T1 (34)
    out_flow_component(20,1 ) = 0; % Box T (20), flow T1 (flow 34) is going to Box 23, which is 1st on list in out_index for Box 20 OK
    in_flow_component(23,1) = T1_new; % Box W (23), flow T1 (flow 34) is coming from Box 20 (T), which is 1st on list in in_index for Box 23 OK
    
    % Flow Q2 (30)
    in_flow_component(20,1) = Q2_original - T1_new; % Box T (20), flow Q2 (flow 30) is coming from Box 17 (Q), which is 1st on list in in_index for Box 23
    
    %% Special case: flow 11 must be <50% flow of Box C and <= 16.5% CAGR from 2021
    
    %Calculate flow at time t based on initial value with 16.5% CAGR; set
    %this as max allowed flow at this t
    
    E1_original = interaction_flow_at_t(11,1);
    half_total_Box_C_output = (interaction_flow_at_t(5) + interaction_flow_at_t(6)) * 0.5;
    
    if E1_original > half_total_Box_C_output
        E1_new = half_total_Box_C_output;
    else
        E1_new = E1_original;
    end
    
    if baseline_run == 0
        
        if fix(t)<=4
            flow_limit_11 = baseline_individual_flow(11, plastic_type) * 1.02^( fix(t)-1 );
        elseif fix(t)>4
            flow_limit_11 = (baseline_individual_flow(11, plastic_type) * 1.02^5) * 1.165^( fix(t)-4 );
        end
        
        if E1_new >= flow_limit_11
            E1_new = flow_limit_11;
        end
    end
    
    % update flows
    interaction_flow_at_t(11,1) = E1_new;
    
    out_flow_component(5,2) = E1_new; % Box E (5),  flow E1 (flow 11) is going to Box K (11), which is 2nd on list in out_index for Box 5
    in_flow_component(11,2) = E1_new; % Box K (11), flow E1 (flow 11) is coming from Box E (5), which is 2nd on list in in_index for Box 11
    
    %% IF FLOWS EXCEED THE MAXIMUM SET CAPACITY OF THE BOX, SCALE DOWN RELEVANT FLOWS
    % FLOWS IN AND FLOWS OUT
    
    in_flow_component_temp = in_flow_component;
    out_flow_component_temp = out_flow_component;
    interaction_flow_temp = interaction_flow_at_t;
    
    for i = 3: 16 %upper limit could be n_stocks, but no need to apply this condition for unmanaged waste
        
        % default to flow not exceeded
        test_box_outflow = false();
        
        % check if total box(i) outflows exceed the maximum_capacity(i)
        if i==test_box
            fprintf('total_outflow(%d)=%f  outflow capacity=%f\n', i, sum(out_flow_component(i,:) ), box_flow_capacity(current_year_number, i, plastic_type));
            fprintf('t=%f\n', t);
            fprintf('%f\n', out_flow_component(i,1:4));
            fprintf('%f\n', ( sum(out_flow_component(i,:) ) - box_flow_capacity(current_year_number, i, plastic_type)  ));
        end
        
        if ( sum(out_flow_component(i,:) ) - box_flow_capacity(current_year_number, i, plastic_type)  ) > 0.001...
                && t >= box_flow_capacity_t_start(i)...
                && t <= box_flow_capacity_t_end(i)
            
            test_box_outflow = true();
        end
        
        
        
        % -------------------- STANDARD VERSION -----------------------------------
        % {
        if test_box_outflow == true()
            
            total_outflow_ratio = sum(out_flow_component(i,:) ) / box_flow_capacity(current_year_number, i, plastic_type);
            if i==test_box
                fprintf('total_outflow(%d)=%f  capacity=%f   total_outflow_ratio=%f\n',...
                    i, sum(out_flow_component(i,:) ), box_flow_capacity(current_year_number,...
                    i, plastic_type), total_outflow_ratio);
            end
            % scale all outflows from box by total_outflow_ratio
            out_flow_component_temp(i,:) = out_flow_component(i,:) ./ total_outflow_ratio;
            
            % find those same flows as inflows
            for j=1:n_outflows(i)
                pos = find(in_interaction_number == out_interaction_number(i,j) );
                
                % scale the flow when it's an inflow
                in_flow_component_temp(pos) = in_flow_component(pos) / total_outflow_ratio;
                if i==test_box
                    fprintf('Infows: %d  pos = %d,  in_flow_component(pos) = %f\n', j, pos, in_flow_component(pos));
                    out_interaction_number(i,j)
                end
                
                % scale the same interaction flow
                interaction_flow_temp(out_interaction_number(i,j)) = interaction_flow_at_t(out_interaction_number(i,j))...
                    / total_outflow_ratio;
            end
        end
        %}
        % -------------------- END: STANDARDL VERSION -----------------------------
        
        
        % -------------------- ALTERNATIVE VERSION --------------------------------
        % This version constrains the period immediately after lifting constraints,
        % to prevent unrealistically large jumps in flow rates
        %{
        post_constraint_period = false();
        if t > box_flow_capacity_t_end(i) && t < box_flow_capacity_t_end(i)+100 %(for rest of simulation)
            post_constraint_period = true();
        end
        
        if test_box_outflow == true() || post_constraint_period == true()
            
            if test_box_outflow == true()
                
                total_outflow_ratio = sum(out_flow_component(i,:) ) / box_flow_capacity(current_year_number, i, plastic_type);
            
            elseif post_constraint_period == true()

                %required_limit = box_flow_capacity(box_flow_capacity_t_end(i), i, plastic_type)...
                %    + ( (t - box_flow_capacity_t_end(i)) * 0.2 * box_flow_capacity(box_flow_capacity_t_end(i), i, plastic_type) );
                required_limit = box_flow_capacity(box_flow_capacity_t_end(i), i, plastic_type) * (1 + 0.15)^((t - box_flow_capacity_t_end(i)));
                
                if sum(out_flow_component(i,:)) > 0 && box_flow_capacity(box_flow_capacity_t_end(i), i, plastic_type)<inf
                    total_outflow_ratio = max( 1, sum(out_flow_component(i,:)) / required_limit);
                else
                    total_outflow_ratio = 1;
                end
                if false %total_outflow_ratio > 1 || total_outflow_ratio <=0
                    fprintf('%d : t=%f\tratio=%f\tsum=%f\treq=%f\tend_capacity=%f\n', ...
                        i, t, total_outflow_ratio, sum(out_flow_component(i,:)),...
                        required_limit, box_flow_capacity(box_flow_capacity_t_end(i), i, plastic_type));
                end
            end
            
            if i==test_box
                fprintf('total_outflow(%d)=%f  capacity=%f   total_outflow_ratio=%f\n',...
                    i, sum(out_flow_component(i,:) ), box_flow_capacity(current_year_number,...
                    i, plastic_type), total_outflow_ratio);
            end
            % scale all outflows from box by total_outflow_ratio
            out_flow_component_temp(i,:) = out_flow_component(i,:) ./ total_outflow_ratio;
            
            % find those same flows as inflows
            for j=1:n_outflows(i)
                pos = find(in_interaction_number == out_interaction_number(i,j) );
                
                % scale the flow when it's an inflow
                in_flow_component_temp(pos) = in_flow_component(pos) / total_outflow_ratio;
                if i==test_box
                    fprintf('Infows: %d  pos = %d,  in_flow_component(pos) = %f\n', j, pos, in_flow_component(pos));
                    out_interaction_number(i,j)
                end
                % scale the same interaction flow
                interaction_flow_temp(out_interaction_number(i,j)) = interaction_flow_at_t(out_interaction_number(i,j))...
                    / total_outflow_ratio;
            end
        end
        %}
        % -------------------- END: Version that constrains period immediately after lifting constraints ------------------------------
        
        in_flow_component = in_flow_component_temp;
        out_flow_component = out_flow_component_temp;
        interaction_flow_at_t = interaction_flow_temp;
        
    end
    
    %% Turn off boxes if mass at limit
    
    year_fraction = t - fix(t);
    
    % Check for new year and reset switches
    %if fix(t) ~= year_counter
    if fix(t) > year_counter
        year_counter = fix(t);
        switched_on = zeros(n_stocks, 1);
        switched_off = zeros(n_stocks, 1);
        %fprintf('\nSwitches reset at: %f\n', t );
    end
    
    
    % Close boxes to further inputs if mass reaches threshold
    for i=3:16
        %  3: Boxes 1 and 2 are defined externally, so these conditions don't apply
        % 16: This condition doesn't apply to unmanaged waste, so can ignore box >16
        
        if year_fraction < 0.10 && switched_on(i) == 0 %turn on only at start of year, if not already on
            
            if call_origin == 1 % integration phase
                
                if ismember(i, finite_sinks)
                    % FINITE SINK, therefore FLOW CAPACITY
                    test_value = (box_total_capacity(current_year_number, i, plastic_type) );
                    
                else
                    test_value = inf; % keep on at the start of the year
                    
                end
                
                switched_on(i) = 1;
                box_active(i) = 1; %turn the box back on as there is spare capacity
                switched_on_times(current_year_number, i, plastic_type) = t;
                if i==test_box
                    fprintf('\nBox %d turned on at t=%f\tM(%d)=%f\ttest_value=%f\n', i, t, i, M(i), test_value);
                end
                
            elseif call_origin ==2 % flows phase
                switched_on(i) = 1;
                box_active(i) = 1; %turn the box back on as there is spare capacity
            else
            end
        end
        
        if (year_fraction >= 0.10) && (year_fraction < 0.90) && (switched_off(i) == 0)  %only turn off later in the year if not already off
            
            % test value compares mass present to flow rate - switch off if
            % mass present is more than can be processed in remainder of the
            % year; confined away from edges of the year to avoid complications
            % due to integrator moving fowards/backwards in time (to improve
            % accuracy, especially around the discontinuity with capacity
            % increases), which complicates registering switch on/off
            
            if call_origin == 1 % integration phase
                
                if ismember(i, finite_sinks)
                    test_value = (box_total_capacity(current_year_number, i, plastic_type) );
                    %set time start/end for finite sinks
                    t_start = box_capacity_t_start(i);
                    t_end = box_capacity_t_end(i);
                else
                    %set time start/end for non-finite sinks
                    t_start = box_flow_capacity_t_start(i);
                    t_end = box_flow_capacity_t_end(i);
                    if box_flow_capacity(current_year_number, i, plastic_type) == inf
                        test_value = inf;
                    else
                        test_value = ( ( box_flow_capacity(fix(t)+1,i, plastic_type) - (box_flow_capacity(fix(t)+1,i, plastic_type)*0.05) )...
                            * (1 - (t - fix(t))) )...
                            + (box_flow_capacity(fix(0)+1,i, plastic_type)*0.05);
                    end
                end
                
                if M(i ) >= test_value && t >= t_start && t <= t_end
                    switched_off(i) = 1;
                    box_active(i) = 0; %turn the box off as it's full
                    if i==test_box2
                        fprintf('Box %d turned off at t=%f\tM(%d)=%f\ttest_value=%f\t plastic type %d\n', i, t, i, M(i), test_value, plastic_type);
                    end
                    
                    if switched_off_times(current_year_number, i, plastic_type) == inf %only update if not already changed from default
                        switched_off_times(current_year_number, i, plastic_type) = t;
                    end
                end
                
            elseif call_origin == 2 % flows phase
                if t >= switched_off_times(current_year_number, i, plastic_type)
                    switched_off(i) = 1;
                    box_active(i) = 0;
                    if i==test_box2
                        fprintf('\nFlows:  Box %d turned off at t=%f\t plastic type %d\n', i, t, plastic_type);
                    end
                end
                
            else
            end
            
        end
    end
    
    % Apply constraints for total flow through box (from box)
    % Note: checking for total flows has to be done at the end because all flows need
    % to be known, so can't do it flow-by-flow above
    
    
    % Scale the recycled component of interaction_flow
    for j = 1 : length(recycling_inflows)
        interaction_flow_at_t(recycling_inflows(j),1) = interaction_flow_at_t(recycling_inflows(j),1) * dt;
    end
    
    
end % (only for macro)


%% Final flow calculations to return
for j = 1:n_stocks
    in_flow(j) = sum(in_flow_component(j,:) );
    out_flow(j) = sum(out_flow_component(j,:) );
end

interaction_flow_at_t = interaction_flow_at_t';


end

