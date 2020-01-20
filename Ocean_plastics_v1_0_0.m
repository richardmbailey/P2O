%Read in config info
if baseline_run == 1
    fprintf('CREATING BASELINE DATA\n\n');
else
    fprintf('RUNNING SIMULATION\n\n');
end

fprintf('\tScenario #%d, Zone #%d\n\n', s, z);
fprintf('\tReading configuration files...');

if micro
    config_prefix = '/config_files/';
end

ff = fullfile(pwd, [config_prefix, 'basic_info.csv']);
basic_info = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'box_conditions.csv']);
box_conditions = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'stock_stock_interactions.csv']);
stock_stock_interactions = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'flow_flow_interactions.csv']);
flow_flow_interactions = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'interaction_parameters.csv']);
interaction_parameters = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'time_series_imported.csv']);
arrow_timeseries_original = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'production_calculated_timeseries.csv']);
production_calculated_timeseries = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'box_flow_capacity.csv']);
box_flow_capacity_table = csvread(ff);

ff = fullfile(pwd, [config_prefix, 'stock_names.txt']);
fileID = fopen(ff,'r');formatSpec = '%s';
S_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
stock_names = vertcat(S_names{:});


if macro % if micro, no econ/GHG/jobs
    
    % Demand
    demand_prefix =[config_prefix2, num2str(1), '_Zone_', num2str(z), '_'];
    % (e.g. HI_Urban_BAUII_Plastic_2_Zone_1_Costs.csv)
    
    ff = fullfile(pwd, [demand_prefix, 'Demand.csv']);
    demand_table_original = csvread(ff);
    
    % Proportions
    ff = fullfile(pwd, [demand_prefix, 'Proportions.csv']);
    zone_proportions = csvread(ff);
    
    % Costs
    ff = fullfile(pwd, [demand_prefix, 'Costs.csv']);
    costs_table_original = csvread(ff);
    
    % GHG & Jobs
    ff = fullfile(pwd, [demand_prefix, 'GHGJobs.csv']);
    GHGJobs_table_original = csvread(ff);
    
    % GHG & Jobs
    ff = fullfile(pwd, [demand_prefix, 'Prices.csv']);
    prices_timeseries_original = csvread(ff);
    
    % CAPEX timeseries
    ff = fullfile(pwd, [demand_prefix, 'CAPEX.csv']);
    CAPEX_timeseries_original = csvread(ff);
    
    % GHG rates timeseries
    ff = fullfile(pwd, [demand_prefix, 'GHGTimeseries.csv']);
    GHG_timeseries_original = csvread(ff);
    
    % Jobs rates timeseries
    ff = fullfile(pwd, [demand_prefix, 'JobsTimeseries.csv']);
    Jobs_timeseries_original = csvread(ff);

end

%basic info
duration = basic_info(1);
n_plastic_types = sum(basic_info(2,:)>0); %numel(basic_info(2,:));
plastic_types = basic_info(2,1:sum(basic_info(2,:)>0));

check_exists = exist('number_MC_iterations', 'var');
if check_exists == 0 && baseline_run == 0
    n_MC_iterations = basic_info(3);
elseif check_exists == 1 && baseline_run == 0
    n_MC_iterations = number_MC_iterations;
else
end

MC_distribution_type = basic_info(4); %1=Normal, 2=Uniform
production_calculation_type = basic_info(5); %1=Parameters, 2=timeseries from Master sheet
create_figures = basic_info(6); %true() or false()
start_date = basic_info(7);
output_resolution = basic_info(8);
n_stocks = basic_info(9);
production_interaction = basic_info(10);
waste_generated_box_number = basic_info(11);
imports_interaction = basic_info(12);
finite_sinks = basic_info(13,1:sum(basic_info(13,:)>0));
pedigree_values = basic_info(14:18);

    
rng('shuffle'); % random seed, for 'on-the-fly' random number generation

if do_one_at_a_time == true()
    % 3 variants used to assess linear sensitivity (gradient)
    n_MC_iterations = 3;
end

if baseline_run == 1 % over-write these variables
    duration = 1;
    n_MC_iterations = 1;
    create_figures = 0;
end

n_datapoints = (duration * datapoints_per_year) + 0; %number of datapoints per year in output
tspan = linspace(0,duration,n_datapoints); % time reference points for ODE solution
dt = duration / n_datapoints;

n_flows = sum(sum(stock_stock_interactions>0));
n_interactions = sum(sum(stock_stock_interactions>0));
rnd_dist_type = MC_distribution_type;

Econ_sensitivity_table = zeros(duration, Econ_table_n_cols, 3);
flow_sensitivity_table = zeros(duration, n_flows, 3);

if baseline_run == 1
    initial_M = zeros(n_stocks, n_plastic_types); %assume all boxes start at zero at t=0, for all plastic types
elseif baseline_run == 0
    %initial_M = zeros(n_stocks, n_plastic_types);
    initial_M = initial_mass_from_baseline_run * 1;
end

global box_active; % to keep an account during integration
global switched_on; % boxes off/on during flow calculations
global switched_off; % boxes off/on during flow calculations
global switched_on_times;
global switched_off_times;
global year_counter;
global rand_num_number;

%global baseline_individual_flow;

%baseline_individual_flow = ones(n_flows, n_plastic_types) .* inf;
box_active = ones(n_stocks,1);
switched_on = zeros(n_stocks, 1);
switched_off = zeros(n_stocks, 1);
switched_on_times = zeros(duration, n_stocks, n_plastic_types);
switched_off_times = ones(duration, n_stocks, n_plastic_types) .* inf;
year_counter = 0;


%% Allocate interaction / flow parameters

intervention_t_start = zeros(n_interactions, n_plastic_types);
intervention_t_end = zeros(n_interactions, n_plastic_types);

plug = zeros(n_flows,n_plastic_types);

time_series_pedigree = zeros(n_interactions, n_plastic_types);

a_mean = zeros(n_interactions, n_plastic_types);
a_error = zeros(n_interactions, n_plastic_types);
a = zeros(n_interactions, n_plastic_types);

b_mean = zeros(n_interactions, n_plastic_types);
b_error = zeros(n_interactions, n_plastic_types);
b = zeros(n_interactions, n_plastic_types);

c_mean = zeros(n_interactions, n_plastic_types);
c_error = zeros(n_interactions, n_plastic_types);
c = zeros(n_interactions, n_plastic_types);

d_mean = zeros(n_interactions, n_plastic_types);
d_error = zeros(n_interactions, n_plastic_types);
d = zeros(n_interactions, n_plastic_types);

function_type = zeros(n_interactions, n_plastic_types);

max_annual_flow_rate = zeros(n_interactions, n_plastic_types);
equation_or_timeseries = zeros(n_interactions, n_plastic_types);

enforced_proportion = zeros(n_interactions, n_plastic_types);
processing_rate = zeros(n_interactions, n_plastic_types);
relative_absolute = zeros(n_interactions, n_plastic_types);

production_rates = zeros(n_datapoints, n_plastic_types);
MC_production_rates = zeros(n_datapoints, n_plastic_types, n_MC_iterations);

Mass = zeros(n_datapoints, n_stocks, n_plastic_types, n_MC_iterations);

imports_rates = zeros(n_datapoints, n_plastic_types);

MC_arrow_timeseries = zeros( (n_plastic_types * n_interaction_arrows), duration, n_MC_iterations);


%% -------- Box mass capacities --------------------------------------------

% Assumption: the multiplicaton factors are the same for all plastic types

if baseline_run == 1
    
    box_flow_capacity =  ones(duration + 1, n_stocks, n_plastic_types) * inf;
    box_total_capacity = ones(duration + 1, n_stocks, n_plastic_types) * inf;
    
    box_capacity_t_start = zeros(n_stocks,1);
    box_capacity_t_end = ones(n_stocks,1) * duration;
    
    box_flow_capacity_t_start = zeros(n_stocks,1);
    box_flow_capacity_t_end = ones(n_stocks,1) * duration;
    
elseif baseline_run == 0
    
    % ------------- Box mass capacities:
    box_total_capacity_multiplier = box_conditions(1:n_stocks,1);
    %enforce 'no limits' if none set (blanks read as zero)
    box_total_capacity_multiplier(box_total_capacity_multiplier == 0) = inf;
    box_capacity_CAGR = box_conditions(1:n_stocks,2);
    box_capacity_t_start = box_conditions(1:n_stocks,3);
    box_capacity_t_end = box_conditions(1:n_stocks,4);
    
    box_total_capacity = zeros(duration, n_stocks, n_plastic_types);
    
    for k = 1:n_plastic_types
        for j = 1:n_stocks
            
            if intersect(j,finite_sinks) > 0 % if j is one of the finite sinks
                if box_capacity_CAGR(j) > 0 %if CAGR>0,
                    
                    % scale the initial mass
                    initial_mass = baseline_mass(1,j,k) /...
                        ( box_capacity_CAGR(j) * (1 + box_capacity_CAGR(j)) );
                    baseline_mass(1,j,k) = initial_mass;
                    
                end
                box_total_capacity(1,j,k) = end_first_yr_mass(1,j,k); %<--
            end
            
        end
    end
    
    for k = 1:n_plastic_types
        for i=1:duration
            for j = 1:n_stocks
                
                t_real = i;
                t_for_calc = t_real;
                
                if isempty( intersect(j,finite_sinks) ) == true() %if it's not a finite sink
                    
                    if baseline_mass(1,j,k) == 0
                        box_total_capacity(i+1,j,k) = inf;
                    elseif baseline_mass(1,j,k) > 0
                        box_total_capacity(i+1,j,k) = baseline_mass(1,j,k)...
                            * box_total_capacity_multiplier(j) * ((1 + box_capacity_CAGR(j))^(t_for_calc-1));
                    end
                    
                elseif isempty( intersect(j,finite_sinks) ) == false() %if it's a finite sink
                    
                    if box_capacity_CAGR(j) == 0 %if no growth in capacity, just add base amount each year
                        if end_first_yr_mass(1,j,k) == 0
                            box_total_capacity(i+1,j,k) = inf;
                        else
                            box_total_capacity(i+1,j,k) = ...
                                end_first_yr_mass(1,j,k) * box_total_capacity_multiplier(j) + ((i) * baseline_mass(1,j,k));
                        end
                        
                    elseif abs(box_capacity_CAGR(j)) > 0
                        
                        box_total_capacity(i+1,j,k) =...
                            box_total_capacity(i,j,k)...
                            + ( ( baseline_mass(1,j,k) * box_total_capacity_multiplier(j) * ((1 + box_capacity_CAGR(j))^(i+1)) )...
                            - ( baseline_mass(1,j,k) * box_total_capacity_multiplier(j) * ((1 + box_capacity_CAGR(j))^(i+0)) ) );
                    end
                else
                end
            end
        end
    end
    
    box_total_capacity(1,:,:) = inf; %no constraints in first year
    
    % ------------- Box_flow_capacites:
    box_flow_capacity_multiplier = box_flow_capacity_table(1:n_stocks,1);
    %enforce 'no limits' if none set (blanks read as zero --> inf)
    box_flow_capacity_multiplier(box_flow_capacity_multiplier == 0) = inf;
    box_flow_capacity_CAGR = box_flow_capacity_table(1:n_stocks,2);
    box_flow_capacity_t_start = box_flow_capacity_table(1:n_stocks,3);
    box_flow_capacity_t_end = box_flow_capacity_table(1:n_stocks,4);
    
    box_flow_capacity = zeros(duration + 1, n_stocks, n_plastic_types);
    box_flow_capacity(1,:,:) = inf; %no constraints in first year
    
    for k = 1:n_plastic_types
        
        for i=1:duration + 1
            
            for j = 1:n_stocks 
                
                t_real = i;
                t_for_calc = t_real;
                
                if baseline_flow(1,j,k) == 0
                    
                    box_flow_capacity(i+1,j,k) = inf;
                    
                elseif baseline_flow(1,j,k) > 0
                    
                    box_flow_capacity(i+0,j,k) = baseline_flow(1,j,k)...
                        * box_flow_capacity_multiplier(j) * ((1 + box_flow_capacity_CAGR(j))^(t_for_calc-1));
                end
            end
        end
    end
    
end

% Load the flow parameters
for i = plastic_types
    
    plug(:,i) = interaction_parameters(1:n_interactions,( 1 + ((i-1)*(n_interaction_parameters+1) )));
    
    relative_absolute(:,i) = interaction_parameters(1:n_interactions,( 2 + ((i-1)*(n_interaction_parameters+1) )));
    
    time_series_pedigree(:,i) = interaction_parameters(1:n_interactions,( 6 + ((i-1)*(n_interaction_parameters+1) )));
    
    max_annual_flow_rate(:,i) = interaction_parameters(1:n_interactions,( 13 + ((i-1)*(n_interaction_parameters+1) )));
    equation_or_timeseries(:,i) = interaction_parameters(1:n_interactions,( 14 + ((i-1)*(n_interaction_parameters+1) )));
    
    enforced_proportion(:,i) = interaction_parameters(1:n_interactions,( 15 + ((i-1)*(n_interaction_parameters+1) )));
    
    processing_rate(:,i) = interaction_parameters(1:n_interactions,( 16 + ((i-1)*(n_interaction_parameters+1) )));
    
    function_type(:,i) = interaction_parameters(1:n_interactions,( 17 + ((i-1)*(n_interaction_parameters+1) )));
end


% MC coefficients for opex, capex, jobs, GHG
% Tables for MC results
MC_final_flows = zeros(n_MC_iterations, n_flows, max(plastic_types));
MC_final_masses = zeros(n_MC_iterations, n_stocks, max(n_plastic_types));
MC_interaction_flows = zeros(n_datapoints, n_interactions,n_plastic_types,n_MC_iterations);

fprintf('complete\n');


% ------------------------ Monte Carlo iterations -----------------------

if macro
    % ---- Make an MC'd version of the demand table -----
    demand_pedigree = demand_table_original(:,1);
    
    % create matrix of coefficients for variation in demand parameter timeseries
    [n_rows_demand_table, ~] = size(demand_table_original);
    demand_table_MC_multipliers = ones(n_rows_demand_table, n_MC_iterations);
    
    % create multiplier for each row of the demand,costs/GHG table, for each required MC run
    % (only do this if a MC, not if it's a one-at-a-time sensitivity run)
    if do_one_at_a_time == false()
        for i = 2:n_MC_iterations-1   % first/last MC iteration run under central values
            if rand_scheme == 1 %individual rnd numbers
            demand_table_MC_multipliers(:,i) = MC_sample(rnd_dist_type, ones(n_rows_demand_table,1),...
                demand_pedigree, pedigree_values);
            elseif rand_scheme ==2 %long list
                demand_table_MC_multipliers(:,i) = MC_sample_list(rnd_dist_type, ones(n_rows_demand_table,1),...
                demand_pedigree, pedigree_values, rand_list);
            end
        end
    end  
end

% Create table for MC arrow timeseries
for dummy_mc = 1: n_MC_iterations
    MC_arrow_timeseries(:,:,dummy_mc) = arrow_timeseries_original(1: n_plastic_types * n_interaction_arrows, 1:duration);
end

%if n_MC_iterations > 1
    production_record_MC = zeros(duration, 33, n_MC_iterations);
%end

for mc = 1:n_MC_iterations

    tic();

    if n_MC_iterations > 1
        fprintf('\n\tMC run %d / %d :\n', mc, n_MC_iterations);
    end
    
    fprintf('\tDefining parameters...');
    
    % MC version of demand parameters
    
    % Macro: use full production calculation
    if macro
        calc_production;
    end

    
    if baseline_run == 0 && macro
    production_record_MC(:,1,mc) = waste_per_capita(1:duration)'; % 1 row (therefore transposed for production_record_MC)
    production_record_MC(:,2:4,mc) = waste_proportion(:,1:duration)'; %3 rows (plastic types)
    production_record_MC(:,5:7,mc) = gros_plastic_demand(:,:,z)'; %3 rows (plastic types)
    
    production_record_MC(:,8:10,mc) = reduce_eliminate'; %3 rows (plastic types)
    production_record_MC(:,13:15,mc) = reduce_reuse'; %3 rows (plastic types)
    production_record_MC(:,16:18,mc) = reduce_new_delivery'; %3 rows (plastic types)
    production_record_MC(:,19:21,mc) = substitute_paper'; %3 rows (plastic types)
    production_record_MC(:,22:24,mc) = substitute_coated_paper'; %3 rows (plastic types)
    production_record_MC(:,25:27,mc) = substitute_compostables'; %3 rows (plastic types)
    
    % Shift from Multi to Rigid
    production_record_MC(:,28,mc) = shift_multi_to_rigid';
    % Shift from Multi to Flexible
    production_record_MC(:,29,mc) = shift_multi_to_flexible';
    % Shift from Flexible to Rigid
    production_record_MC(:,30,mc) = shift_flexible_to_rigid';
    
    % net plastic demand (after all substitution/reduction/shift)
    production_record_MC(:,31,mc) = plastic_demand(1,:)';
    production_record_MC(:,32,mc) = plastic_demand(2,:)';
    production_record_MC(:,33,mc) = plastic_demand(3,:)';
    
    end    
       
    % MC version of the arrow timeseries
    
    if mc == 1
        % use central values
        arrow_timeseries = arrow_timeseries_original;
        
        if macro
            % save a copy of the R&S values under mean forcing values
            R_and_S_table = zeros(n_plastic_types * 6, duration);
            
            R_and_S_table(1:3,   1:duration) = reduce_eliminate(1:3,   1:duration);
            R_and_S_table(4:6,   1:duration) = reduce_reuse(1:3,   1:duration);
            R_and_S_table(7:9,   1:duration) = reduce_new_delivery(1:3,   1:duration);
            R_and_S_table(10:12, 1:duration) = substitute_paper(1:3,   1:duration);
            R_and_S_table(13:15, 1:duration) = substitute_coated_paper(1:3,   1:duration);
            R_and_S_table(16:18, 1:duration) = substitute_compostables(1:3,   1:duration);
            
            % Save econ table for this zone (e.g. A1_S1_Z1_ECON.csv)
            filename = [results_filename_prefix, '_RandS', '.csv'];
            file_path = fullfile(pwd, 'Output_files', filename);
            csvwrite(file_path, R_and_S_table );
        end
        
        if micro
            production_rates = arrow_timeseries(production_interaction+1,:)';
        end
        
    elseif mc > 1
        % Create MC variants of each timeseries (parallel scaled to proportional error)
        % Option to run one-at-a-time MC : do_one_at_a_time = true();
        
        for k = plastic_types
            
            if do_one_at_a_time == true()
                
                arrow_MC_list = x_analysis_flow; %one_at_a_time_flow_number;
                %MC_arrow_timeseries(:,:,mc) = arrow_timeseries_original(1: n_plastic_types * n_interaction_arrows, 1:duration);
                
            elseif do_one_at_a_time == false()
                
                arrow_MC_list = 1 : n_interactions;
                
            else
            end
                       
            jLB = 2 + (k - 1) * n_interaction_arrows;
                        
            for j = arrow_MC_list
                temp = arrow_timeseries_variant_list(...
                    arrow_timeseries_original(j + jLB - 1,1:duration),...
                    time_series_pedigree(j, k),...
                    pedigree_values,...
                    MC_distribution_type,...
                    do_one_at_a_time,...
                    mc, rand_list);
                
                MC_arrow_timeseries(j + jLB - 1,:,mc) = temp;
                
                % Special case for flows 1-4, as these are initial prescribed flows
                if j == 1 || j == 3
                    %calculate appropriately adjusted version of flow j+1 for its variant
                    flow_diff = arrow_timeseries_original(j + jLB - 1,1:duration) - temp;                  
                    MC_arrow_timeseries((j+1) + jLB - 1,:,mc) = arrow_timeseries_original((j+1) + jLB - 1,1:duration) + flow_diff;
                end
                if j == 2 || j == 4
                    %calculate appropriately adjusted version of flow j+1 for its variant
                    flow_diff = arrow_timeseries_original(j + jLB - 1,1:duration) - temp;                  
                    MC_arrow_timeseries((j-1) + jLB - 1,:,mc) = arrow_timeseries_original((j-1) + jLB - 1,1:duration) + flow_diff;
                end
                             
            end
        end
        
        arrow_timeseries = MC_arrow_timeseries(:,:,mc);
        
        if micro
            production_rates = arrow_timeseries(production_interaction+1,:)';
        end
        
    else
    end
     
    % in/out flows for each stock
    %(assuming system structure same for each Macro plastic type, or Micro source type)
    n_inflows = zeros(n_stocks,1);
    in_index = zeros(n_stocks, n_stocks);
    in_interaction_number = zeros(n_stocks, n_stocks);
    
    n_outflows = zeros(n_stocks,1);
    out_index = zeros(n_stocks, n_stocks);
    out_interaction_number = zeros(n_stocks, n_stocks);
    
    
    for i = 1: n_stocks
        
        %inputs: vertical components (columns) of adjacency matrix
        n_inflows(i) = sum(stock_stock_interactions(:,i)>0);
        in_index(i,1:n_inflows(i)) = find(stock_stock_interactions(:,i)>0);
        in_interaction_number(i,1:n_inflows(i)) = stock_stock_interactions(in_index(i,1:n_inflows(i)),i);
        
        %outputs: horizontal components (rows) of adjacency matrix
        n_outflows(i) = sum(stock_stock_interactions(i,:)>0);
        out_index(i,1:n_outflows(i)) = find(stock_stock_interactions(i,:)>0);
        out_interaction_number(i,1:n_outflows(i)) = stock_stock_interactions(i,out_index(i,1:n_outflows(i)));
    end
    
    in_source = in_index;
    out_source = out_index;
    
    for p=1:n_plastic_types  % ==> set plugs to eqn (to make plug values consistent)
        
        plugs = plug(:,p) == 1;
        
        %set all plugs to eqn
        equation_or_timeseries(plugs) = 1;
        
        % ==> set plug interactions as same rel/abs as non-plugs (per box), to make plug values consistent
        for i = 1 : n_stocks
            
            flows_out = out_interaction_number(i, 1:n_outflows(i));
            
            %check if any of these are plug numbers
            n_plugs = sum(plug(flows_out, p));
            
            plug_flows = [];
            
            %if plugs present...
            if n_plugs > 0
                
                %identify the plug flows
                p_indices = plug(flows_out, p) == 1;
                plug_flows = flows_out(p_indices);
                
                %find and count the non-plug flows
                non_p_indices = plug(flows_out, p) == 0;
                n_non_plugs = sum(non_p_indices);
                non_plug_flows = flows_out(non_p_indices);
                
                %check for user input error
                assert(n_non_plugs >0, ['All flows cannot be plugs, check Box ', num2str(i)] );
                
                % find non-plug rel/abs type
                rel_abs_type = relative_absolute(non_plug_flows(1), p);
                
                % set rel/abs type of plugs as the same as the first of the non-plugs
                %relative_absolute(plug_flows, p) = rel_abs_type;
                
            end
            
        end
        
    end
    
    fprintf('complete\n');
    
    % ---------------------------- Integration ----------------------------
    fprintf('\n\tNumerical integration:\n');
    options = odeset('AbsTol',1e-04, 'RelTol',1e-04, 'NonNegative', 1:n_stocks);
    scaling = 1;
    call_origin = 1;
    
    switched_on_times = zeros(duration, n_stocks, n_plastic_types);
    switched_off_times = ones(duration, n_stocks, n_plastic_types) .* inf;
    
    for plastic_type = plastic_types
        
        year_counter = 0;
        
        fprintf('\tPlastic type %d...', plastic_type);
        
        [t,M] = ode45(@(t,M) solve_ODEs_network_v0_2(t,tspan,n_datapoints,dt,M,scaling,plastic_type,...
            intervention_t_start, intervention_t_end, equation_or_timeseries, arrow_timeseries,...
            relative_absolute, a,b,c,d,box_total_capacity, box_flow_capacity, box_flow_capacity_t_start, box_flow_capacity_t_end,...
            box_capacity_t_start, box_capacity_t_end,...
            production_interaction, production_rates, imports_interaction, exports_interaction, imports_rates, waste_generated_box_number,...
            n_stocks,n_flows,max_annual_flow_rate,...
            n_inflows, n_outflows,...
            in_interaction_number,in_index,in_source,...
            out_interaction_number,out_index,out_source,...
            enforced_proportion, processing_rate,...
            stock_stock_interactions,flow_flow_interactions,...
            plug,function_type, mc, call_origin, zone_proportions, z, finite_sinks, micro, macro, baseline_run, baseline_individual_flow),...
            tspan,initial_M(:, plastic_type),options);
        
        Mass(:,:,plastic_type,mc) = M;
        MC_final_masses(mc,:,plastic_type) = Mass(end,:,plastic_type,mc);
        
        fprintf('complete\n');
    end
    
    
    % --------------- Calculate in and out flows --------------------------
    
    outflow = zeros(n_datapoints, n_stocks,n_plastic_types);
    inflow = zeros(n_datapoints, n_stocks,n_plastic_types);
    interaction_flow = zeros(n_datapoints, n_interactions,n_plastic_types);
    box_active = ones(n_stocks,1);
    
    
    scaling = 1; %to delete in later version
    call_origin = 2;
    
    fprintf('\n\tCalculating flows:\n');
    
    for j = plastic_types
        
        year_counter = 0;
        
        fprintf('\tPlastic type %d...', j);
        
        for i=1:n_datapoints
            
            [inflow(i,:,j), outflow(i,:,j), interaction_flow(i,:,j)] = calc_flows(Mass(i,:,j,mc),t(i),...
                tspan,n_datapoints,dt,scaling,...
                n_stocks,n_flows,j,stock_stock_interactions,...
                intervention_t_start, intervention_t_end, equation_or_timeseries,arrow_timeseries,...
                relative_absolute, a,b,c,d,box_total_capacity, box_flow_capacity, box_flow_capacity_t_start, box_flow_capacity_t_end,...
                box_capacity_t_start, box_capacity_t_end,...
                production_interaction, production_rates, imports_interaction, exports_interaction, imports_rates, waste_generated_box_number,...
                plug, function_type, enforced_proportion, processing_rate,max_annual_flow_rate,...
                flow_flow_interactions, in_interaction_number, out_interaction_number,...
                n_inflows, n_outflows, in_index, in_source, out_source, out_index, mc, call_origin,...
                zone_proportions, z, finite_sinks, micro, macro, baseline_run, baseline_individual_flow);
            
        end
        
        MC_final_flows(mc,:,j) = interaction_flow(end-1,:,j); % each final flows
        MC_interaction_flows(:,:,j,mc) = interaction_flow(:,:,j); % all flows
        
        fprintf('complete\n');
        
    end
    
    %over-write last values to avoid zeros
    interaction_flow(end,:,:) = interaction_flow(end-1,:,:);
    MC_interaction_flows(end,:,:,mc) = interaction_flow(end-1,:,:);
    
    % Define imports flows, for use in mass conservation calculation
    imports_rates(:, 1) = interaction_flow(:,imports_interaction);
    
    initial_mass_correction = sum(initial_M(:,1));
    
    %Check for conservation of mass for non-MC'd Macro runs
    %Note: micro mass conservation checked externally for present version
    if macro
        if mc == 1
            
            expected_Mass = ...
                + ( sum(interaction_flow(:,1,1)) * dt)...
                + ( sum(interaction_flow(:,2,1)) * dt)...
                + (sum(imports_rates(:,1,1)) * dt)...
                + initial_mass_correction;
            
            sum_Mass = sum(Mass(end,:,1,mc));
            check_Mass = sum_Mass / expected_Mass;
            
            
            if check_Mass > 1.01 || check_Mass < 0.99
                fprintf('\n\n\tMass conservation warning: Mass ratio = %f', check_Mass);
            else
                fprintf('\n\n\tMass conservation better than %.2g %%', 100*abs(1-check_Mass));
            end
        end
    end

    fprintf('\n\tCalculation time: %.1f seconds\n', toc());
    
    % ---------- If a sensitivity run, calculate Econ table for each MC run ----------
    if baseline_run == 0
        if do_sensitivity_analysis==true() && do_one_at_a_time==true()
            
            MC_stdev_masses = zeros(n_datapoints, n_stocks,n_plastic_types);
            MC_stdev_flows  = zeros(n_datapoints, n_interactions,n_plastic_types);
            
            %MC_interaction_flows_for_sensitivity = MC_interaction_flows(:,:,:,mc);
            
            % Compress data so it can be used by calc_econ2
            [output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
                compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
                Mass, MC_interaction_flows(:,:,:,mc), MC_stdev_flows, MC_stdev_masses);
            
            if macro
                % Run calc_econ2 to get the costs associated with each of the sensitivity runs
                calc_econ2;
                
                % Store econ table of results for fitting later
                Econ_sensitivity_table(:,:, mc) = Econ_table_MC(:,:,1);
            end
            
            % Store flows and mass summaries for fitting later
            flow_sensitivity_table(:,:, mc) = sum( output_interaction_flow(2:end,:,:), 3 );
            
        end
    end
    
end %...of MC iterations



%% Calculate stdev of MC runs for each plastic type

if baseline_run == 0
    
    fprintf('\n\nDATA SUMMARIES & OUTPUT FILES\n');
    
    fprintf('\n\tCalculating variance of MC iterations...');
    
    MC_stdev_masses = zeros(n_datapoints, n_stocks,n_plastic_types);
    MC_stdev_flows  = zeros(n_datapoints, n_interactions,n_plastic_types);
    
    for p = 1:n_plastic_types
        for j = 1:n_stocks
            for i = 1:n_datapoints
                MC_stdev_masses(i,j,p) = std( Mass(i,j,p,:) );
            end
        end
    end
    
    for p = 1:n_plastic_types
        for j = 1:n_interactions
            for i = 1:n_datapoints
                MC_stdev_flows(i,j,p) = std( MC_interaction_flows(i,j,p,:) );
            end
        end
    end
    
    fprintf(' complete.\n');
    
end


%% Compress to desired output resolution

if baseline_run==0 && do_sensitivity_analysis==false()
    % Compress data so it can be used by calc_econ2 (and set to 1st MC iteration for mean conditions)
    [output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
        compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
        Mass, MC_interaction_flows(:,:,:,1), MC_stdev_flows, MC_stdev_masses);
end


% store outputs for all MC runs

if n_MC_iterations > 1
    if baseline_run==0 && do_sensitivity_analysis==false()
        % Compress data so it can be used by calc_econ2 (and set to 1st MC iteration for mean conditions)
        MC_flows_plastic1 = zeros( (duration+1)*n_MC_iterations, n_interactions);
        if macro
            MC_flows_plastic2 = zeros( (duration+1)*n_MC_iterations, n_interactions);
            MC_flows_plastic3 = zeros( (duration+1)*n_MC_iterations, n_interactions);
        end
        
        for iMC = 1:n_MC_iterations
            [output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
                compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
                Mass, MC_interaction_flows(:,:,:,iMC), MC_stdev_flows, MC_stdev_masses);
            
            MC_flows_plastic1( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,1);
            MC_flows_plastic1( (iMC-1)*(duration+1)+26, : ) = NaN;
            
            if macro
                MC_flows_plastic2( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,2);
                MC_flows_plastic2( (iMC-1)*(duration+1)+26, : ) = NaN;
                
                MC_flows_plastic3( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,3);
                MC_flows_plastic3( (iMC-1)*(duration+1)+26, : ) = NaN;
            end
            
        end
    end
end


% ---------------------------- Save files --------------------------------

if baseline_run == 0
    
    fprintf('\n\n\tSaving timeseries data to files:\n\n')
    
    calculate_econ = true(); %could be added to conig files if necessary
    
    if calculate_econ == true() && macro
 
        calc_econ2; % produces Econ_table_MC            
        
        % Save econ table for this zone (e.g. A1_S1_Z1_ECON.csv)
        filename = [results_filename_prefix, '_ECON', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, Econ_table_MC(1:duration,1:length(Econ_table), 1) ); %central values
        
        % Save econ MC iterations
        econMC = zeros(n_econ_MC_iterations * Econ_table_rows, Econ_table_cols);
        for i = 1:n_econ_MC_iterations
            econMC( (i-1)*Econ_table_rows+1 : (i-1)*Econ_table_rows+25, :) = Econ_table_MC(:,:,i);
        end
        filename = [results_filename_prefix, '_ECON_MC_runs', '.csv'];
        file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
        csvwrite(file_path, econMC );
        
        
        fprintf('\tEcon/GHG/Jobs mean values: %s\n', filename);
        
        % Calculate stdev for the econ results
        Econ_table_stdev = zeros(Econ_table_rows, Econ_table_cols);
        for i = 1: Econ_table_rows
            for j = 1: Econ_table_cols
                Econ_table_stdev(i, j) = std(Econ_table_MC(i, j, :));
            end
        end
        
        % Save econ stdev table for this zone (e.g. A1_S1_Z1_ECON_SD.csv)
        filename = [results_filename_prefix, '_ECON_SD', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, Econ_table_stdev );
        
        fprintf('\tEcon/GHG/Jobs stdev values: %s\n', filename);
    
    elseif micro
                
        [output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
        compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
        Mass, MC_interaction_flows(:,:,:,:), MC_stdev_flows, MC_stdev_masses);
        
        % production of summed flows happens in cal_econ2 for macro
        % calc_econ2 not called for micro, so do summed_flows here
        summed_flow_types          = sum(output_interaction_flow(2:end,:,:), 3);
        summed_flow_types_stdev    = sum(output_interaction_flow_stdev(2:end,:,:), 3);
        summed_mass_types          = sum(output_mass(2:end,:,:), 3);
        summed_mass_types_stdev    = sum(output_mass_stdev(2:end,:,:), 3);
        
        filename = [results_filename_prefix, '_summed_flows', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, summed_flow_types );
        fprintf('\n\tSummed flows: %s\n', filename);

        filename = [results_filename_prefix, '_summed_flows_SD', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, summed_flow_types_stdev );
        fprintf('\tSummed flows stdev: %s\n', filename);

        filename = [results_filename_prefix, '_summed_mass', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, summed_mass_types );
        fprintf('\n\tSummed mass: %s\n', filename);

        filename = [results_filename_prefix, '_summed_mass_SD', '.csv'];
        file_path = fullfile(pwd, 'Output_files', filename);
        csvwrite(file_path, summed_mass_types_stdev );
        fprintf('\tSummed mass stdev: %s\n\n', filename);
        
    end
    
    
    for i = 1:n_plastic_types
        
        fprintf('\n');
%{        
        % Raw resolution flows for MC==1
        FlowsTable = array2table(MC_interaction_flows(:,:,i,1));
        filename = [results_filename_prefix, '_P', num2str(i), '_flow_RAW', '.txt'];
        file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
        writetable(FlowsTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
        
        %----------
        % Raw resolution mass for MC==1
        MassTable = array2table(Mass(:,:,i,1),...
            'VariableNames',stock_names);
        filename = [results_filename_prefix, '_P', num2str(i), '_mass_RAW', '.txt'];
        file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
        writetable(MassTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
%}        
        %----------
        % Summary flows
        SummaryFlowsTable = array2table(output_interaction_flow(2:end,:,i));
        filename = [results_filename_prefix, '_P', num2str(i), '_flow', '.txt'];
        file_path = fullfile(pwd, 'Output_files', filename);
        writetable(SummaryFlowsTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
        
        %----------
        % Summary flows stdev
        SummaryFlowsStdevTable = array2table(output_interaction_flow_stdev(2:end,:,i));
        filename = [results_filename_prefix, '_P', num2str(i), '_flow_SD', '.txt'];
        file_path = fullfile(pwd, 'Output_files', filename);
        writetable(SummaryFlowsStdevTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
        
        %----------
        % Summary mass
        SummaryMassTable = array2table(output_mass(2:end,:,i),...
            'VariableNames',stock_names);
        filename = [results_filename_prefix, '_P', num2str(i), '_mass', '.txt'];
        file_path = fullfile(pwd, 'Output_files', filename);
        writetable(SummaryMassTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
        
        %----------
        % Summary mass stdev
        SummaryMassStdevTable = array2table(output_mass_stdev(2:end,:,i),...
            'VariableNames',stock_names);
        filename = [results_filename_prefix, '_P', num2str(i), '_mass_SD', '.txt'];
        file_path = fullfile(pwd, 'Output_files', filename);
        writetable(SummaryMassStdevTable, file_path, 'Delimiter','\t');
        fprintf('\t%s\n', filename);
        
    end
%{    
    % Save box mass capacity data (e.g. A1_S1_Z1_Box_mass capacity.csv)
    filename = [results_filename_prefix, '_Box_mass_capacity', '.csv'];
    file_path = fullfile(pwd, 'Output_files', filename);
    csvwrite(file_path, box_total_capacity );
    
    % Save box flow capacity data (e.g. A1_S1_Z1_Box_flow_capacity.csv)
    filename = [results_filename_prefix, '_Box_flow_capacity', '.csv'];
    file_path = fullfile(pwd, 'Output_files', filename);
    csvwrite(file_path, box_flow_capacity );
 %}   
    
end

% Sensitivity calculations --> fill table

if baseline_run == 0
    if do_sensitivity_analysis == true()
        
        sensitivity_table_flow_x = zeros(n_plastic_types, 3, 1);
        sensitivity_table_flow_y = zeros(n_plastic_types, 3, n_interactions + 5);
        sensitivity_matrix_flow = zeros(n_plastic_types, n_interactions + 5); %[3,49]
        
        sensitivity_table_econ_x = zeros(n_plastic_types, 3, 1);
        sensitivity_table_econ_y = zeros(n_plastic_types, 3, Econ_table_n_cols + 3);
        sensitivity_matrix_econ = zeros(n_plastic_types, Econ_table_n_cols + 3);
        
               
        if do_one_at_a_time
            
            % ---------- Sensitivity of various flows ---------------------
            
            extra_y_interaction_flows = 25;
            
            y_analysis_flows = zeros(n_interactions + extra_y_interaction_flows, 5);
            y_analysis_flows(1:n_interactions,1) = 1:n_interactions;
            
            y_analysis_flows(n_interactions + 1, 1:3) = [34, 35, 38]; %ocean leakage: 45
            y_analysis_flows(n_interactions + 2, 1:2) = [29, 37]; %open burning: 46
            y_analysis_flows(n_interactions + 3, 1:2) = [19, 20]; %closed-loop MR: 47
            y_analysis_flows(n_interactions + 4, 1:2) = [21, 44]; %open-loop MR: 48
            y_analysis_flows(n_interactions + 5, 1:4) = [14, 15, 16, 17]; %formal sorting: 49
            
            
            for y_sensitivity_flow = 1: n_interactions + 5
                
               
                [sensitivity_x, sensitivity_y] = calc_sensitivity_one_at_a_time(...
                    analysis_year, x_analysis_flow, y_analysis_flows(y_sensitivity_flow, :),...
                    plastic_types, n_plastic_types, MC_interaction_flows, n_MC_iterations, dt,...
                    flow_sensitivity_table, flow_sensitivity_table, output_mass(2:end,:,:));


                sensitivity_table_flow_x(:, :, 1) = sensitivity_x;
                sensitivity_table_flow_y(:, :, y_sensitivity_flow) = sensitivity_y;
                
                % Calculate gradients on sensitivity data
                for k = 1 %plastic_types amalgamated in to total plastic

                    if (sensitivity_y(k,1) > sensitivity_y(k,2)) && ( sensitivity_y(k,1) > sensitivity_y(k,3) )
                        sensitivity_matrix_flow(k, y_sensitivity_flow) = 0;
                    elseif (sensitivity_y(k,1) < sensitivity_y(k,2)) && ( sensitivity_y(k,1) < sensitivity_y(k,3) )
                        sensitivity_matrix_flow(k, y_sensitivity_flow) = 0;
                    else
                        fit1 = polyfit(sensitivity_x(k,:), sensitivity_y(k,:), 1);
                        sensitivity_matrix_flow(k, y_sensitivity_flow) = fit1(1);
                    end
                    
                    
                end
            end
            
            % ------------- Sensitivity for econ parameters ---------------
            
            y_analysis_flows = zeros(Econ_table_n_cols, 25);
            y_analysis_flows(1:Econ_table_n_cols,1) = 1: Econ_table_n_cols;
            
            % 'Total' CAPEX (122)
            y_analysis_flows(Econ_table_n_cols + 1, 1:15) = [45,46,47,48,49,50,51,52,53,54,55,56,57,58,59];  
            % All except virgin plastic production and plastic conversion
            
            % 'Total' OPEX (123)
            y_analysis_flows(Econ_table_n_cols + 2, 1:19) = [24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]; 
            % All except virgin plastic production and plastic conversion
            
            % 'Total' REQ_INVESTMENT (124)
            y_analysis_flows(Econ_table_n_cols + 3, 1:15) = [73,74,75,76,77,78,79,80,81,82,83,84,85,86,87]; 
            % All except virgin plastic production and plastic conversion
            
            for y_sensitivity_econ = 1: Econ_table_n_cols + 3
                
                [sensitivity_x, sensitivity_y] = calc_sensitivity_one_at_a_time(...
                    analysis_year, x_analysis_flow, y_analysis_flows(y_sensitivity_econ, :),...
                    plastic_types, n_plastic_types, MC_interaction_flows, n_MC_iterations, dt,...
                    flow_sensitivity_table, Econ_sensitivity_table, output_mass(2:end,:,:));
                
                sensitivity_table_econ_x(:, :, 1) = sensitivity_x;
                sensitivity_table_econ_y(:, :, y_sensitivity_econ) = sensitivity_y;
                
                % Calculate gradients on sensitivity data
                for k = [1] % Note: total plastic
                    if (sensitivity_y(k,1) > sensitivity_y(k,2)) && ( sensitivity_y(k,1) > sensitivity_y(k,3) )
                        sensitivity_matrix_econ(k, y_sensitivity_econ) = 0;
                    elseif (sensitivity_y(k,1) < sensitivity_y(k,2)) && ( sensitivity_y(k,1) < sensitivity_y(k,3) )
                        sensitivity_matrix_econ(k, y_sensitivity_econ) = 0;
                    else
                        fit1 = polyfit(sensitivity_x(k,:), sensitivity_y(k,:), 1);
                        sensitivity_matrix_econ(k, y_sensitivity_econ) = fit1(1);
                    end
                end
            end
            
            [~,rank_econ]=ismember(sensitivity_matrix_econ(1,:),sort(sensitivity_matrix_econ(1,:),'descend'));
            
            % Plot sensitivity results
            if create_figures == 1
                figure(); bar(sensitivity_matrix_flow(1,:)); title('Flows sensitivity');
                figure(); bar(sensitivity_matrix_econ(1,:)); title('Econ sensitivity');
            end
            
        end
        
    end
end


if baseline_run == 0
     fprintf('\n\tMODEL CALCULATIONS COMPLETE\n\n');
end


%% save outputs for MC runs: flows, masses & production 

if n_MC_iterations > 1
    if baseline_run==0 && do_sensitivity_analysis==false()
        
        if macro
        % Save Production MC iterations
        [pr_rows, pr_cols, ~] = size(production_record_MC);
        
        productionMC_table = zeros(n_MC_iterations * pr_rows, pr_cols);
        for i = 1:n_MC_iterations
            productionMC_table( (i-1)*pr_rows+1 : (i-1)*pr_rows+25, :) = production_record_MC(:,:,i);
        end
        filename = [results_filename_prefix, '_PROD_RandS_MC_runs', '.csv'];
        file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
        csvwrite(file_path, productionMC_table );
        end
        
        
        % Compress data so it can be used by calc_econ2; set to 1st MC iteration for mean conditions
        MC_flows_plastic1 = zeros( (duration+1)*n_MC_iterations, n_interactions);
        MC_mass_plastic1 = zeros( (duration+1)*n_MC_iterations, n_stocks);
        if macro
            MC_flows_plastic2 = zeros( (duration+1)*n_MC_iterations, n_interactions);
            MC_mass_plastic2 = zeros( (duration+1)*n_MC_iterations, n_stocks);
            
            MC_flows_plastic3 = zeros( (duration+1)*n_MC_iterations, n_interactions);
            MC_mass_plastic3 = zeros( (duration+1)*n_MC_iterations, n_stocks);
        end
        
        for iMC = 1:n_MC_iterations
            
            fprintf('\n\tCreating output at required resolution...');
            
            [output_interaction_flow, output_mass, ~, ~] = ...
                compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
                Mass(:,:,:,iMC), MC_interaction_flows(:,:,:,iMC), MC_stdev_flows, MC_stdev_masses);
            
            MC_flows_plastic1( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,1);
            MC_flows_plastic1( (iMC-1)*(duration+1)+26, : ) = NaN;
            
            MC_mass_plastic1( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_mass(2:end,:,1);
            MC_mass_plastic1( (iMC-1)*(duration+1)+26, : ) = NaN;
            
            if macro
                MC_flows_plastic2( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,2);
                MC_flows_plastic2( (iMC-1)*(duration+1)+26, : ) = NaN;
                MC_mass_plastic2( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_mass(2:end,:,2);
                MC_mass_plastic2( (iMC-1)*(duration+1)+26, : ) = NaN;
                
                MC_flows_plastic3( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_interaction_flow(2:end,:,3);
                MC_flows_plastic3( (iMC-1)*(duration+1)+26, : ) = NaN;
                MC_mass_plastic3( (iMC-1)*(duration+1)+1 : (iMC-1)*(duration+1)+25, : ) = output_mass(2:end,:,3);
                MC_mass_plastic3( (iMC-1)*(duration+1)+26, : ) = NaN;
            end
            
        end
        
        % Save to file
        % (e.g. A1_S1_Z1_Plastic1_MC_flows.csv)
        
        if macro
            
            filename = [results_filename_prefix, '_Plastic1_MC_flows', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_flows_plastic1 );
            
            filename = [results_filename_prefix, '_Plastic1_MC_mass', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_mass_plastic1 );
                   
            filename = [results_filename_prefix, '_Plastic2_MC_flows', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_flows_plastic2 );
            
            filename = [results_filename_prefix, '_Plastic2_MC_mass', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_mass_plastic2 );
            
            filename = [results_filename_prefix, '_Plastic3_MC_flows', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_flows_plastic3 );
            
            filename = [results_filename_prefix, '_Plastic3_MC_mass', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_mass_plastic3 );
            
        end
        
        if micro
            
            filename = [results_filename_prefix, '_', char(micro_type), '_MC_flows', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_flows_plastic1 );
            
            filename = [results_filename_prefix, '_', char(micro_type), '_MC_mass', '.csv'];
            file_path = fullfile(pwd, 'Output_files', 'Raw', filename);
            csvwrite(file_path, MC_mass_plastic1 );
            
        end
        
    end
end

if baseline_run == 0
     fprintf('\n\n\tFinished.\n\n');
end



