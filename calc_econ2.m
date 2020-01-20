
if do_sensitivity_analysis && do_one_at_a_time
    n_econ_MC_iterations = 1;
    range_MC = 0;
else
    n_econ_MC_iterations = 500;
    range_MC = 1;
end

Econ_table = zeros(duration, Econ_table_n_cols);
Econ_table_MC = zeros(duration, Econ_table_n_cols, n_econ_MC_iterations);
    

%% MC versions of the Costs, GHGJobs tables

% Set up empty matraces for the MC versions of the costs and GHGJobs table
[n_rows_costs_table, n_cols_costs_table] = size(costs_table_original);
costs_table_MC = ones(n_rows_costs_table, n_cols_costs_table, n_econ_MC_iterations);

[n_rows_GHGJobs_table, n_cols_GHGJobs_table] = size(GHGJobs_table_original);
GHGJobs_table_MC = ones(n_rows_GHGJobs_table, n_cols_GHGJobs_table, n_econ_MC_iterations);

% Price rates timeseries
[n_rows_prices_timeseries, n_cols_prices_timeseries] = size(prices_timeseries_original);
n_cols_prices_timeseries = n_cols_prices_timeseries - 1; %1st column is pedigree number
prices_timeseries_MC = ones(n_rows_prices_timeseries, n_cols_prices_timeseries, n_econ_MC_iterations);
prices_pedigree = prices_timeseries_original(:,1);

% CAPEX rates timeseries
[n_rows_CAPEX_timeseries, n_cols_CAPEX_timeseries] = size(CAPEX_timeseries_original);
n_cols_CAPEX_timeseries = n_cols_CAPEX_timeseries - 1; %1st column is pedigree number
CAPEX_timeseries_MC = ones(n_rows_CAPEX_timeseries, n_cols_CAPEX_timeseries, n_econ_MC_iterations);
CAPEX_timeseries_pedigree = CAPEX_timeseries_original(:,1);
CAPEX_timeseries_type = costs_table_original(1:18, 17);

% GHG rates timeseries
[n_rows_GHG_timeseries, n_cols_GHG_timeseries] = size(GHG_timeseries_original);
n_cols_GHG_timeseries = n_cols_GHG_timeseries - 1; %1st column is pedigree number
GHG_timeseries_MC = ones(n_rows_GHG_timeseries, n_cols_GHG_timeseries, n_econ_MC_iterations);
GHG_timeseries_pedigree = GHG_timeseries_original(:,1);

% Jobs rates timeseries
[n_rows_Jobs_timeseries, n_cols_Jobs_timeseries] = size(Jobs_timeseries_original);
n_cols_Jobs_timeseries = n_cols_Jobs_timeseries - 1; %1st column is pedigree number
Jobs_timeseries_MC = ones(n_rows_Jobs_timeseries, n_cols_Jobs_timeseries, n_econ_MC_iterations);
Jobs_timeseries_pedigree = Jobs_timeseries_original(:,1);


for MCi = 1:n_econ_MC_iterations
    costs_table_MC(:,:,MCi) = costs_table_original; 
    GHGJobs_table_MC(:,:,MCi) = GHGJobs_table_original;
    prices_timeseries_MC(:,:,MCi) = prices_timeseries_original(:, 2:end);
    CAPEX_timeseries_MC(:,:,MCi) = CAPEX_timeseries_original(:, 2:end);
    GHG_timeseries_MC(:,:,MCi) = GHG_timeseries_original(:, 2:end);
    Jobs_timeseries_MC(:,:,MCi) = Jobs_timeseries_original(:, 2:end);
end

% Make changes for the MC'd versions of the costs and GHGJobs tables
for MCi = 2:n_econ_MC_iterations-range_MC %first/last one is the standard (mean) values

    % Costs table
    costs_table_MC(:, 1, MCi) = MC_sample(1, costs_table_original(:,1), ...
        costs_table_original(:,3), pedigree_values);
    costs_table_MC(:, 2, MCi) = MC_sample(1, costs_table_original(:,2), ...
        costs_table_original(:,3), pedigree_values);
    costs_table_MC(:, 4, MCi) = MC_sample(1, costs_table_original(:,4), ...
        costs_table_original(:,5), pedigree_values);
    costs_table_MC(:, 6, MCi) = MC_sample(1, costs_table_original(:,6), ...
        costs_table_original(:,7), pedigree_values);
    costs_table_MC(:, 8, MCi) = MC_sample(1, costs_table_original(:,8), ...
        costs_table_original(:,9), pedigree_values);
    costs_table_MC(:, 10, MCi) = MC_sample(1, costs_table_original(:,10), ...
        costs_table_original(:,12), pedigree_values);
    costs_table_MC(:, 11, MCi) = MC_sample(1, costs_table_original(:,11), ...
        costs_table_original(:,12), pedigree_values);
    costs_table_MC(:, 13, MCi) = MC_sample(1, costs_table_original(:,13), ...
        costs_table_original(:,15), pedigree_values);
    costs_table_MC(:, 14, MCi) = MC_sample(1, costs_table_original(:,14), ...
        costs_table_original(:,15), pedigree_values); 
    costs_table_MC(:, 16, MCi) = MC_sample(1, costs_table_original(:,16), ...
        costs_table_original(:,12), pedigree_values);
    
    % GHGJobs_table
    GHGJobs_table_MC(:, 1, MCi) = MC_sample(1, GHGJobs_table_original(:,1), ...
        GHGJobs_table_original(:,3), pedigree_values);
    GHGJobs_table_MC(:, 2, MCi) = MC_sample(1, GHGJobs_table_original(:,2), ...
        GHGJobs_table_original(:,3), pedigree_values);
    
    GHGJobs_table_MC(:, 4, MCi) = MC_sample(1, GHGJobs_table_original(:,4), ...
        GHGJobs_table_original(:,6), pedigree_values);
    GHGJobs_table_MC(:, 5, MCi) = MC_sample(1, GHGJobs_table_original(:,5), ...
        GHGJobs_table_original(:,6), pedigree_values);
    
    % Prices_timeseries
    for pt = 1:n_rows_prices_timeseries
        prices_timeseries_MC(pt, :, MCi) = stdev_scaling(prices_timeseries_original(pt,2:end), prices_pedigree(pt), MCi);
    end
    
    % CAPEX_timeseries
    for pt = 1:6
        CAPEX_timeseries_MC(pt, :, MCi) = stdev_scaling(CAPEX_timeseries_original(pt,2:end), CAPEX_timeseries_pedigree(pt), MCi);
    end
    
    % GHG_timeseries
    for pt = 1:17
        GHG_timeseries_MC(pt, :, MCi) = stdev_scaling(GHG_timeseries_original(pt,2:end), GHG_timeseries_pedigree(pt), MCi);
    end
    
    % Jobs_timeseries
    for pt = 1:17
        Jobs_timeseries_MC(pt, :, MCi) = stdev_scaling(Jobs_timeseries_original(pt,2:end), Jobs_timeseries_pedigree(pt), MCi);
    end

     
    
end

% FLOWS & stdev, MASSES & stdev (summed over all plastic types)
 
 % use MC==1 if a MC run; use MC-specific value (one of 3) if a sensitivity run
 if do_sensitivity_analysis==false()
    [output_interaction_flow, output_mass, output_interaction_flow_stdev, output_mass_stdev] = ...
            compress_data(baseline_run, duration, dt, output_resolution, datapoints_per_year, n_stocks, n_flows, n_plastic_types,...
            Mass, MC_interaction_flows(:,:,:,1), MC_stdev_flows, MC_stdev_masses);
 end
 
        
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


Virgin_plastic_production  = summed_flow_types(:,1) + summed_flow_types(:,2)...
    - summed_flow_types(:,19) - summed_flow_types(:,22);  %- summed_flow_types(:,44)

Virgin_plastic_production_stdev = combine_errors(summed_flow_types_stdev(:,1), summed_flow_types_stdev(:,2),...
    summed_flow_types_stdev(:,19), summed_flow_types_stdev(:,22) ); %- summed_flow_types_stdev(:,44)

Plastic_conversion         = summed_flow_types(:,1) + summed_flow_types(:,2); % total waste generated (post R&S)
Plastic_conversion_stdev   = combine_errors(summed_flow_types_stdev(:,1), summed_flow_types_stdev(:,2) );

Formal_collection          = summed_flow_types(:,3); %flow B1 (flow in to Box C)
Formal_collection_stdev    = summed_flow_types_stdev(:,3); %flow B1 (flow in to Box C)

Informal_collection        = summed_flow_types(:,4) + summed_flow_types(:,36); %
Informal_collection_stdev  = combine_errors(summed_flow_types_stdev(:,4), summed_flow_types_stdev(:,36) ); %

Formal_sorting             = summed_flow_types(:,5) + summed_flow_types(:,13) + summed_flow_types(:,18); %
Formal_sorting_stdev       = combine_errors( summed_flow_types_stdev(:,5), summed_flow_types_stdev(:,13), summed_flow_types_stdev(:,18) );

Closed_loop_MR             = summed_flow_types(:,7) + summed_flow_types(:,14); %
Closed_loop_MR_stdev       = combine_errors( summed_flow_types_stdev(:,7), summed_flow_types_stdev(:,14) ); %

Open_loop_MR               = summed_flow_types(:,8) + summed_flow_types(:,15); %
Open_loop_MR_stdev         = combine_errors( summed_flow_types_stdev(:,8), summed_flow_types_stdev(:,15) ); %

Chemical_conversion_P2P    = summed_flow_types(:,22); % arrow K1
Chemical_conversion_P2P_stdev = summed_flow_types_stdev(:,22);

Chemical_conversion_P2F    = summed_flow_types(:,23); % arrow K2
Chemical_conversion_P2F_stdev = summed_flow_types_stdev(:,23); % arrow K2

Total_chem_conversion = summed_flow_types(:,9) + summed_flow_types(:,11);
Total_chem_conversion_stdev = combine_errors( summed_flow_types_stdev(:,9), summed_flow_types_stdev(:,11) );


Thermal_treatment          = summed_flow_types(:,27); % arrow M1
Thermal_treatment_stdev    = summed_flow_types_stdev(:,27); % arrow M1

Engineered_landfills       = summed_flow_types(:,28); % arrow M2
Engineered_landfills_stdev = summed_flow_types_stdev(:,28); % arrow M2

Import_sorting             = summed_flow_types(:,18); % arrow H1
Import_sorting_stdev       = summed_flow_types_stdev(:,18); % arrow H1

Chemical_conversion        = summed_flow_types(:,22) + summed_flow_types(:,23) + summed_flow_types(:,24); % Total going out of box 11
Chemical_conversion_stdev  = combine_errors( summed_flow_types_stdev(:,22), summed_flow_types_stdev(:,23), summed_flow_types_stdev(:,24) );


% These are time series that need to be scaled by a single pedigree value
Reduce_Eliminate           = sum(reduce_eliminate, 1); % calculated in calc_production
Reduce_Reuse               = sum(reduce_reuse, 1); % calculated in calc_production
Reduce_New_delivery_models = sum(reduce_new_delivery, 1); % calculated in calc_production

Substitute_Paper           = sum(substitute_paper, 1); % calculated in calc_production
Substitute_Coated_paper    = sum(substitute_coated_paper, 1); % calculated in calc_production
Substitute_Compostables    = sum(substitute_compostables, 1); % calculated in calc_production


Actual_closed_loop_MR      = summed_flow_types(:,19);  % arrow I1
Actual_closed_loop_MR_stdev = summed_flow_types_stdev(:,19);  % arrow I1

Actual_open_loop_MR        = summed_flow_types(:,44);  % arrow J0
Actual_open_loop_MR_stdev  = summed_flow_types_stdev(:,44);  % arrow J0

Open_burning               = summed_flow_types(:,29) + summed_flow_types(:,37);
Open_burning_stdev         = combine_errors( summed_flow_types_stdev(:,29), summed_flow_types_stdev(:,37) );

% Create table of summed_flow_types_MC, for all the MC variants
[n_summed_flow_rows, n_summed_flow_cols] = size(summed_flow_types);
summed_flow_types_MC = zeros(n_summed_flow_rows, n_summed_flow_cols);


%% MC iterations
for MCi = 1:n_econ_MC_iterations
    
% set the correct costs and GHGJobs (MC) version of tables
costs_table = costs_table_MC(:, :, MCi);
GHGJobs_table = GHGJobs_table_MC(:, :, MCi);
prices_timeseries = prices_timeseries_MC(:, :, MCi);
CAPEX_timeseries = CAPEX_timeseries_MC(:, :, MCi);
GHG_timeseries = GHG_timeseries_MC(:, :, MCi);
Jobs_timeseries = Jobs_timeseries_MC(:, :, MCi);

Econ_table(:,1) = stdev_scaling( Virgin_plastic_production, Virgin_plastic_production_stdev, MCi );      
Econ_table(:,2) = stdev_scaling( Plastic_conversion, Plastic_conversion_stdev, MCi );                    
Econ_table(:,3) = stdev_scaling( Formal_collection, Formal_collection_stdev, MCi );                      
Econ_table(:,4) = stdev_scaling( Informal_collection, Informal_collection_stdev, MCi );
Econ_table(:,5) = stdev_scaling( Formal_sorting, Formal_sorting_stdev, MCi );
Econ_table(:,6) = stdev_scaling( Closed_loop_MR, Closed_loop_MR_stdev, MCi );
Econ_table(:,7) = stdev_scaling( Open_loop_MR, Open_loop_MR_stdev, MCi );
Econ_table(:,8) = stdev_scaling( Chemical_conversion_P2P, Chemical_conversion_P2P_stdev, MCi );
Econ_table(:,9) = stdev_scaling( Chemical_conversion_P2F, Chemical_conversion_P2F_stdev, MCi );
Econ_table(:,10) = stdev_scaling( Thermal_treatment, Thermal_treatment_stdev, MCi );
Econ_table(:,11) = stdev_scaling( Engineered_landfills, Engineered_landfills_stdev, MCi );
Econ_table(:,12) = stdev_scaling( Import_sorting, Import_sorting_stdev, MCi );

Econ_table(:,13) = Reduce_Eliminate; % Note: * these are varied (scaled) in Ocean_plastics_vx_x_x
Econ_table(:,14) = Reduce_Reuse; % *
Econ_table(:,15) = Reduce_New_delivery_models; % *
Econ_table(:,16) = Substitute_Paper; % *
Econ_table(:,17) = Substitute_Coated_paper; % *
Econ_table(:,18) = Substitute_Compostables; % *

Econ_table(:,19) = stdev_scaling( Actual_closed_loop_MR, Actual_closed_loop_MR_stdev, MCi );
Econ_table(:,20) = stdev_scaling( Actual_open_loop_MR, Actual_open_loop_MR_stdev, MCi );
Econ_table(:,21) = stdev_scaling( Open_burning, Open_burning_stdev, MCi );

column_counter = 21;

% OPEX costs
OPEX_Virgin_plastic_production  = opex_cost_timeseries( Virgin_plastic_production, costs_table(1,1), costs_table(1,2) );
OPEX_Plastic_conversion         = opex_cost_timeseries( Plastic_conversion, costs_table(2, 1), costs_table(2, 2) );
OPEX_Formal_collection          = opex_cost_timeseries( Formal_collection , costs_table(3, 1), costs_table(3, 2) );
OPEX_Informal_collection        = opex_cost_timeseries( Informal_collection, costs_table(4,1), costs_table(4,2) );
OPEX_Formal_sorting             = opex_cost_timeseries( Formal_sorting, costs_table(5,1), costs_table(5,2) );
OPEX_Closed_loop_MR             = opex_cost_timeseries( Closed_loop_MR, costs_table(6,1), costs_table(6,2) );
OPEX_Open_loop_MR               = opex_cost_timeseries( Open_loop_MR, costs_table(7,1), costs_table(7,2) );

[OPEX_Chemical_conversion_P2P, ~]...
    = capex_cost_timeseries( Chemical_conversion_P2P,   -1,  -1,  -1,...
    costs_table(8,2),   costs_table(8,1),  4, -1, -1, Total_chem_conversion  );

[OPEX_Chemical_conversion_P2F, ~]...
    = capex_cost_timeseries( Chemical_conversion_P2F,   -1,  -1,  -1,...
    costs_table(9,2),   costs_table(9,1),  4, -1, -1, Total_chem_conversion  );


OPEX_Thermal_treatment          = opex_cost_timeseries( Thermal_treatment, costs_table(10,1), costs_table(10,2) );
OPEX_Engineered_landfills       = opex_cost_timeseries( Engineered_landfills, costs_table(11,1), costs_table(11,2) );
OPEX_Import_sorting             = opex_cost_timeseries( Import_sorting, costs_table(12,1), costs_table(12,2) );
OPEX_Reduce_Eliminate           = opex_cost_timeseries( Reduce_Eliminate, costs_table(13,1), costs_table(13,2) );
OPEX_Reduce_Reuse               = opex_cost_timeseries( Reduce_Reuse, costs_table(14,1), costs_table(14,2) );
OPEX_Reduce_New_delivery_models = opex_cost_timeseries( Reduce_New_delivery_models, costs_table(15,1), costs_table(15,2) );
OPEX_Sub_Paper_Prod             = opex_cost_timeseries( Substitute_Paper, costs_table(16,1), costs_table(16,2) );
OPEX_Sub_Coated_paper_Prod      = opex_cost_timeseries( Substitute_Coated_paper, costs_table(17,1), costs_table(17,2) );
OPEX_Sub_Compostables_Prod      = opex_cost_timeseries( Substitute_Compostables, costs_table(18,1), costs_table(18,2) );
OPEX_Sub_Paper_Waste            = opex_cost_timeseries( Substitute_Paper, costs_table(19,1), costs_table(19,2) );
OPEX_Sub_Coated_paper_Waste     = opex_cost_timeseries( Substitute_Coated_paper, costs_table(20,1), costs_table(20,2) );
OPEX_Sub_Compostables_Waste     = opex_cost_timeseries( Substitute_Compostables, costs_table(21,1), costs_table(21,2) );

Econ_table(:,column_counter + 1) = OPEX_Virgin_plastic_production; %22
Econ_table(:,column_counter + 2) = OPEX_Plastic_conversion; %23
Econ_table(:,column_counter + 3) = OPEX_Formal_collection; %24
Econ_table(:,column_counter + 4) = OPEX_Informal_collection; %25
Econ_table(:,column_counter + 5) = OPEX_Formal_sorting; %26
Econ_table(:,column_counter + 6) = OPEX_Closed_loop_MR; %27
Econ_table(:,column_counter + 7) = OPEX_Open_loop_MR; %28
Econ_table(:,column_counter + 8) = OPEX_Chemical_conversion_P2P; %29
Econ_table(:,column_counter + 9) = OPEX_Chemical_conversion_P2F; %30
Econ_table(:,column_counter + 10) = OPEX_Thermal_treatment; %31
Econ_table(:,column_counter + 11) = OPEX_Engineered_landfills; %32
Econ_table(:,column_counter + 12) = OPEX_Import_sorting; %33
Econ_table(:,column_counter + 13) = OPEX_Reduce_Eliminate; %34
Econ_table(:,column_counter + 14) = OPEX_Reduce_Reuse; %35
Econ_table(:,column_counter + 15) = OPEX_Reduce_New_delivery_models; %36
Econ_table(:,column_counter + 16) = OPEX_Sub_Paper_Prod; %37
Econ_table(:,column_counter + 17) = OPEX_Sub_Coated_paper_Prod; %38
Econ_table(:,column_counter + 18) = OPEX_Sub_Compostables_Prod; %39
Econ_table(:,column_counter + 19) = OPEX_Sub_Paper_Waste; %40
Econ_table(:,column_counter + 20) = OPEX_Sub_Coated_paper_Waste; %41
Econ_table(:,column_counter + 21) = OPEX_Sub_Compostables_Waste; %42

column_counter = 42;

% CAPEX -------------------------------------------------------------------

dummy_CAPEXann_Formal_sorting_cost_rate = ones(duration,1);
dummy_CAPEX_timeseries_MC = ones(duration,1);
dummy_plastic_flow = ones(duration,1);
% dummy values for last argument to CAPEX function, when arg not needed

% CAPEX_timeseries_type=1 for all of these CAPEX flows
[CAPEXann_Virgin_plastic_production, CAPEXann_Virgin_plastic_production_cost_rate]...
    = capex_cost_timeseries( Virgin_plastic_production, costs_table(1,4),  costs_table(1,6),  costs_table(1,8),...
    costs_table(1,11),   costs_table(1,16),  CAPEX_timeseries_type(1),  dummy_CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Virgin_plastic_production  );

[CAPEXann_Plastic_conversion, CAPEXann_Plastic_conversion_cost_rate]...
    = capex_cost_timeseries( Plastic_conversion,        costs_table(2,4),  costs_table(2,6),  costs_table(2,8),...
    costs_table(2,11),   costs_table(2,16),  CAPEX_timeseries_type(2), dummy_CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Plastic_conversion  );

[CAPEXann_Formal_collection, CAPEXann_Formal_collection_cost_rate]...
    = capex_cost_timeseries( Formal_collection,         costs_table(3,4),  costs_table(3,6),  costs_table(3,8),...
    costs_table(3,11),   costs_table(3,16),  CAPEX_timeseries_type(3), dummy_CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Formal_collection  );

[CAPEXann_Informal_collection, CAPEXann_Informal_collection_cost_rate]...
    = capex_cost_timeseries( Informal_collection,       costs_table(4,4),  costs_table(4,6),  costs_table(4,8),...
    costs_table(4,11),   costs_table(4,16),  CAPEX_timeseries_type(4), dummy_CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Informal_collection  ); 

[CAPEXann_Formal_sorting, CAPEXann_Formal_sorting_cost_rate]...
    = capex_cost_timeseries( Formal_sorting,            costs_table(5,4),  costs_table(5,6),  costs_table(5,8),...
    costs_table(5,11),   costs_table(5,16),  CAPEX_timeseries_type(5), dummy_CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Formal_sorting  );

[CAPEXann_Closed_loop_MR, CAPEXann_Closed_loop_MR_cost_rate]...
    = capex_cost_timeseries( Closed_loop_MR,            costs_table(6,4),  costs_table(6,6),  costs_table(6,8),...
    costs_table(6,11),   costs_table(6,16),  CAPEX_timeseries_type(6), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Closed_loop_MR  );

[CAPEXann_Open_loop_MR, CAPEXann_Open_loop_MR_cost_rate]...
    = capex_cost_timeseries( Open_loop_MR,              costs_table(7,4),  costs_table(7,6),  costs_table(7,8),...
    costs_table(7,11),   costs_table(7,16),  CAPEX_timeseries_type(7), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Open_loop_MR  );

[CAPEXann_Chemical_conversion_P2P, CAPEXann_Chemical_conversion_P2P_cost_rate]...
    = capex_cost_timeseries( Chemical_conversion_P2P,   costs_table(8,4),  costs_table(8,6),  costs_table(8,8),...
    costs_table(8,11),   costs_table(8,16),  CAPEX_timeseries_type(8), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Total_chem_conversion  );

[CAPEXann_Chemical_conversion_P2F, CAPEXann_Chemical_conversion_P2F_cost_rate]...
    = capex_cost_timeseries( Chemical_conversion_P2F,   costs_table(9,4),  costs_table(9,6),  costs_table(9,8),...
    costs_table(9,11),   costs_table(9,16),  CAPEX_timeseries_type(9), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Total_chem_conversion  );

[CAPEXann_Thermal_treatment, CAPEXann_Thermal_treatment_cost_rate]...
    = capex_cost_timeseries( Thermal_treatment,         costs_table(10,4), costs_table(10,6), costs_table(10,8),...
    costs_table(10,11),  costs_table(10,16), CAPEX_timeseries_type(10), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Thermal_treatment  );

[CAPEXann_Engineered_landfills, CAPEXann_Engineered_landfills_cost_rate]...
    = capex_cost_timeseries( Engineered_landfills,      costs_table(11,4), costs_table(11,6), costs_table(11,8),...
    costs_table(11,11),  costs_table(11,16), CAPEX_timeseries_type(11), CAPEXann_Formal_sorting_cost_rate, dummy_CAPEX_timeseries_MC, Engineered_landfills  );

% no Importing (12)

[CAPEXann_Reduce_Eliminate, CAPEXann_Reduce_Eliminate_cost_rate]...
    = capex_cost_timeseries( Reduce_Eliminate,          costs_table(13,4), costs_table(13,6), costs_table(13,8),...
    costs_table(13,11),  costs_table(13,16), CAPEX_timeseries_type(13), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(1,:), Reduce_Eliminate  );

[CAPEXann_Reduce_Reuse, CAPEXann_Reduce_Reuse_cost_rate]...
    = capex_cost_timeseries( Reduce_Reuse,              costs_table(14,4), costs_table(14,6), costs_table(14,8),...
    costs_table(14,11),  costs_table(14,16), CAPEX_timeseries_type(14), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(2,:), Reduce_Reuse  );

[CAPEXann_Reduce_New_delivery_models, CAPEXann_Reduce_New_delivery_models_cost_rate]...
    = capex_cost_timeseries( Reduce_New_delivery_models,costs_table(15,4), costs_table(15,6), costs_table(15,8),...
    costs_table(15,11),  costs_table(15,16), CAPEX_timeseries_type(15), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(3,:), Reduce_New_delivery_models  );

[CAPEXann_Substitute_Paper, CAPEXann_Substitute_Paper_cost_rate]...
    = capex_cost_timeseries( Substitute_Paper,          costs_table(16,4), costs_table(16,6), costs_table(16,8),...
    costs_table(16,11),  costs_table(16,16), CAPEX_timeseries_type(16), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(4,:), Substitute_Paper  );

[CAPEXann_Substitute_Coated_paper, CAPEXann_Substitute_Coated_paper_cost_rate]...
    = capex_cost_timeseries( Substitute_Coated_paper,   costs_table(17,4), costs_table(17,6), costs_table(17,8),...
    costs_table(17,11),  costs_table(17,16), CAPEX_timeseries_type(17), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(5,:), Substitute_Coated_paper  );

[CAPEXann_Substitute_Compostables, CAPEXann_Substitute_Compostables_cost_rate]...
    = capex_cost_timeseries( Substitute_Compostables,   costs_table(18,4), costs_table(18,6), costs_table(18,8),...
    costs_table(18,11),  costs_table(18,16), CAPEX_timeseries_type(18), CAPEXann_Formal_sorting_cost_rate, CAPEX_timeseries(6,:), Substitute_Compostables  );

Econ_table(:,column_counter + 1) = CAPEXann_Virgin_plastic_production; %43
Econ_table(:,column_counter + 2) = CAPEXann_Plastic_conversion; %44
Econ_table(:,column_counter + 3) = CAPEXann_Formal_collection; %45
Econ_table(:,column_counter + 4) = CAPEXann_Informal_collection; %46
Econ_table(:,column_counter + 5) = CAPEXann_Formal_sorting; %47
Econ_table(:,column_counter + 6) = CAPEXann_Closed_loop_MR; %48
Econ_table(:,column_counter + 7) = CAPEXann_Open_loop_MR; %49
Econ_table(:,column_counter + 8) = CAPEXann_Chemical_conversion_P2P; %50
Econ_table(:,column_counter + 9) = CAPEXann_Chemical_conversion_P2F; %51
Econ_table(:,column_counter + 10) = CAPEXann_Thermal_treatment; %52
Econ_table(:,column_counter + 11) = CAPEXann_Engineered_landfills; %53
Econ_table(:,column_counter + 12) = CAPEXann_Reduce_Eliminate; %54
Econ_table(:,column_counter + 13) = CAPEXann_Reduce_Reuse; %55
Econ_table(:,column_counter + 14) = CAPEXann_Reduce_New_delivery_models; %56
Econ_table(:,column_counter + 15) = CAPEXann_Substitute_Paper; %57
Econ_table(:,column_counter + 16) = CAPEXann_Substitute_Coated_paper; %58
Econ_table(:,column_counter + 17) = CAPEXann_Substitute_Compostables; %59

column_counter = 59;


% REVENUES ----------------------------------------------------------------

Closed_loop_MR_price                        = prices_timeseries(1, :);
Open_loop_MR_price                          = prices_timeseries(2, :);
Chemical_conversion_P2P_price               = prices_timeseries(3, :);
Chemical_conversion_P2F_price               = prices_timeseries(4, :);
Thermal_treatment_energy_sale_price         = prices_timeseries(5, :);

Econ_table(:,column_counter + 1) = Closed_loop_MR_price';
Econ_table(:,column_counter + 2) = Open_loop_MR_price';
Econ_table(:,column_counter + 3) = Chemical_conversion_P2P_price';
Econ_table(:,column_counter + 4) = Chemical_conversion_P2F_price';
Econ_table(:,column_counter + 5) = Thermal_treatment_energy_sale_price';

column_counter = 65;

REVENUE_Closed_loop_MR                = Actual_closed_loop_MR   .* Closed_loop_MR_price';
REVENUE_Open_loop_MR                  = Actual_open_loop_MR     .* Open_loop_MR_price';
REVENUE_Chemical_conversion_P2P       = Chemical_conversion_P2P .* Chemical_conversion_P2P_price';
REVENUE_Chemical_conversion_P2F       = Chemical_conversion_P2F .* Chemical_conversion_P2F_price';
REVENUE_Thermal_treatment_energy_sale = Thermal_treatment       .* Thermal_treatment_energy_sale_price';

Econ_table(:,column_counter + 1) = REVENUE_Closed_loop_MR;
Econ_table(:,column_counter + 2) = REVENUE_Open_loop_MR;
Econ_table(:,column_counter + 3) = REVENUE_Chemical_conversion_P2P;
Econ_table(:,column_counter + 4) = REVENUE_Chemical_conversion_P2F;
Econ_table(:,column_counter + 5) = REVENUE_Thermal_treatment_energy_sale;

column_counter = 70;


% REQUIRED INVESTMENT -----------------------------------------------------
% CAPEX at year 15 x asset duration

REQ_INVESTMENT_Virgin_plastic_production =  CAPEXann_Virgin_plastic_production_cost_rate(15)      * costs_table(1,8) * zone_proportions(z);
REQ_INVESTMENT_Plastic_conversion =         CAPEXann_Plastic_conversion_cost_rate(15)             * costs_table(2,8) * zone_proportions(z);
REQ_INVESTMENT_Formal_collection =          CAPEXann_Formal_collection_cost_rate(15)              * costs_table(3,8) * zone_proportions(z);
REQ_INVESTMENT_Informal_collection =        CAPEXann_Informal_collection_cost_rate(15)            * costs_table(4,8) * zone_proportions(z);
REQ_INVESTMENT_Formal_sorting =             CAPEXann_Formal_sorting_cost_rate(15)                 * costs_table(5,8) * zone_proportions(z);
REQ_INVESTMENT_Closed_loop_MR =             CAPEXann_Closed_loop_MR_cost_rate(15)                 * costs_table(6,8) * zone_proportions(z);
REQ_INVESTMENT_Open_loop_MR =               CAPEXann_Open_loop_MR_cost_rate(15)                   * costs_table(7,8) * zone_proportions(z);
REQ_INVESTMENT_Chemical_conversion_P2P =    CAPEXann_Chemical_conversion_P2P_cost_rate(15)        * costs_table(8,8) * zone_proportions(z);
REQ_INVESTMENT_Chemical_conversion_P2F =    CAPEXann_Chemical_conversion_P2F_cost_rate(15)        * costs_table(9,8) * zone_proportions(z);
REQ_INVESTMENT_Thermal_treatment =          CAPEXann_Thermal_treatment_cost_rate(15)              * costs_table(10,8) * zone_proportions(z);
REQ_INVESTMENT_Engineered_landfills =       CAPEXann_Engineered_landfills_cost_rate(15)                               * zone_proportions(z);

REQ_INVESTMENT_Reduce_Eliminate =           CAPEXann_Reduce_Eliminate_cost_rate(15)               * costs_table(13,8) * zone_proportions(z);
REQ_INVESTMENT_Reduce_Reuse =               CAPEXann_Reduce_Reuse_cost_rate(15)                   * costs_table(14,8) * zone_proportions(z);
REQ_INVESTMENT_Reduce_New_delivery_models = CAPEXann_Reduce_New_delivery_models_cost_rate(15)     * costs_table(15,8) * zone_proportions(z);
REQ_INVESTMENT_Substitute_Paper =           CAPEXann_Substitute_Paper_cost_rate(15)               * costs_table(16,8) * zone_proportions(z);
REQ_INVESTMENT_Substitute_Coated_paper =    CAPEXann_Substitute_Coated_paper_cost_rate(15)        * costs_table(17,8) * zone_proportions(z);
REQ_INVESTMENT_Substitute_Compostables =    CAPEXann_Substitute_Compostables_cost_rate(15)        * costs_table(18,8) * zone_proportions(z);

Econ_table(:,column_counter + 1)  = REQ_INVESTMENT_Virgin_plastic_production; %71
Econ_table(:,column_counter + 2)  = REQ_INVESTMENT_Plastic_conversion; %72
Econ_table(:,column_counter + 3)  = REQ_INVESTMENT_Formal_collection; %73
Econ_table(:,column_counter + 4)  = REQ_INVESTMENT_Informal_collection; %74
Econ_table(:,column_counter + 5)  = REQ_INVESTMENT_Formal_sorting; %75
Econ_table(:,column_counter + 6)  = REQ_INVESTMENT_Closed_loop_MR; %76
Econ_table(:,column_counter + 7)  = REQ_INVESTMENT_Open_loop_MR; %77
Econ_table(:,column_counter + 8)  = REQ_INVESTMENT_Chemical_conversion_P2P; %78
Econ_table(:,column_counter + 9)  = REQ_INVESTMENT_Chemical_conversion_P2F; %79
Econ_table(:,column_counter + 10) = REQ_INVESTMENT_Thermal_treatment; %80
Econ_table(:,column_counter + 11) = REQ_INVESTMENT_Engineered_landfills; %81
Econ_table(:,column_counter + 12) = REQ_INVESTMENT_Reduce_Eliminate; %82
Econ_table(:,column_counter + 13) = REQ_INVESTMENT_Reduce_Reuse; %83
Econ_table(:,column_counter + 14) = REQ_INVESTMENT_Reduce_New_delivery_models; %84
Econ_table(:,column_counter + 15) = REQ_INVESTMENT_Substitute_Paper; %85
Econ_table(:,column_counter + 16) = REQ_INVESTMENT_Substitute_Coated_paper; %86
Econ_table(:,column_counter + 17) = REQ_INVESTMENT_Substitute_Compostables; %87

column_counter = 87;


% GHG

GHG_Virgin_plastic_production   = Virgin_plastic_production .* GHG_timeseries(1,:)';
GHG_Plastic_conversion          = Plastic_conversion        .* GHG_timeseries(2,:)';
GHG_Formal_collection           = Formal_collection         .* GHG_timeseries(3,:)';
GHG_Formal_sorting              = Formal_sorting            .* GHG_timeseries(4,:)';
GHG_Closed_loop_MR              = Closed_loop_MR            .* GHG_timeseries(5,:)';
GHG_Open_loop_MR                = Open_loop_MR              .* GHG_timeseries(6,:)';
GHG_Chemical_conversion_P2P     = Chemical_conversion_P2P   .* GHG_timeseries(7,:)';
GHG_Chemical_conversion_P2F     = Chemical_conversion_P2F   .* GHG_timeseries(8,:)';
GHG_Thermal_treatment           = Thermal_treatment         .* GHG_timeseries(9,:)';
GHG_Engineered_landfills        = Engineered_landfills      .* GHG_timeseries(10,:)';
GHG_Open_burning                = Open_burning              .* GHG_timeseries(11,:)';
GHG_Reduce_Eliminate            = Reduce_Eliminate'          .* GHG_timeseries(12,:)';
GHG_Reduce_Reuse                = Reduce_Reuse'              .* GHG_timeseries(13,:)';
GHG_Reduce_New_delivery_models  = Reduce_New_delivery_models'.* GHG_timeseries(14,:)';
GHG_Substitute_Paper            = Substitute_Paper'          .* GHG_timeseries(15,:)';
GHG_Substitute_Coated_paper     = Substitute_Coated_paper'   .* GHG_timeseries(16,:)';
GHG_Substitute_Compostables     = Substitute_Compostables'   .* GHG_timeseries(17,:)';
%}

Econ_table(:,column_counter + 1)  = GHG_Virgin_plastic_production;
Econ_table(:,column_counter + 2)  = GHG_Plastic_conversion;
Econ_table(:,column_counter + 3)  = GHG_Formal_collection;
Econ_table(:,column_counter + 4)  = GHG_Formal_sorting;
Econ_table(:,column_counter + 5)  = GHG_Closed_loop_MR;
Econ_table(:,column_counter + 6)  = GHG_Open_loop_MR;
Econ_table(:,column_counter + 7)  = GHG_Chemical_conversion_P2P;
Econ_table(:,column_counter + 8)  = GHG_Chemical_conversion_P2F;
Econ_table(:,column_counter + 9)  = GHG_Thermal_treatment;
Econ_table(:,column_counter + 10) = GHG_Engineered_landfills;
Econ_table(:,column_counter + 11) = GHG_Open_burning;
Econ_table(:,column_counter + 12) = GHG_Reduce_Eliminate;
Econ_table(:,column_counter + 13) = GHG_Reduce_Reuse;
Econ_table(:,column_counter + 14) = GHG_Reduce_New_delivery_models;
Econ_table(:,column_counter + 15) = GHG_Substitute_Paper;
Econ_table(:,column_counter + 16) = GHG_Substitute_Coated_paper;
Econ_table(:,column_counter + 17) = GHG_Substitute_Compostables;

column_counter = 104;

% JOB CREATION

JOBS_Virgin_plastic_production  = (Virgin_plastic_production/1000)      .* Jobs_timeseries(1,:)';
JOBS_Plastic_conversion         = (Plastic_conversion/1000)             .* Jobs_timeseries(2,:)';
JOBS_Formal_collection          = (Formal_collection/1000)              .* Jobs_timeseries(3,:)';
JOBS_Informal_collection        = (Informal_collection/1000)            .* Jobs_timeseries(4,:)';
JOBS_Formal_sorting             = (Formal_sorting/1000)                 .* Jobs_timeseries(5,:)';
JOBS_Closed_loop_MR             = (Actual_closed_loop_MR/1000)          .* Jobs_timeseries(6,:)';
JOBS_Open_loop_MR               = (Actual_open_loop_MR/1000)            .* Jobs_timeseries(7,:)';
JOBS_Chemical_conversion_P2P    = (Chemical_conversion_P2P/1000)        .* Jobs_timeseries(8,:)';
JOBS_Chemical_conversion_P2F    = (Chemical_conversion_P2F/1000)        .* Jobs_timeseries(9,:)';
JOBS_Thermal_treatment          = (Thermal_treatment/1000)              .* Jobs_timeseries(10,:)';
JOBS_Engineered_landfills       = (Engineered_landfills/1000)           .* Jobs_timeseries(11,:)';
JOBS_Reduce_Eliminate           = (Reduce_Eliminate/1000)'              .* Jobs_timeseries(12,:)';
JOBS_Reduce_Reuse               = (Reduce_Reuse/1000)'                  .* Jobs_timeseries(13,:)';
JOBS_Reduce_New_delivery_models = (Reduce_New_delivery_models/1000)'    .* Jobs_timeseries(14,:)';
JOBS_Substitute_Paper           = (Substitute_Paper/1000)'              .* Jobs_timeseries(15,:)';
JOBS_Substitute_Coated_paper    = (Substitute_Coated_paper/1000)'       .* Jobs_timeseries(16,:)';
JOBS_Substitute_Compostables    = (Substitute_Compostables/1000)'       .* Jobs_timeseries(17,:)';

Econ_table(:,column_counter + 1) = JOBS_Virgin_plastic_production;
Econ_table(:,column_counter + 2) = JOBS_Plastic_conversion;
Econ_table(:,column_counter + 3) = JOBS_Formal_collection;
Econ_table(:,column_counter + 4) = JOBS_Informal_collection;
Econ_table(:,column_counter + 5) = JOBS_Formal_sorting;
Econ_table(:,column_counter + 6) = JOBS_Closed_loop_MR;
Econ_table(:,column_counter + 7) = JOBS_Open_loop_MR;
Econ_table(:,column_counter + 8) = JOBS_Chemical_conversion_P2P;
Econ_table(:,column_counter + 9) = JOBS_Chemical_conversion_P2F;
Econ_table(:,column_counter + 10) = JOBS_Thermal_treatment;
Econ_table(:,column_counter + 11) = JOBS_Engineered_landfills;
Econ_table(:,column_counter + 12) = JOBS_Reduce_Eliminate;
Econ_table(:,column_counter + 13) = JOBS_Reduce_Reuse;
Econ_table(:,column_counter + 14) = JOBS_Reduce_New_delivery_models;
Econ_table(:,column_counter + 15) = JOBS_Substitute_Paper;
Econ_table(:,column_counter + 16) = JOBS_Substitute_Coated_paper;
Econ_table(:,column_counter + 17) = JOBS_Substitute_Compostables;

Econ_table_MC(1:duration,1:length(Econ_table), MCi) = Econ_table;

[Econ_table_rows, Econ_table_cols] = size(Econ_table);

end