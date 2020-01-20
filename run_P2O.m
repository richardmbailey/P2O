clear all; clc;

% ----------- Instructions ------------------------------------------------------------------------------------------

% Copy relevant configuration files in to ~/P2O/config_files from ~/P2O/config_files_store
% (use the Archetype and Scenario numbers to specify the model run in 'Set-up choices' below)
% Set relevant choices below
% Run this script
% Results found in ~/P20/Output_files

% ----------- Archetype and Scenario numbers ------------------------------------------------------------------------

% ARCHETYPES:               SCENARIOS:
% 1. HI_Urban_batch         1. Baseline
% 2. HI_Rural_batch         2. Current commitments
% 3. UMI_Urban_batch        3. Linear Scenario
% 4. UMI_Rural_batch        4. Recycling scenario
% 5. LMI_Urban_batch        5. R&S scenario
% 6. LMI_Rural_batch        6. Scenario X
% 7. LI_Urban_batch
% 8. LI_Rural_batch

% ----------- Set-up choices ----------------------------------------------------------------------------------------

micro = false();            % running micro model: true() / false()
macro = true();             % running macro model: true() / false() - choose either macro OR micro to be true()

batch_run = false();        % running as batch (true()) or individual run (false())
archetypes_to_run = [1];    % if not running in batch mode, specify archetype to run [1,2,3,4,5,6,7,8]
scenarios_to_run =  [1];    % if not running in batch mode, specify scenario to run [1,2,3,4,5,6]
do_sensitivity_analysis = false(); %running sensitivity analysis (macro only): true() / false()

% Comment out 'micro_type' definition below if a micro batch run is being
% set up using set-up scripts which write 'micro_type_to_run.mat' to the code directory
micro_type = 'TIRE'; % set micro source type to match chosen config files: 'PCP', 'TIRE', 'SYNTHTEX', 'PELLETS'

% Comment out 'number_MC_iterations' below to revert to reading this value in from configuration files
number_MC_iterations = 1;

% ----------- Sensitivity analysis set-up ---------------------------------------------------------------------------

% These are only used during sensitivity runs, ignored otherwise
model_year_to_analyse = 2;  %from 1 to 'duration' (25 yrs in downloaded version)
flow_numbers_for_analysis = [5]; %[5, 17, 18, 20, 21, 22, 28, 7, 9, 10, 11, 13, 14, 16, 24, 25 ];
n_variants = 500;           %for random numbers stored; needs to be at least as high as number of MC runs

% ----------- Additional parameters ---------------------------------------------------------------------------------

MC_distribution_type = 2;   % Sampling distribution type for MC analysis (overides basic_info.csv)
                            % 1=Normal; 2=Uniform
rand_scheme = 2;            % random numbers: 1=generated on-the-fly, 2=read from list created [default = 2]

% (parameters below do not need to be changed unless a new model is defined)
duration = 25;              % 25 is correct value unless new model is defined, then update here, and in basic_info.csv
macro_exports_flow = 17;    % 17 is correct value (macro-plastic only) unless new model is defined, then update here
n_interaction_parameters = 17; % number of columns in each parameter table
n_interaction_arrows = 48;  % total rows per plastic type in 'Imported_data' tab
datapoints_per_year = 250;  % speed/accuracy trade-off; 250 is a good compromise
Econ_table_n_cols = 121;    % Included for flexibility: 121 is correct value unless new model is defined
sensitivity_matrix_columns = 49; % Included for flexibility: 49 is correct value unless new model is defined

% -------------------------------------------------------------------------------------------------------------------





if batch_run
    arch_sen_to_run = csvread('arch_sen_to_run.csv');
    batch_archetype_list = arch_sen_to_run(1);
    batch_scenario_list = arch_sen_to_run(2);
else
    batch_archetype_list = archetypes_to_run;
    batch_scenario_list  = scenarios_to_run;
end

if macro
    zones_to_run = [1,2];
    n_flows = 44;
    
elseif micro
    zones_to_run = 1; % (only 1 zone in micro set-up)
    zone_proportions = 1; %(set in config files for macro)
    
    % Default: set micro plastic type from file (created by batch set-up scripts)
    check_exists = exist('micro_type', 'var');
    if check_exists == 0
        load('micro_type_to_run.mat');
    end
    
    config_prefix = '/config_files/';
    ff = fullfile(pwd, [config_prefix, 'stock_stock_interactions.csv']);
    stock_stock_interactions = csvread(ff);
    n_flows = sum(sum(stock_stock_interactions>0));
    
else
end


if do_sensitivity_analysis
    do_one_at_a_time = true();
    analysis_year = model_year_to_analyse;
    x_analysis_flow_list = flow_numbers_for_analysis;  

    n_sensitivity_runs = length(x_analysis_flow_list);
    sensitivity_matrix_flow_table = zeros(n_sensitivity_runs, sensitivity_matrix_columns);
    sensitivity_matrix_econ_table = zeros(n_sensitivity_runs, Econ_table_n_cols + 3);
    
else
    n_sensitivity_runs = 1;
    do_one_at_a_time = false();
end

micro_exports_flow = 0; %no exports in micro model

if macro
    n_plastic_types = 3;
    exports_interaction = macro_exports_flow;
end

if micro
    n_plastic_types = 1;
    exports_interaction = micro_exports_flow;
end

rng(batch_scenario_list(1));

if MC_distribution_type == 1 %Normal
    rand_list = random('Normal', 0, 1, (n_variants * 500), 1);
elseif MC_distribution_type == 2 % Uniform
    rand_list = random('Uniform',-1, 1, (n_variants * 500), 1);
else
    fprintf('\n\n>> ERROR: MC distribution type not recognized <<\n\n');
end

n_batch_archetypes = length(batch_archetype_list);
n_batch_scenarios  = length(batch_scenario_list);
n_zones = length(zones_to_run);
gros_plastic_demand = zeros(n_plastic_types, duration, n_zones);
baseline_individual_flow = zeros(n_flows, n_plastic_types);

global rand_num_number;
rand_num_number = 0;

for batch_archetype_i = 1:n_batch_archetypes
    
    for batch_scenario_i = 1:n_batch_scenarios
        
        % set archetype
        archetype_to_run = batch_archetype_list(batch_archetype_i);
        
        % set scenario
        n_scenarios = 1; %to run in each loop
        scenarios_to_run = batch_scenario_list(batch_scenario_i);
        
        % run model
        for sr = 1: n_sensitivity_runs
            
            if macro
                % change pwd to config_files directory
                cd config_files
                
                FileName='archetype_names.txt';
                fileID = fopen(FileName,'r');formatSpec = '%s';
                S_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
                archetype_names = vertcat(S_names{:});
                
                FileName='scenario_names.txt';
                fileID = fopen(FileName,'r');formatSpec = '%s';
                S_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
                scenario_names = vertcat(S_names{:});
                
                FileName='zone_names.txt';
                fileID = fopen(FileName,'r');formatSpec = '%s';
                S_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
                zone_names = vertcat(S_names{:});
                
                FileName='plastic_names.txt';
                fileID = fopen(FileName,'r');formatSpec = '%s';
                S_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
                plastic_names = vertcat(S_names{:});
                
                cd ..
                
            end
            
            if do_sensitivity_analysis
                x_analysis_flow = x_analysis_flow_list(sr);
            end
            
            
            %loop model runs over scenarios and zones
            for sx = 1:n_scenarios
                
                production_rates_table = zeros(duration, n_plastic_types, n_zones);
                s = scenarios_to_run(sx);
                
                if do_sensitivity_analysis
                    n_zones = 1;
                    zones_to_run = 1;
                end
                
                for z = 1:n_zones
                    
                    if macro
                        % needed: e.g. HI_Urban_BAUI_Z2_
                        config_prefix =['/config_files/', archetype_names{archetype_to_run}, '_', scenario_names{s},...
                            '_Z', num2str(z), '_'];
                        
                        % needed: e.g. HI_Urban_R&Sscenario_Plastic_
                        config_prefix2 =['/config_files/', archetype_names{archetype_to_run}, '_', scenario_names{s},...
                            '_Plastic_'];
                        %for results
                        results_filename_prefix = ['A',num2str(archetype_to_run), '_S', num2str(s),...
                            '_Z', num2str(z)];
                    elseif micro
                        config_prefix = '/config_files/';
                        config_prefix2 = '/config_files/';
                        results_filename_prefix = ['A',num2str(archetype_to_run), '_S', num2str(s),...
                            '_Z', num2str(z)];
                    end
                    
                    % RUN MODEL
                    
                    % Run model for one year to create baseline conditions
                    baseline_run = 1;
                    fprintf('\n\n-------- P2O Ocean Plastics Model v1.0.0 --------\n\n\n');
                    Ocean_plastics_v1_0_0;
                    
                    % Determine flow and mass capacities
                    end_first_yr_mass = MC_final_masses;
                    baseline_mass = zeros(1, n_stocks, n_plastic_types);
                    initial_mass_from_baseline_run = zeros(1, n_stocks, n_plastic_types);
                    
                    for k = 1:n_plastic_types
                        
                        for i = 1:n_stocks
                            
                            if isempty( intersect(i,finite_sinks) ) == false() % it's a finite sink
                                
                                baseline_mass(1,i,k) = ( Mass(end,i,k) - Mass(round(datapoints_per_year * 0.98, 0),i,k) ) * 50;
                                initial_mass_from_baseline_run(i, k) = 0;
                                
                            elseif isempty( intersect(i,finite_sinks) ) == true() % it's not a finite sink
                                
                                baseline_mass(1,i,k) = Mass(round(datapoints_per_year * 0.5, 0),i,k); %take mid-year flow as baseline
                                initial_mass_from_baseline_run(i, k) = Mass(round(datapoints_per_year * 1.0, 0),i,k);
                            end
                        end
                    end
                    
                    baseline_flow = zeros(1, n_stocks, n_plastic_types);
                      
                    for k = 1:n_plastic_types
                        for i = 1:n_stocks
                            baseline_flow(1,i,k) = sum(MC_final_flows(1, out_interaction_number(i,1:n_outflows(i)),k));
                        end
                    end
                    
                    %baseline individual flows
                    for k = 1:n_plastic_types
                        baseline_individual_flow(:,k) = interaction_flow( round(datapoints_per_year * 0.5, 0), : ,k);
                    end
                                       
                    % Run model for the full duration
                    baseline_run = 0;
                    fprintf('\n\n');
                    Ocean_plastics_v1_0_0;
                    
                end  % z: zones
                
            end % s: scenarios
            
            % Store sensitivity data in matrix
            if do_sensitivity_analysis == true()
                sensitivity_matrix_flow_table(sr, :) = sensitivity_matrix_flow(1,:);
                sensitivity_matrix_econ_table(sr, :) = sensitivity_matrix_econ(1,:);
            end
        end % sr: sensitivity runs
        
        
        % Save sensitivity data to file
        if do_sensitivity_analysis && macro    
        
            filename = ['A',num2str(archetype_to_run), '_S', num2str(scenarios_to_run),'_Sensitivity_Flow.csv'];
            file_path = fullfile(pwd, 'Output_files', filename);
            csvwrite(file_path, sensitivity_matrix_flow_table );
            
            filename = ['A',num2str(archetype_to_run), '_S', num2str(scenarios_to_run),'_Sensitivity_Econ.csv'];
            file_path = fullfile(pwd, 'Output_files', filename);
            csvwrite(file_path, sensitivity_matrix_econ_table );
            
        end
        
        % Save workspace for current archetype run
        if do_sensitivity_analysis == false()
            
            filename = ['A',num2str(archetype_to_run), '_S', num2str(s),...
                '_All_variables'];
            % 'A1_S1_All_variables.mat'
            full_path = fullfile(pwd, 'Output_files', filename);
            save(full_path);
            
        end
        
        
        % Accumulate data at Archetype level
        if do_sensitivity_analysis == false()
            
            archetype_econ_table = zeros(duration, Econ_table_n_cols, n_zones);
            archetype_econ_SD_table = zeros(duration, Econ_table_n_cols, n_zones);
            archetype_scenario_econ_table = zeros(duration, Econ_table_n_cols);
            archetype_flows_table = zeros(duration, n_interactions, n_zones);
            archetype_flows_SD_table = zeros(duration, n_interactions, n_zones);
            archetype_mass_table = zeros(duration, n_stocks, n_zones);
            archetype_mass_SD_table = zeros(duration, n_stocks, n_zones);
            
            for s = 1:n_scenarios
                
                sen = scenarios_to_run(s);
                
                for z = 1: n_zones
                    
                    % Summed flows
                    % file to read : e.g. A1_S1_Z2_summed_flows
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                        '_Z', num2str(z), '_summed_flows.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    summed_flows_temp = csvread(file_path);
                    archetype_flows_table(:,:,z) = summed_flows_temp;
                    
                    % Summed flows stdev
                    % file to read : e.g. A1_S1_Z2_summed_flows_SD
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                        '_Z', num2str(z), '_summed_flows_SD.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    summed_flows_SD_temp = csvread(file_path);
                    archetype_flows_SD_table(:,:,z) = summed_flows_SD_temp;
                    
                    % Summed masses
                    % file to read : e.g. A1_S1_Z2_summed_mass
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                        '_Z', num2str(z), '_summed_mass.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    summed_mass_temp = csvread(file_path);
                    archetype_mass_table(:,:,z) = summed_mass_temp;
                    
                    % Summed flows stdev
                    % file to read : e.g. A1_S1_Z2_summed_mass_SD
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                        '_Z', num2str(z), '_summed_mass_SD.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    summed_mass_SD_temp = csvread(file_path);
                    archetype_mass_SD_table(:,:,z) = summed_mass_SD_temp;
                    
                    if macro
                        % Econ
                        % file to read : e.g. A1_S1_Z3_ECON
                        filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                            '_Z', num2str(z), '_ECON.csv'];
                        file_path = fullfile(pwd, 'Output_files', filename);
                        econ_temp = csvread(file_path);
                        archetype_econ_table(:,:,z) = econ_temp;
                        
                        % Econ stdev
                        % file to read : e.g. A1_S1_Z3_ECON_SD
                        filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),...
                            '_Z', num2str(z), '_ECON_SD.csv'];
                        file_path = fullfile(pwd, 'Output_files', filename);
                        econ_SD_temp = csvread(file_path);
                        archetype_econ_SD_table(:,:,z) = econ_SD_temp;
                    end
                    
                end
                
                archetype_scenario_econ_table = sum(archetype_econ_table, 3);
                if macro
                    for i_sum = 1:Econ_table_n_cols
                        archetype_scenario_econ_SD_table(:, i_sum) = combine_errors( archetype_econ_SD_table(:,i_sum,1),...
                            archetype_econ_SD_table(:,i_sum,2) );
                    end
                end
                
                archetype_scenario_flows_table = sum(archetype_flows_table, 3);
                for i_sum = 1:n_flows
                    if macro
                        archetype_scenario_flows_SD_table(:, i_sum) = combine_errors( archetype_flows_SD_table(:,i_sum,1),...
                            archetype_flows_SD_table(:,i_sum,2) );
                    elseif micro
                        archetype_scenario_flows_SD_table(:, i_sum) = archetype_flows_SD_table(:,i_sum);
                    end
                    
                end
                
                archetype_scenario_mass_table = sum(archetype_mass_table, 3);
                for i_sum = 1:n_stocks
                    if macro
                        archetype_scenario_mass_SD_table(:, i_sum) = combine_errors( archetype_mass_SD_table(:,i_sum,1),...
                            archetype_mass_SD_table(:,i_sum,2) );
                    elseif micro
                        archetype_scenario_mass_SD_table(:, i_sum) = archetype_mass_SD_table(:,i_sum);
                    end
                end
                
                if macro
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_ECON.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    csvwrite(file_path, archetype_scenario_econ_table );
                    
                    filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_ECON_SD.csv'];
                    file_path = fullfile(pwd, 'Output_files', filename);
                    csvwrite(file_path, archetype_scenario_econ_SD_table );
                end
                
                filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_FLOWS.csv'];
                file_path = fullfile(pwd, 'Output_files', filename);
                csvwrite(file_path, archetype_scenario_flows_table );
                
                filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_FLOWS_SD.csv'];
                file_path = fullfile(pwd, 'Output_files', filename);
                csvwrite(file_path, archetype_scenario_flows_SD_table );
                
                filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_MASS.csv'];
                file_path = fullfile(pwd, 'Output_files', filename);
                csvwrite(file_path, archetype_scenario_mass_table );
                
                filename = ['A',num2str(archetype_to_run), '_S', num2str(sen),'_MASS_SD.csv'];
                file_path = fullfile(pwd, 'Output_files', filename);
                csvwrite(file_path, archetype_scenario_mass_SD_table );
                
            end
            
        end
        
        
        if do_sensitivity_analysis
            
            filename = ['A',num2str(archetype_to_run), '_S', num2str(s),...
                '_All_variables_Sensitivity'];
            
            full_path = fullfile(pwd, 'Output_files', filename);
            save(full_path);
            
        end
        
        
    end % Scenario batch loop
    
end %Archetype batch loop


% For running batch jobs, write empty 'CONFIRM...' file to internally confirm run completed successfully
if batch_run
    confirmation_name = ['CONFIRM_a',num2str(batch_archetype_list),'_s',num2str(batch_scenario_list), '.csv'];
    file_path = fullfile(pwd, 'Output_files', confirmation_name);
    csvwrite(file_path, NaN);
end

