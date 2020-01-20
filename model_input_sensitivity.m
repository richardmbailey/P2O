clc;

file_location = '/Users/richardbailey/Documents/GitHub/Plastics_model_ARCHIVE/Report_and_Science_paper_Jan_8_2020/Macro/OUTPUTS_files_MACRO';

input_rows = [391, 392, 436, 437, 481, 482];
OL_row = 556;
years = 2;
data_store = zeros(8, 6, numel(years) );

for k = years
    fprintf('Year %d\n ', k);
    
    for i = 1:8
        fprintf('arch %d: ', i);
        
        for j = 1:6
            
            fprintf('%d ', j);
            
            %load data
            filename = [ 'OUTPUT_A',num2str(i), '_S', num2str(j), '.csv' ];
            
            full_path = fullfile(file_location, filename);
            data = csvread(full_path);
            
            %total ocean leakage
            total_ocean_leakage = sum( data(input_rows, k) );
            OL = data( OL_row, k );
            
            %sensitivity
            s = OL / total_ocean_leakage;
            
            %store sensitivity
            data_store(i,j,k) = s;
        end
        fprintf('\n');
    end
    
end  

%% plots
arch = 1; 
sen = 1;
dat = reshape( data_store(arch,sen,:), 1, numel(data_store(arch,sen,:)) );
plot(dat, '-o');


