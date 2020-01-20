function sampled_version = MC_sample_list(rnd_distribution_type,...
    mean_val, error_val, pedigree_values, rand_list)


global rand_num_number

n = length(mean_val);
sampled_version = zeros(n, 1);

if rnd_distribution_type == 1 %Normal
    
    for i = 1:n
        
        rand_num_number = rand_num_number + 1;

        stdev = mean_val(i) * pedigree_values(error_val(i)+1) * 0.01;

        sampled_version(i) = mean_val(i) + (rand_list(rand_num_number) * stdev);
        
        if mean_val(i)>=0
            sampled_version(i) = max(0, sampled_version(i) );
        end
    end
    
elseif rnd_distribution_type == 2 % Uniform
    
    for i = 1:n
        rand_num_number = rand_num_number + 1;

        if mean_val(i) <= 0
            lower_val = mean_val(i)-(mean_val(i) * pedigree_values(error_val(i)+1) * 0.01);
        else
            lower_val = max(0, mean_val(i)-(mean_val(i) * pedigree_values(error_val(i)+1) * 0.01) );
        end
        
        upper_val = mean_val(i)+(mean_val(i) .* pedigree_values(error_val(i)+1) * 0.01);
        
        %sampled_version(i) = lower_val + (upper_val - lower_val) * rand_list(rand_num_number); %works for rnd~U[0,1]
        sampled_version(i) = mean_val(i) + (pedigree_values(error_val(i)+1) * 0.01) * mean_val(i) * rand_list(rand_num_number); %works for rnd~U[-1,1]
        
    end
    
else
end

end