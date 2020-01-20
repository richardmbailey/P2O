function sampled_version = MC_sample(rnd_distribution_type,...
    mean_val, error_val, pedigree_values)


n = length(mean_val);
sampled_version = zeros(n, 1);

if rnd_distribution_type == 1 %Normal
    
    for i = 1:n
        stdev = mean_val(i) * pedigree_values(error_val(i)+1) * 0.01;
        sampled_version(i) = random('Normal', mean_val(i), stdev);
        if mean_val(i)>=0
            sampled_version(i) = max(0, sampled_version(i) );
        end
    end
    
elseif rnd_distribution_type == 2 % Uniform
    
    for i = 1:n
        
        if mean_val(i) <= 0
            lower_val = mean_val(i)-(mean_val(i) * pedigree_values(error_val(i)+1) * 0.01);
        else
            lower_val = max(0, mean_val(i)-(mean_val(i) * pedigree_values(error_val(i)+1) * 0.01) );
        end
        
        sampled_version(i) = random('Uniform', ...
            lower_val,...
            mean_val(i)+(mean_val(i) .* pedigree_values(error_val(i)+1) * 0.01)...
            );
    end
    
else
end

end