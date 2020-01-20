function variant = arrow_timeseries_variant_list(timeseries,...
    pedigree_categorey, pedigree_values, rnd_distribution_type, do_one_at_a_time, mc, rand_list)


global rand_num_number
rand_num_number = rand_num_number + 1;

%fprintf('arrow: #%d\trnd: %f\n', rand_num_number, rand_list(rand_num_number));

n = length(timeseries);
variant = zeros(1, n);

one_at_a_time_scaler = [1.0, 0.975, 1.025];


if do_one_at_a_time == false()
    
if rnd_distribution_type == 1 %Normal (rel errors defined 0 ± 1)
    
    scaling = rand_list(rand_num_number) * (pedigree_values(pedigree_categorey + 1) * 0.01);
    
    for i = 1:n
        variant(i) = timeseries(i) +  timeseries(i) * scaling ; 
    end
    
elseif rnd_distribution_type == 2 % Uniform (rel errors defined from -1 to +1)
    
    scaling = rand_list(rand_num_number) * (pedigree_values(pedigree_categorey + 1) * 0.01);
    
    for i = 1:n
        variant(i) = timeseries(i) +  timeseries(i) * scaling; 
    end
    
else
end

elseif do_one_at_a_time == true()
    
    variant = timeseries * one_at_a_time_scaler(mc);
    
else
end

variant = max(0, variant);

end