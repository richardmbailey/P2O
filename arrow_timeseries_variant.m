function variant = arrow_timeseries_variant(timeseries,...
    pedigree_categorey, pedigree_values, rnd_distribution_type, do_one_at_a_time, mc)

n = length(timeseries);
variant = zeros(1, n);

one_at_a_time_scaler = [1.0, 0.975, 1.025];


if do_one_at_a_time == false()
    
if rnd_distribution_type == 1 %Normal
    
    scaling = random('Normal', 0, pedigree_values(pedigree_categorey + 1) * 0.01);
    
    for i = 1:n
        variant(i) = timeseries(i) +  timeseries(i) * scaling ;
    end
    
elseif rnd_distribution_type == 2 % Uniform
    
    scaling = random('Uniform', -pedigree_values(pedigree_categorey + 1), pedigree_values(pedigree_categorey + 1));
    
    for i = 1:n
        variant(i) = timeseries(i) +  timeseries(i) * scaling * 0.01;
    end
    
else
end

elseif do_one_at_a_time == true()
    
    variant = timeseries * one_at_a_time_scaler(mc);
    
else
end

variant = max(0, variant);

end