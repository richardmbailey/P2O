
function variant = stdev_scaling(timeseries, timeseries_stdev, MC_variant_number)

if MC_variant_number > 1
    
    rnd_scaler = normrnd(0,1);
    variant = timeseries + rnd_scaler .* timeseries_stdev;
    
    variant = max(0, variant);
    
else
    
    variant = timeseries;
    
end

