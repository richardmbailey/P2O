function  interpolated_timeseries_value = interpolate_arrow_at_t(equation_or_timeseries, plastic_type, n_flows,...
    tspan, arrow_timeseries, t_real, imports_interaction, exports_interaction)

timeseries_to_interpolate = find(equation_or_timeseries(:,plastic_type) == 2);
n_to_interpolate = sum(equation_or_timeseries(:,plastic_type) == 2);
interpolated_timeseries_value = zeros(n_flows,1);
years = 0:max(tspan);

% Note: t_real is the 'real' full precision time, rather than the year number

row_offset = ((plastic_type - 1) * 48) + 1;


for i = 1:n_to_interpolate
    
    if t_real == max(years)
        interpolated_timeseries_value(timeseries_to_interpolate(i)) =...
            arrow_timeseries(row_offset + timeseries_to_interpolate(i), max(years));
        
    elseif t_real < 1
        
        if ismember(timeseries_to_interpolate(i), [1,2,3,4, imports_interaction, exports_interaction])
            interpolated_timeseries_value(timeseries_to_interpolate(i)) = max(...
                0,...
                arrow_timeseries(row_offset + timeseries_to_interpolate(i), 1) );
            
        else
            interpolated_timeseries_value(timeseries_to_interpolate(i)) = ...
                (interp1([1, 2],...
                [arrow_timeseries(row_offset + timeseries_to_interpolate(i), 1),...
                arrow_timeseries(row_offset + timeseries_to_interpolate(i), 2)],...
                t_real-0.0, 'linear', 'extrap'));
            if interpolated_timeseries_value(timeseries_to_interpolate(i)) < 0
                interpolated_timeseries_value(timeseries_to_interpolate(i)) = ...
                    arrow_timeseries(row_offset + timeseries_to_interpolate(i), 1);
                    %catch the case where second value is large and
                    %extrapolation gives a negative value
            end
            
        end
    else
        
        first_t_index = fix(t_real);
        second_t_index = fix(t_real) + 1;
        
        
        if ismember(timeseries_to_interpolate(i), [1,2,3,4,9,11, imports_interaction, exports_interaction])
            interpolated_timeseries_value(timeseries_to_interpolate(i)) =...
                arrow_timeseries(row_offset + timeseries_to_interpolate(i), second_t_index);
        else
            interpolated_timeseries_value(timeseries_to_interpolate(i)) = max(...
                0,...
                (interp1([first_t_index, second_t_index],...
                [arrow_timeseries(row_offset + timeseries_to_interpolate(i), first_t_index),...
                arrow_timeseries(row_offset + timeseries_to_interpolate(i), second_t_index)],...
                t_real+0.5, 'linear', 'extrap'))   );
            
            if interpolated_timeseries_value(timeseries_to_interpolate(i)) < 0
                interpolated_timeseries_value(timeseries_to_interpolate(i)) =...
                    arrow_timeseries(row_offset + timeseries_to_interpolate(i), 1);
            end
            
        end
        
    end
    
end

end