function flow_rate = calculate_flow_rate(a,b,c,d,t,function_type)


%switch function_type
    
if (function_type == 1)
        %linear dependence
        flow_rate = a + t * b;
end
    

if (function_type == 2)
        %exponential dependence
        flow_rate = a * exp(t * b);
end
        

if (function_type == 3)
        %compound effect
        if a==0
            flow_rate = 0;
        elseif a~=0
        flow_rate = a * ((1 + b)^t);
        else
        end
end


if (function_type == 4)
        %compound effect with minumum value set to c
        if a==0
            flow_rate = 0;
        elseif a~=0
        
        flow_rate = max( a * ((1 + b)^t), c);
        else
        end
end


if (function_type == 5)
        %compound effect with limits of [0,1]
        if a==0
            flow_rate = 0;
        elseif a~=0
            test_val = a * ((1 + b)^t) ;
            if test_val > 1
                flow_rate = 1;
            elseif test_val < 0
                flow_rate = 0;
            else
                flow_rate = test_val;
            end
        else
        end
end


        
    % more as necessary...
    c; %#ok<VUNUS>
    d; %#ok<VUNUS>
    
    
end

