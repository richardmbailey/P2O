function total_error = combine_errors(varargin)

% detect if a timeseries - assuming all timeseries are the same length
if numel(varargin{1}) > 1 && nargin > 1
    is_timeseries = true();
    
elseif numel(varargin{1}) == 1 && nargin > 1
    is_timeseries = false();

else
    disp('>>> Problem with combining uncertainties - only 1 in list? <<<');
end


if is_timeseries
    
    squares_table = zeros( numel(varargin{1}), nargin-1 );
    
    for i=1:nargin
        z=varargin{i};
        z_as_col = reshape(z, numel(z), 1);
        squares_table(:,i) = z_as_col.^2;
    end
    
    sum_squares = sum(squares_table,2);
    total_error = sqrt(sum_squares);
        
elseif ~is_timeseries
    
    z2 = zeros(nargin,1);
    
    for i=1:nargin
        z=varargin{i};
        assert(z >= 0 && isreal(z),'>>> Check errors - all errors must be positive real numbers <<<');
        z2(i) = z^2;
    end
    
    sum_squares = sum(z2);
    total_error = sqrt(sum_squares);
    
else
end
