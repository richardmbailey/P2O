function total_error = combine_errors2(varargin)

% error combination for addition/subtraction

is_timeseries = false();
is_single_values = false();
is_matrix = false();
is_single_timeseries = false();
is_single_value = false();

[n_rows, n_cols] = size(varargin{1});

% Assumes all timeseries are the same length if multiple timeseries or
% matrix (assumes matrix data in multiple rows)

% detect if a timeseries
if numel(varargin{1}) > 1 && nargin > 1 && n_rows == 1
    is_timeseries = true();
    %disp('multiple time series inputted');
    
elseif numel(varargin{1}) == 1 && nargin > 1 && n_rows == 1
    is_single_values = true();
    %disp('multiple single values inputted');
    
elseif numel(varargin{1}) > 1 && nargin == 1 && n_rows > 1
    is_matrix = true();
    %disp('matrix inputted');
    
elseif numel(varargin{1}) > 1 && nargin == 1 && n_rows == 1
    is_single_timeseries = true();
    %disp('single timeseries inputted');
    
elseif numel(varargin{1}) == 1 && nargin == 1 && n_rows == 1
    is_single_value = true();
    %disp('single value inputted');
else
    disp('>>> Problem combining uncertainties....');
    disp('>>> Check format of input data to "combine_errors2.m"');
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
        
elseif is_single_values
    
    z2 = zeros(nargin,1);
    
    for i=1:nargin
        z=varargin{i};
        assert(z >= 0 && isreal(z),'>>> Check errors - all errors must be positive real numbers <<<');
        z2(i) = z^2;
    end
    
    sum_squares = sum(z2);
    total_error = sqrt(sum_squares);
    
elseif is_matrix % (assumes matrix data in columns)
    
    squares_table = zeros( n_rows, n_cols );
    
    for i=1:n_rows
        z=varargin{1};
        squares_table = z.^2;
    end
    
    sum_squares = sum(squares_table,1);
    total_error = sqrt(sum_squares);

elseif is_single_timeseries %Just return the input as nothing to combine
    total_error = varargin{1};
    
elseif is_single_value %Just return the input as nothing to combine
    total_error = varargin{1};
    
else
end

end

