function all_combinations = get_all_combinations(options)
%     all_combinations = get_all_combinations(options)
%     list all combinations of given options / Iddo Weiner
% 
%     NOTE: this functions is built for listing string combinations
% 
%     IN:
%     options = cellarray of cells. The big cellarray represents the indices of
%     the final output, and each inner cell holds all the optional strings
%     for its position. 
%     INPUT EXAMPLE: {{'a','b','c'},{'#','$'},{'1','2'}} 
%     This input means that in position1 we allow 'a' / 'b' / 'c', 
%     in position2 we allow '#' / '$' 
%     and in position3 we allow '1' / '2'. 
% 
%     OUT:
%     simply a cellarray holding all possible combinations
% 
%     HOW IT WORKS:
%     The function creates n columns, matching n positions (3 in the example above). 
%     Each one is N rows long, matching the total # of possibilities. 
%     The function then concatenates all columns to form the output, 
%     so the crucial part is building the columns proper order, 
%     to avoid receiving non-unique outputs, or (another way to look at it) -  
%     miss some combinations.
%     Columns are built on the following combinatorical rule:
%     #(reps of each unique char within a position) = ...
%                                   N / L / (product of previous options)
%     Where:
%     N = total number of options in output (simply the product of possibilities
%     per position).
%     L = number of options in position i.
%     product of previous options = multiplication of number of options in
%     previous positions. 
%     EXAMPLE: when you're in position1, this value is 1 by default. 
%     In our example (from above), in position 2, this value is 3, bcs there are 3
%     options in position 1. In position 3, the value is 6, and so on. 
%     It is easy to see that #(reps of each unique char within a position) in
%     the last position is always 1.
% 
%     OUTPUT EXAMPLE:
%     get_all_combinations({{'a','b','c'},{'#','$'},{'1','2'}})
%     ans = 
%         'a#1'
%         'a#2'
%         'a$1'
%         'a$2'
%         'b#1'
%         'b#2'
%         'b$1'
%         'b$2'
%         'c#1'
%         'c#2'
%         'c$1'
%         'c$2'

%% get parameters for setup and initiation
% setup
op_len = cellfun(@length, options);
N = prod(op_len); %total options
% initiation 
prev_prod = 1;
all_combinations = cell(N,1);

%% create output
for i = 1:(length(op_len)) %index over columns
    
    % general iteration data
    temp_data = options{i};
    L = op_len(i);
    
    % figure out how many times each digit needs to be repeated
    % this is #(reps of each unique char within a position)
    Reps = N / L / prev_prod;
    
    % create temp column elementary subunit
    temp_col = cell(0,0);
    for j = 1:L
        temp_col = vertcat(temp_col, repmat({temp_data{j}},Reps,1));
    end
    % expand it to length N
    Multiplier = N / length(temp_col);
    temp_col = repmat(temp_col, Multiplier, 1);
    
    % concat temp column with output column
    all_combinations = cellfun(@(x,y) [x y], all_combinations, temp_col,...
        'un', 0);
    
    % update prev_prod for next iteration
    prev_prod = prev_prod * L;

end

%% verify that all went well
if length(unique(all_combinations)) - N ~= 0
    error(sprintf('The function output contained non-unique rows\nCheck the input data carefully for non-unique options within a position\nother problems could be non-recognized strings specified as options'))
end
end