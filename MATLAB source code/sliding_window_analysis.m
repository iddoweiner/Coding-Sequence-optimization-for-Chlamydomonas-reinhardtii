function [All_Profiles,mean_Profile] = ...
    sliding_window_analysis(Data,window_size,your_func,number_of_wanted_output,varargin)
% [All_Profiles,mean_Profile] = sliding_window_analysis(Data,window_size,your_func,number_of_wanted_output,varargin)
% Data = cell of sequences on which some data is to be calculated
% window_size = self explanitory
% your_func = calls the calculated function (e.g. *@rnafold*)
% number_of_wanted_output = 1 / 2 / 3 - which of the outputs do you want?
% additional: 
% (...,'figure',0) will supress the figure from appearing
% ('threshold',0.75) will plot data until less than 75% of the rows are
% included in the mean profile. Assuming the data has differnt lengths.
% Default here is 0.5
% (...,'Y_label','FE') will write *FE* on the y axis label

%% create output matrices (columns as maximum length, rows as number of seqs)
Rows = size(Data,1);
Cols = max(cellfun(@length,Data)) - window_size +1;
All_Profiles = zeros(Rows,Cols);
All_Profiles(:) = nan;

%% fill matrix with 2 loops - i for rows (big) and j for cols (quick)
for i=1:Rows
    % point at 1 row
    temp_row = Data{i,1};
    for j = 1:(length(temp_row)-window_size+1)
        % calc whatever you want to calculate for this window
        if number_of_wanted_output==1
            [All_Profiles(i,j)] = your_func( temp_row( j:(j+window_size-1) ) );
        elseif number_of_wanted_output==2
            [~,All_Profiles(i,j)] = your_func( temp_row( j:(j+window_size-1) ) );
        elseif number_of_wanted_output==3
            [~,~,All_Profiles(i,j)] = your_func( temp_row( j:(j+window_size-1) ) );
        end
    end 
    
    % report progress
    disp(['finished row ' num2str(i)])
end

% fill mean matrix
mean_Profile = nanmean(All_Profiles);

%% get how many genes were calculated per position
% get how many nans per position
nan_per_position = sum( isnan(All_Profiles) );
pct_relevant_genes = (Rows - nan_per_position) / Rows;

% define stopping criteria
if find(strcmpi(varargin,'threshold'))
    threshold = varargin{find(strcmpi(varargin,'threshold'))+1};
else
    threshold = 0.5;
end

last_X = find(pct_relevant_genes<=threshold,1,'first');


%% plot all data

%look for ylabel argument
if find(strcmpi(varargin,'Y_label'))
    Y_label = varargin{find(strcmpi(varargin,'Y_label'))+1};
else
    Y_label = 'Feature';
end

%check if requested and plot unless requested not to
if find(strcmpi(varargin,'figure')) %if something said about figure
    if varargin{find(strcmpi(varargin,'figure'))+1}==1 %check if it was asked for
        plot_profile(last_X,mean_Profile,window_size,Y_label) %and plot
    end %and if it wasn't - just don't show
    
else %if nothing said about figure
    plot_profile(last_X,mean_Profile,window_size,Y_label) %simply plot

end

end


%% aid func
function plot_profile(last_X,mean_Profile,window_size,Y_label)

        figure
        X_data = [1:last_X];
        plot(X_data,mean_Profile(1:last_X))
        title('Mean profile')
        ylabel(Y_label)
        xlabel('Position of window start')
        legend(['profile for window size ' num2str(window_size)])
        set(gca,'fontsize',18)

end