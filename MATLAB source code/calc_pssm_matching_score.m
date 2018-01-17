function score = calc_pssm_matching_score(pssm, seq, varargin)
% score = calc_pssm_matching_score(pssm, seq, varargin)
% pssm = probability matrix / table (heavy time-wise)
% seq  = DNA seq
% var:
% (...,'log')
% the function gives a score to each position according to the pssm
% then geomeans over all position scores
% NOTE: take care of zeros! add epsilon in that case

%%
% transform to matrix if nec 
% (*this is a costly step time-wise, better do this outside the function once)
if istable(pssm)
    pssm = table2array(pssm);
end
% transform to log if asked to
if sum(find(strcmpi(varargin,'log')))
    pssm = log(pssm);
end
% general stuff
seq = upper(seq); %protect from lower case errors
epsilon = 0.001;  %define eps to be added to zeros
Cols = size(pssm,2); %get # of colunms in seq

%% match
scoring_vector = zeros(1,Cols);
for idx = 1:Cols
   temp_char = seq(idx);
   % build scenarios
   switch temp_char
       case 'A'
           scoring_vector(1,idx) = pssm(1,idx);
       case 'C'
           scoring_vector(1,idx) = pssm(2,idx);
       case 'G'
           scoring_vector(1,idx) = pssm(3,idx);
       case {'T', 'U'}
           scoring_vector(1,idx) = pssm(4,idx);
       otherwise
           % warning(['Fuck! letter number ' num2str(idx) ' is not A/C/G/T'])
           scoring_vector(1,idx) = nan;
   end %of switch case for the i nt
end %of scoring for this seq

%% finalize
% get rid of zeros by adding epsilon, this enables calculation of geomean 
% (if the input comes from HOMER - this isn't supposed to happen)
scoring_vector(scoring_vector==0) = epsilon;
% calc mean
if sum(find(strcmpi(varargin,'log')))
    % math mean on logs if required
    score = mean(scoring_vector);
    score = exp(score);
else %this is the default
    score = geomean(scoring_vector);
end

end