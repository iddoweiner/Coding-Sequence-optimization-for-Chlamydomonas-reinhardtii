function [clean_seq, mismatches] = erase_splicing(seq, data_file)
%     seq = to be cleaned from splicing
%     data_file = name of file to load
%     The file should contain UTR5_pssm, UTR3_pssm, UTR5_threshold and UTR3_threshold
%     The function goes over all windows, finds problematic ones, checks all options
%     for synonymous codons, picks the one with score under threshold and closest
%     to original in sequence. After fixing 5' and then 3', the function rechecks 5'
%     to see that the 3' fixing didin't recreate a 5' signal, and so on..
%     If this casues an endless loop - thee will starta a shift and less similar 
%     sequences will be taken, until a double valid occurs
%
% go over problematic positions, get all combs, find best one:
% 1. expand range of problematic window to match a complete RF
% 2. extract relevant DNA sequence, translate to AA
% 3. get all combinations for this sequence
% 4. calc PSSM for the 1 *relevant window* in each combination
% 5. throw away all possibilities above threshold
% 6. rank according to similarity to original DNA seq
% 7. Choose highest
% * if iterations is high, start selecting smaller ranks
% 8. Do 1-7 for 3' signal too
% if both 5' and 3' return empty, break out of while loop
% * gather analytics at the end: how many mismatches were created between
% the original input sequence and the output given sequence?

%% initiation
load(data_file)
working_seq = seq;
iterations = 1;

ss5 = sp_sg5; ss3 = sp_sg3;
ss5_threshold = Thres5;
ss3_threshold = Thres3;

%% work    
% detection 5'
problematic_5 = find_problematic_positions(working_seq, ss5, ss5_threshold);
if ~isempty(problematic_5) %i.e. if there are problems
    % solve them
    working_seq = fix_problems(working_seq, problematic_5, ss5, ss5_threshold);
end
% show progress
disp(['cleaned ' num2str(length(problematic_5)) ' problems with 5'' signal'])

% detection 3'
problematic_3 = find_problematic_positions(working_seq, ss3, ss3_threshold);
if ~isempty(problematic_3) %i.e. if there are problems
    % solve the problems
    working_seq = fix_problems(working_seq, problematic_3, ss3, ss3_threshold);
elseif isempty(problematic_3)
    go_3 = 0;
end
% show progress
disp(['cleaned ' num2str(length(problematic_3)) ' problems with 5'' signal'])

%% chcek how many mismatches with original seq and give output
mismatches = sum(seq ~= working_seq);
% report
disp(['splicing deletion resulted in ' num2str(1 - (mismatches/length(seq))) ' similarity to original sequence'])
clean_seq = working_seq;

end

%% helper functions:
% 1. locate problematic windows in a sequence
% (called from main)
function problematic_positions = find_problematic_positions(seq, PSSM, threshold)
    scores = nan(length(seq) - size(PSSM,2) +1, 1);
    for i = 1:(length(scores))
        scores(i) = calc_pssm_matching_score( PSSM, seq(i : i+size(PSSM,2)-1), 'log' );
    end
    problematic_positions = find(scores > threshold);
end

%% 2. synonymously mutate a problematic window to solve the problem
% (called from main)
function new_seq = fix_problems(seq, prob_window_idx, PSSM, threshold)

for i = 1:(length(prob_window_idx))
    
    % get general window indices
    temp_cords(1) = prob_window_idx(i);
    temp_cords(2) = temp_cords(1) + size(PSSM,2) - 1;

    % adjust to reading frame (handle 1st index)
    beginning_rem = rem(temp_cords(1),3);
    switch beginning_rem
        case 1 %already in frame
            beg_idx = temp_cords(1);
            five_dent = 0;
        case 2 %xepand by 1
            beg_idx = temp_cords(1) - 1;
            five_dent = 1;
        case 0 %expand by 2
            beg_idx = temp_cords(1) - 2;
            five_dent = 2;
    end

    % adjust to reading frame (handle last index)
    end_rem = rem(temp_cords(2),3);
    switch end_rem
        case 0 %perfect, it's already in frame
            end_idx = temp_cords(2);
            three_dent = 0;
        case 1 %expand by 2
            end_idx = temp_cords(2) + 2;
            three_dent = 2;
        case 2 %expand by 1
            end_idx = temp_cords(2) + 1;
            three_dent = 1;
    end

    % extract the RF
    original_window_DNA = seq(beg_idx : end_idx);
    window_AA = nt2aa(original_window_DNA,'AlternativeStartCodons',0);
    % generate all possibilities
    allcombs = generate_all_codon_combinations(window_AA);
    % check score for each option:
    allcombs_pssm_score = nan(size(allcombs));
    for comb = 1:(length(allcombs))
        temp_comb = allcombs{comb};
        % 1. trim to look just at relevant window
        temp_comb = temp_comb(1+five_dent : end-three_dent);
        % 2. get pssm score
        allcombs_pssm_score(comb) = ...
        calc_pssm_matching_score( PSSM, temp_comb, 'log' );
    end
    
    % get rid of options above threshold
    del_idx = find(allcombs_pssm_score > threshold);
    
    if length(del_idx) == length(allcombs) 
        %meaning there are no vlaid options, so take the one with lowest
        %pssm:
        % 1. get index
        [~,best_comb] = min(allcombs_pssm_score);
        % 2. and now actual seq
        best_comb = allcombs{best_comb};
    else % if there are some options to pick from
        allcombs(del_idx) = []; 
        % calc similarity to original DNA seq and pick winner:
        % 1. this gives only the index
        best_comb = cellfun(@(x) sum(x~=original_window_DNA), allcombs);
        [~,best_comb] = min(best_comb); 
        % 2. now get the actual sequence
        best_comb = allcombs{best_comb};
    end
    
    % update the sequence (make replacement)
    if sum(nt2aa(seq(beg_idx : end_idx),'AlternativeStartCodons',0) ~= nt2aa(best_comb,'AlternativeStartCodons',0) ) > 0
        error('while fixing a splice site - we were about to create a non-synonymous mutation')
    end
    seq(beg_idx : end_idx) = best_comb;
    
end %of comb for loop

% finalize:
new_seq = seq;

end