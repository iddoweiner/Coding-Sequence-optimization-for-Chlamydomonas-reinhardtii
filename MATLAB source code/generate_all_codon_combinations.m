function all_codon_combinations = generate_all_codon_combinations(aa_seq)
%     all_codon_combinations = generate_all_codon_combinations(aa_seq)
%     get all codon combinations for an AA sequence / Iddo Weiner
%     
%     IN:
%     1 letter amino acid sequence (string)
%     
%     OUT:
%     a cellarray containing all synonymous DNA sequences (i.e. codon combinations)
%     that code the specified amino acid sequence
%
%     EXAMPLE:
%     generate_all_codon_combinations('MCHQ*')
%     ans = 
%         'ATGTGCCACCAATAA'
%         'ATGTGCCACCAATAG'
%         'ATGTGCCACCAATGA'
%         'ATGTGCCACCAGTAA'
%         'ATGTGCCACCAGTAG'
%         'ATGTGCCACCAGTGA'
%         'ATGTGCCATCAATAA'
%         'ATGTGCCATCAATAG'
%         'ATGTGCCATCAATGA'
%         'ATGTGCCATCAGTAA'
%         'ATGTGCCATCAGTAG'
%         'ATGTGCCATCAGTGA'
%         'ATGTGTCACCAATAA'
%         'ATGTGTCACCAATAG'
%         'ATGTGTCACCAATGA'
%         'ATGTGTCACCAGTAA'
%         'ATGTGTCACCAGTAG'
%         'ATGTGTCACCAGTGA'
%         'ATGTGTCATCAATAA'
%         'ATGTGTCATCAATAG'
%         'ATGTGTCATCAATGA'
%         'ATGTGTCATCAGTAA'
%         'ATGTGTCATCAGTAG'
%         'ATGTGTCATCAGTGA'

%% verify/fix input
% fix case
aa_seq = upper(aa_seq);
% load codon map (.mat file)
load('AA_name_conversion.mat');
% verify valid AA chars
if length(find(isKey(one_letter_key_2_all_codons, cellstr(aa_seq')) == 0)) > 0
    error('at least one of the chars in the input sequence is not a valid 1 letter amino acid')
end

%% get all codons for each amino acid      
% get length of input sequence
N = length(aa_seq);
% list relevant options
options = cell(1,N);
for i = 1:N
    options{i} = one_letter_key_2_all_codons(aa_seq(i));
end

%% send options to the combinatorial function
all_codon_combinations = get_all_combinations(options);

end