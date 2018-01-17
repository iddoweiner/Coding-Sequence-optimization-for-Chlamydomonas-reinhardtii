function [Codon_map, accumulative_CUB, normal_CUB] = get_codon_bias(CDSs)

aa_seq = 'ARNDCEQGHILKMFPSTWYV*';

% get gene string
gene_string = strjoin(CDSs,'');
% create special bias - accumulative probability
accumulative_CUB = accumulative_probability(gene_string);
% create map for connecting optimal codon with 1 AA letter
normal_CUB = codonbias(gene_string);
% create a map of keys(1 letter AA name) and values (right codon) for the relevant codons (only those in the sequence)
all_AAs = cellstr(aa_seq')';
% get only unique entries
uniq_AA_1letter = unique(all_AAs)';
% turn into 3 letter AA
uniq_AA_3letter = cellfun(@aa1_convert_aa3, uniq_AA_1letter, 'un', 0);
% create list for optimal codons
Optimal_codon = cell(size(uniq_AA_3letter));
% fill the list
for i = 1:(length(Optimal_codon))
    temp_freq = normal_CUB.(uniq_AA_3letter{i}).Freq;
    [~,IDX] = max(temp_freq);
    Optimal_codon{i} = normal_CUB.(uniq_AA_3letter{i}).Codon{IDX};
end
% create the map
Codon_map = containers.Map(uniq_AA_1letter, Optimal_codon);

end