function optimal_codons_DNA_seq = optimal_codons_distribution(aa_seq, Codon_map)
% optimal_codons_DNA_seq = optimal_codons_distribution(aa_seq, normal_CUB)
% input: 
% aa_seq = string of aa being optimized 
% normal_CUB = codon usage struct 
% output: 
% optimal_codons_DNA_seq = a DNA seq where only the optimal
%   codons are used, it codes of course for the AA seq given

%% go over the sequence and asign optimal codons
% open empty string
optimal_codons_DNA_seq = '';
% go over seq
for i = 1:(length(aa_seq))
    optimal_codons_DNA_seq = [optimal_codons_DNA_seq, Codon_map(aa_seq(i))];
end

end