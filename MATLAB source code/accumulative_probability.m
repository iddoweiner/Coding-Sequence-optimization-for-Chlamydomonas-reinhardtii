function Bias_for_probability_algorithm=accumulative_probability(seq)

%     this function receives a sequence as a string. a common example would be 
%     a long string of concatenated genes for which we want to calculate 
%     codon bias. Instead of giving the regular probability for each codon
%     within an amino acid, so that the sum of probabilities equals 1, this 
%     func' gives accumulative probabilities, so that the last codon for each aa
%     equals 1 itself. This is good for algorithms in which you want to copy the
%     codon frequency of genes in a new sequence

%% create special bias - accumulative probability
Bias=codonbias(seq);
aa_names=fieldnames(Bias);
for i=1:(length(aa_names))
   for j=2:(length(Bias.(aa_names{i}).Freq))
       Bias.(aa_names{i}).Freq(j)=Bias.(aa_names{i}).Freq(j-1)+...
           Bias.(aa_names{i}).Freq(j);
   end
end

Bias_for_probability_algorithm=Bias;

end