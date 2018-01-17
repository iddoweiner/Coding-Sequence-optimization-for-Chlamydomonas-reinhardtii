function codon_distrubution_optimized_gene = copy_codon_distribution(aa_seq,codon_bias)

%     this function copies the codon bias of a reference set. It is originally
%     used as the first part of optimizing genes for expression in chlamy

%     VARIABLES:
%     aa_seq = amino acid sequence of wanted protein, including "stop"
%     1) codon_bias = struct with special bias (output of the function 
%        accumulative_probability*) in the 100 most transcribed genes 
%     2) codon_distrubution_optimized_gene = DNA string with codons that match 
%        the bias of the reference set

N = length(aa_seq);
probability_vec = rand(1,N);

aa_seq_cell=cell(1,N);
for i=1:N
    aa_seq_cell{i}=aa_seq(i);
end

% convert each letter to a triplet
aa_seq_cell=aa1_convert_aa3(aa_seq_cell);

% create cell for matching codons
codon_seq_cell=cell(1,N); 

% fill cell with the right codons
for i=1:N
    % create easy handle names
    temp_aa_name=aa_seq_cell{i}; %name of current aa
    temp_probability=probability_vec(i); %current probability number
    accumulative_probabilities=codon_bias.(temp_aa_name).Freq; %get all frequencies for current aa
    
    % find the first number that is smaller than 
    % temp_probability (called  idx_first_bigger_number)
    go=1;
    j=1;
    while go
        if temp_probability<accumulative_probabilities(j)
            idx_first_bigger_number=j;
            go=0;
        else
            j=j+1;
        end    
    end %while
    
    % use the index to pick the right codon and put it in the cell
    codon_seq_cell{i}=codon_bias.(temp_aa_name).Codon{idx_first_bigger_number};
end

% concatenate strings from all cells to form a gene with codon
% distribution similar to the mean distribution in 100 highly
% transcribed chlamy genes
codon_distrubution_optimized_gene='';
for i=1:N
    codon_distrubution_optimized_gene=strcat(codon_distrubution_optimized_gene,...
        codon_seq_cell{i});
end

end