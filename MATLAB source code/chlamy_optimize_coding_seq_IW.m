function [optimized_dna_seq, best_mean_FE] =...
    chlamy_optimize_coding_seq_IW(aa_seq, fourty_nt_upstream_atg, ...
    CDSs, output_file_name, gene_title)
%   This is a master function that organizes the data and sends it
%   to several other optimization functions, each focusing on a 
%   specific aspect of gene expression, and working on the output of the
%   former function 
%   Note: the stop codon is regarded as any other codon
%
%     VARIABLES, INPUT:
%     *aa_seq = amino acid sequence of wanted protein
%     *fourty_nt_upstream_atg = the 40 nt immidiately before first ATG 
%     *CDSs = CDS of reference genes
%
%     VARIABLES, OUTPUT:
%     *dna_seq = final sequence chosen for expression
%     *mean_FE = mean folding energy of mRNA in ATG vicinity
%
%     The function does the following:
%     1) Calculate codon usage bias
%
%     2) Create M variants of the gene 
%
%     3) delete splicing signals
%
%     4) find the sequence with the best folding
%
%     5) Validataion
%
%     6) Create output file

%% 1) get codon bias

[Codon_map, accumulative_CUB] = get_codon_bias(CDSs);
disp('got codon bias')

%% 2) create M variants of the CDS

% find barrier (30 for long genes, 50% for really short genes)
N = length(aa_seq);
if ceil(N * 0.5) < 30
    Barrier = ceil(length(aa_seq) * 0.5);
elseif ceil(N * 0.5) >= 30
    Barrier = 30;
end

% create the genes
M = 100; % decide how many variants to receive in this stage
optimized_CDS = cell(M,1); 

% get chunks
for i = 1:M
    
    first_chunk = copy_codon_distribution(aa_seq(1 : Barrier), accumulative_CUB);
    second_chunk = optimal_codons_distribution(aa_seq(Barrier+1 : end), Codon_map);
    optimized_CDS{i} = [first_chunk second_chunk];
    
    % show progress
    disp(['gene # ' num2str(i) ' was created (optimal codon bias)'])
end

% chuck un-unique variants, if there are any (waste of time down the road)
optimized_CDS = unique(optimized_CDS);

% verify that the AA seq was not altered
if find(cellfun(@(x) sum(nt2aa(x,'AlternativeStartCodons',0)~=aa_seq), optimized_CDS))
    error('non-synonymous sequences were created')
end



%% 3) delete splicing signals in each sequence

splicing_data = 'splicing_signals.mat'; %name of file to load in function
for i = 1:(length(optimized_CDS))
    optimized_CDS{i} = erase_splicing(optimized_CDS{i}, splicing_data);
    disp(['----------seq ' num2str(i) ' cleaned---------------'])
end

% chuck un-unique variants, if there are any (waste of time down the road)
%sound(randn(4096, 1), 8192) %make some sound
optimized_CDS = unique(optimized_CDS);


%% 4) find best folding in ATG area

%progress update
disp('mRNA FE calculations have just begun')

% window size
window_size = 39;

% create sequences of the gene with unique attributes
upstream_chunk = fourty_nt_upstream_atg;
upstream_chunk = upstream_chunk((end-window_size+1) : end);

% here use *Barrier* from section (2)
Sequences_for_FE_check = cellfun(@(x) [upstream_chunk x(1:Barrier*3)],optimized_CDS,...
    'uniformoutput',0);

% find seq with best folding
FE_vals = sliding_window_analysis(Sequences_for_FE_check, window_size, @rnafold, ...
    2, 'figure', 0);

% get mean of each row
FE_vals = mean(FE_vals,2);

% get highest FE (most open)
[best_mean_FE, best_folding_idx] = max(FE_vals);
best_CDS = optimized_CDS{best_folding_idx};

%% 5) validation and context concatination 

if sum(nt2aa(best_CDS,'AlternativeStartCodons',0) ~= aa_seq) == 0
    disp('all went well')
else
    error('The AA seq of the output DNA did not match the input aa seq, look for bug...')
end

optimized_dna_seq = best_CDS;

%% 6) create output file

disp('preparing output file')
gene_length = length(optimized_dna_seq);
gene_GC = GC_content(optimized_dna_seq);

Headers = {gene_title;...
    'Total length';...
    'GC content';...
    'Mean mRNA FE at ATG vicinity';
    '40 NT upstream of gene'};

Content = {optimized_dna_seq;...
    num2str(gene_length);...
    num2str(gene_GC);...
    num2str(best_mean_FE);...
    fourty_nt_upstream_atg};

fastawrite(output_file_name, Headers, Content);

% update
disp(['Run finished for ' gene_title ' --------------------------'])

end