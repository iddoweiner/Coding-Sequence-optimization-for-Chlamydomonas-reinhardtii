% run_optimization_script
% this is an example for runnig the sequence optimization script

% NOTE: make sure all the functions and data files are in your working path

% 1) load the CDS data
% (the example here is with the real data, no modification required)
load('ref_CDSs.mat')

% 2) give your amino acid sequence
% (modify the example)
aa_seq = 'M*';

% 3) give the 40 NT which will appear upstream from your gene 
% (modify the example)
fourty_nt_upstream_atg = randseq(40);

% 4) give the title for your gene and the file name
% (modify the example)
output_file_name = 'optimized_seq.fa';
gene_title = 'best_chlamy_gene';

% run the function
[optimized_dna_seq, best_mean_FE] = chlamy_optimize_coding_seq_IW(aa_seq, ...
    fourty_nt_upstream_atg, ref_CDSs, output_file_name, gene_title);