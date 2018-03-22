function logical = is_valid_AA_seq(Seq);
logical = 1;
aa_seq = 'ARNDCEQGHILKMFPSTWYV*';
allowed_bases = cellstr(aa_seq')';
for i = 1:(length(Seq))
    if ismember(Seq(i),allowed_bases)
        continue
    else
        logical = 0;
        break
    end
end
end