function logical = is_valid_DNA_seq(Seq);
logical = 1;
allowed_bases = {'A','a','C','c','G','g','T','t'};
for i = 1:(length(Seq))
    if ismember(Seq(i),allowed_bases)
        continue
    else
        logical = 0;
        break
    end
end
end