% function for calculating GC content of a sequence
function [y]=GC_content(sequence)
if isnan(sequence)
    y = NaN;
else
y = (length(regexpi(sequence,'g')) + length(regexpi(sequence,'c'))) / ...
    (length(regexpi(sequence,'g')) + length(regexpi(sequence,'c')) + length(regexpi(sequence,'a')) + length(regexpi(sequence,'t')) );
end
end
