function aa3 = aa1_convert_aa3(aa1)

   % get data
load('AA_name_conversion.mat')  

if iscell(aa1)
   L = length(aa1);
   aa3 = cell(L,1);
   for i = 1:L
      aa3{i} = one_letter_key(aa1{i});
   end

elseif isstr(aa1) %this means there's only 1 AA
    aa3 = one_letter_key(aa1);
    
end %of if

end