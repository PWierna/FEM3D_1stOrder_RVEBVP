function isequal = checkAisequal2B( refmat , checkmat )

%Quick script to compare matrices with dim up to 3%

%Hard setting: relative tolerance
reltol = 1e-10;

%Get max abs value of ref mat:
maxrefval = max(max(max( abs(refmat) )));
if maxrefval <= reltol %avoid NaN
    maxrefval = 1; 
end

%Compute element-wise relative difference:
reldiff = abs(refmat-checkmat) / maxrefval;

%Check:
isequal = all(all(all( reldiff<=reltol )));

%If doesnt pass:
if isequal == 0
    maxerr = max(max(max(reldiff)));
    disp(maxerr);
end

end