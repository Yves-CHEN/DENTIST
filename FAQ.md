## Why error msg "LD was set to 0, because var_i = 0.334057, var_j  0.000000"?
Let LD = cov(i,j)/sqrt( var_i * var_j ). The  var_i or var_j should be not 0.
The above case (var_j) can happen when the MAF of the SNP is 0.


## Why error msg "Divividing zero : Rsq = 1.028902 "
This could occur when genetic data has missing values (NAs).
Try run with --with-NA-geno flag.


