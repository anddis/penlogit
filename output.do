*  Paper: Approximate Bayesian logistic regression
*          via penalized likelihood estimation with data augmentation
*  Authors: A.Discacciati, N.Orsini, S.Greenland


sjlog using penlogit1, replace
clear
input nomonit deaths n
0 3 694     
1 14 2298   
end 
penlogit deaths nomonit, binomial(n) ppl(nomonit) 
sjlog close , replace
  
 
sjlog using penlogit2, replace
penlogit deaths nomonit, binomial(n) ppl(nomonit) ///
 nprior(nomonit ln(2) 0.5) or
sjlog close , replace  


sjlog using penlogit3, replace
penlogit deaths nomonit, binomial(n) ppl(nomonit) ///
 lfprior(nomonit ln(2) 2000 2 1) or
sjlog close , replace  


sjlog using penlogit4, replace
penlogit deaths nomonit, binomial(n) ppl(nomonit) ///
 lfprior(nomonit ln(2) 2000 2 0.5) or
sjlog close , replace 


sjlog using penlogit5, replace
quietly mlexp (ln(invlogit({b0}+{nomonit}*nomonit))*deaths + ///
        ln(1-(invlogit({b0}+{nomonit}*nomonit)))*(n-deaths) + ///
        (1000*(({nomonit}-ln(2))/0.5+ln(2000/2)) - ///
        1001*ln(1+exp(({nomonit}-ln(2))/0.5+ln(2000/2))))/2)
lincom [nomonit]_cons, or
sjlog close , replace


sjlog using penlogit6, replace
use http://www.imm.ki.se/biostatistics/data/neutra1978.dta, clear
penlogit death nomonit teenages gestage abort dyslab ward malpres ///
 nonwhite nullip isoimm hydram placord twint prerupt, ppl(hydram) or
sjlog close , replace 


sjlog using penlogit7, replace
penlogit death nomonit teenages gestage abort dyslab ward malpres ///
 nonwhite nullip isoimm hydram placord twint prerupt, ///
 nprior(nomonit ln(2) 0.5 teenages ln(2) 0.5 gestage ln(4) 0.5 abort 0 0.5 ///
 dyslab ln(2) 0.5 ward ln(2) 0.5 malpres ln(4) 0.5 ///
 nonwhite ln(2) 0.5 nullip ln(2) 0.5 isoimm ln(2) 0.5 ///
 placord ln(2) 0.5 twint ln(4) 0.5 hydram ln(4) 0.5 ///
 prerupt ln(2) 0.5) ///
 ppl(nomonit teenages gestage abort dyslab ward malpres ///
 nonwhite nullip isoimm hydram placord twint prerupt) ///
 nppl(50) or
sjlog close , replace 


sjlog using penlogit8, replace
clear
set obs 500
local intercept = -11.6
local beta = ln(4)
local xbeta ""
foreach i of numlist 1/10 {
	generate x`i' = rbinomial(1, 0.5)
	local xbeta "`xbeta' + `beta' * x`i'"
}
generate xb = `intercept' `xbeta'
generate y = rbinomial(1, invlogit(xb))
sjlog close, replace

