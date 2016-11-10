clear all

*------- PROGRAMS USED IN SECTION 4 -------*

*--- mcsim ---*
capture program drop mcsim
program define mcsim, rclass
syntax [anything] [, xn(string) ///
beta(string) intercept(string) p(string) NOBS(integer 100)]
clear

set obs `nobs'

foreach g of numlist 1/`xn' {
	gen x`g' = runiform() < `p'
}

local model ""
foreach i of numlist 1/`xn' {
	local model "`model' + `beta' * x`i'"
}

gen logitp = `intercept' `model'

gen py = invlogit(logitp)
gen y = runiform() < py 

tab y

drop py logitp
contract x* y, f(fw)

capture noi logit y x* [fw = fw], or noheader nolog asis
mat coef = e(b)
return scalar beta1 = coef[1,1]
return scalar beta2 = coef[1,2]
return scalar beta3 = coef[1,3]
return scalar beta4 = coef[1,4]
return scalar beta5 = coef[1,5]
return scalar beta6 = coef[1,6]
return scalar beta7 = coef[1,7]
return scalar beta8 = coef[1,8]
return scalar beta9 = coef[1,9]
return scalar beta10 = coef[1,10]
return scalar converged = e(converged)

penlogit y x* [fw = fw], nprior(x1 0 4 x2 0 4 x3 0 4 ///
x4 0 4 x5 0 4 x6 0 4 ///
x7 0 4 x8 0 4 x9 0 4 x10 0 4)  or nolist
mat coefp = e(b)
return scalar betapn1 = coefp[1,1]
return scalar betapn2 = coefp[1,2]
return scalar betapn3 = coefp[1,3]
return scalar betapn4 = coefp[1,4]
return scalar betapn5 = coefp[1,5]
return scalar betapn6 = coefp[1,6]
return scalar betapn7 = coefp[1,7]
return scalar betapn8 = coefp[1,8]
return scalar betapn9 = coefp[1,9]
return scalar betapn10 = coefp[1,10]
return scalar convergedpn = e(converged)


penlogit y x* [fw = fw], nprior(x1 ln(2) 1 x2 ln(2) 1 x3 ln(2) 1 ///
x4 ln(2) 1 x5 ln(2) 1 x6 ln(2) 1 ///
x7 ln(2) 1 x8 ln(2) 1 x9 ln(2) 1 x10 ln(2) 1)  or nolist
mat coefp = e(b)
return scalar betapi1 = coefp[1,1]
return scalar betapi2 = coefp[1,2]
return scalar betapi3 = coefp[1,3]
return scalar betapi4 = coefp[1,4]
return scalar betapi5 = coefp[1,5]
return scalar betapi6 = coefp[1,6]
return scalar betapi7 = coefp[1,7]
return scalar betapi8 = coefp[1,8]
return scalar betapi9 = coefp[1,9]
return scalar betapi10 = coefp[1,10] 
return scalar convergedpi = e(converged)

end
*------*

*--- report ---*
capture program drop report
program define report
foreach v of varlist beta* {
	gen or`v' = exp(`v')
}
tab1 converged convergedpn convergedpi
tabstat orbeta1-orbeta10 if converged  , s(p50 p5 p95 n) format(%4.1f)
tabstat orbetapn1-orbetapn10 if convergedpn, s(p50 p5 p95 n) format(%4.1f)
tabstat orbetapi1-orbetapi10 if convergedpi, s(p50 p5 p95 n) format(%4.1f)
end
*------*

*--- intercept_bisect ---*
capture program drop intercept_bisect
program intercept_bisect, rclass
	args or p
	local x1 = -25
	local x2 = +10
while abs(`x1' - `x2') > 1e-3 {
	local t = (`x1' + `x2') / 2
	local a = binomialp(10,0,.5)*invlogit(`t') + ///
				binomialp(10,1,.5)*invlogit((`t'+1*log(`or'))) + ///
				binomialp(10,2,.5)*invlogit((`t'+2*log(`or'))) + ///
				binomialp(10,3,.5)*invlogit((`t'+3*log(`or'))) + ///
				binomialp(10,4,.5)*invlogit((`t'+4*log(`or'))) + ///
				binomialp(10,5,.5)*invlogit((`t'+5*log(`or'))) + ///
				binomialp(10,6,.5)*invlogit((`t'+6*log(`or'))) + ///
				binomialp(10,7,.5)*invlogit((`t'+7*log(`or'))) + ///
				binomialp(10,8,.5)*invlogit((`t'+8*log(`or'))) + ///
				binomialp(10,9,.5)*invlogit((`t'+9*log(`or'))) + ///
				binomialp(10,10,.5)*invlogit((`t'+10*log(`or')))
	local diff = `p' - `a'
	if `diff' > 0 local x1 = `t'
	else local x2 = `t'
}
return scalar intercept = round(`t',.01)
end
*------*

exit




*------- SIMULATIONS -------*

cls

log using simulation, smcl replace name(sim)

*n=500, beta=log(4)
intercept_bisect 4 .05
local intercept = r(intercept)

simulate ///
beta1 = r(beta1) ///
beta2 = r(beta2) ///
beta3 = r(beta3) ///
beta4 = r(beta4) ///
beta5 = r(beta5) ///
beta6 = r(beta6) ///
beta7 = r(beta7) ///
beta8 = r(beta8) ///
beta9 = r(beta9) ///
beta10 = r(beta10) ///
converged = r(converged) ///
betapn1 = r(betapn1) ///
betapn2 = r(betapn2) ///
betapn3 = r(betapn3) ///
betapn4 = r(betapn4) ///
betapn5 = r(betapn5) ///
betapn6 = r(betapn6) ///
betapn7 = r(betapn7) ///
betapn8 = r(betapn8) ///
betapn9 = r(betapn9) ///
betapn10 = r(betapn10) ///
convergedpn = r(convergedpn) ///
betapi1 = r(betapi1) ///
betapi2 = r(betapi2) ///
betapi3 = r(betapi3) ///
betapi4 = r(betapi4) ///
betapi5 = r(betapi5) ///
betapi6 = r(betapi6) ///
betapi7 = r(betapi7) ///
betapi8 = r(betapi8) ///
betapi9 = r(betapi9) ///
betapi10 = r(betapi10) ///
convergedpi = r(convergedpi) ///
, ///
reps(1000)  seed(4321) saving(s1, replace) nodots: mcsim , xn(10) beta(log(4)) p(.5) ///
intercept(`intercept') nobs(500)

report


*n=500, beta=log(10)
intercept_bisect 10 .05
local intercept = r(intercept)

simulate ///
beta1 = r(beta1) ///
beta2 = r(beta2) ///
beta3 = r(beta3) ///
beta4 = r(beta4) ///
beta5 = r(beta5) ///
beta6 = r(beta6) ///
beta7 = r(beta7) ///
beta8 = r(beta8) ///
beta9 = r(beta9) ///
beta10 = r(beta10) ///
converged = r(converged) ///
betapn1 = r(betapn1) ///
betapn2 = r(betapn2) ///
betapn3 = r(betapn3) ///
betapn4 = r(betapn4) ///
betapn5 = r(betapn5) ///
betapn6 = r(betapn6) ///
betapn7 = r(betapn7) ///
betapn8 = r(betapn8) ///
betapn9 = r(betapn9) ///
betapn10 = r(betapn10) ///
convergedpn = r(convergedpn) ///
betapi1 = r(betapi1) ///
betapi2 = r(betapi2) ///
betapi3 = r(betapi3) ///
betapi4 = r(betapi4) ///
betapi5 = r(betapi5) ///
betapi6 = r(betapi6) ///
betapi7 = r(betapi7) ///
betapi8 = r(betapi8) ///
betapi9 = r(betapi9) ///
betapi10 = r(betapi10) ///
convergedpi = r(convergedpi) ///
, ///
reps(1000)  seed(4321) saving(s2, replace) nodots: mcsim , xn(10) beta(log(10)) p(.5) ///
intercept(`intercept') nobs(500)

report


*n=5000, beta=log(4)
intercept_bisect 4 .005
local intercept = r(intercept)

simulate ///
beta1 = r(beta1) ///
beta2 = r(beta2) ///
beta3 = r(beta3) ///
beta4 = r(beta4) ///
beta5 = r(beta5) ///
beta6 = r(beta6) ///
beta7 = r(beta7) ///
beta8 = r(beta8) ///
beta9 = r(beta9) ///
beta10 = r(beta10) ///
converged = r(converged) ///
betapn1 = r(betapn1) ///
betapn2 = r(betapn2) ///
betapn3 = r(betapn3) ///
betapn4 = r(betapn4) ///
betapn5 = r(betapn5) ///
betapn6 = r(betapn6) ///
betapn7 = r(betapn7) ///
betapn8 = r(betapn8) ///
betapn9 = r(betapn9) ///
betapn10 = r(betapn10) ///
convergedpn = r(convergedpn) ///
betapi1 = r(betapi1) ///
betapi2 = r(betapi2) ///
betapi3 = r(betapi3) ///
betapi4 = r(betapi4) ///
betapi5 = r(betapi5) ///
betapi6 = r(betapi6) ///
betapi7 = r(betapi7) ///
betapi8 = r(betapi8) ///
betapi9 = r(betapi9) ///
betapi10 = r(betapi10) ///
convergedpi = r(convergedpi) ///
, ///
reps(1000)  seed(4321) saving(s3, replace) nodots: mcsim , xn(10) beta(log(4)) p(.5) ///
intercept(`intercept') nobs(5000)

report


*n=5000, beta=log(10)
intercept_bisect 10 .005
local intercept = r(intercept)

simulate ///
beta1 = r(beta1) ///
beta2 = r(beta2) ///
beta3 = r(beta3) ///
beta4 = r(beta4) ///
beta5 = r(beta5) ///
beta6 = r(beta6) ///
beta7 = r(beta7) ///
beta8 = r(beta8) ///
beta9 = r(beta9) ///
beta10 = r(beta10) ///
converged = r(converged) ///
betapn1 = r(betapn1) ///
betapn2 = r(betapn2) ///
betapn3 = r(betapn3) ///
betapn4 = r(betapn4) ///
betapn5 = r(betapn5) ///
betapn6 = r(betapn6) ///
betapn7 = r(betapn7) ///
betapn8 = r(betapn8) ///
betapn9 = r(betapn9) ///
betapn10 = r(betapn10) ///
convergedpn = r(convergedpn) ///
betapi1 = r(betapi1) ///
betapi2 = r(betapi2) ///
betapi3 = r(betapi3) ///
betapi4 = r(betapi4) ///
betapi5 = r(betapi5) ///
betapi6 = r(betapi6) ///
betapi7 = r(betapi7) ///
betapi8 = r(betapi8) ///
betapi9 = r(betapi9) ///
betapi10 = r(betapi10) ///
convergedpi = r(convergedpi) ///
, ///
reps(1000)  seed(4321) saving(s4, replace) nodots: mcsim , xn(10) beta(log(10)) p(.5) ///
intercept(`intercept') nobs(5000)

report

log close sim
