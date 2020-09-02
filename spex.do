clear all
set more off
// adopath ++ ~/Dropbox/Cattaneo-Drukker-Holland/code
// Set adopath outside of script
set linesize 80
mata: mata mlib index
set scheme sj

sjlog using ex0, replace
use spmdata
bfit logit w pindex eindex, corder(3) base(0) sort(aic)
sjlog close, replace

sjlog using ex1, replace
mlogit
sjlog close, replace

sjlog using ex2, replace
predict double (phat0 phat1 phat2), pr
sort w
by w: summarize phat0 phat1 phat2
sjlog close, replace

// overlap plot
sjlog using ex3, replace
kdensity phat0 if w==0, generate(xp00 den00) nograph n(5000) kernel(triangle)
kdensity phat0 if w==1, generate(xp01 den01) nograph n(5000) kernel(triangle)
kdensity phat0 if w==2, generate(xp02 den02) nograph n(5000) kernel(triangle)
twoway line den00 xp00 || line den01 xp01 || line den02 xp02 , 	///
	legend(label(1 "w==0") label(2 "w==1") label(3 "w==2"))	///
	title("Conditional densities for probability of treatment level 0") ///
	name(pw0)
sjlog close, replace
graph export pw0.eps, replace

sjlog using ex4, replace
kdensity phat1 if w==0, generate(xp10 den10) nograph n(5000) kernel(triangle)
kdensity phat1 if w==1, generate(xp11 den11) nograph n(5000) kernel(triangle)
kdensity phat1 if w==2, generate(xp12 den12) nograph n(5000) kernel(triangle)

kdensity phat2 if w==0, generate(xp20 den20) nograph n(5000) kernel(triangle)
kdensity phat2 if w==1, generate(xp21 den21) nograph n(5000) kernel(triangle)
kdensity phat2 if w==2, generate(xp22 den22) nograph n(5000) kernel(triangle)
sjlog close, replace

sjlog using ex5, replace
twoway line den10 xp10 || line den11 xp11 || line den12 xp12 , 	///
	legend(label(1 "w==0") label(2 "w==1") label(3 "w==2"))	///
	title("Conditional densities for probability of treatment level 1") ///
	name(pw1)
twoway line den20 xp20 || line den21 xp21 || line den22 xp22 , 	///
	legend(label(1 "w==0") label(2 "w==1") label(3 "w==2"))	///
	title("Conditional densities for probability of treatment level 2") ///
	name(pw2)
sjlog close, replace
graph export pw1.eps , name(pw1) replace
graph export pw2.eps , name(pw2) replace

sjlog using ex6, replace
bfit regress spmeasure pindex eindex, corder(3) 
sjlog close, replace

sjlog using ex7, replace
regress
sjlog close, replace

sjlog using ex7b, replace
quietly bfit regress spmeasure pindex eindex, corder(3) 
regress spmeasure `r(bvlist)'
sjlog close, replace

sjlog using ex8, replace
poparms (w c.(pindex eindex)##c.(pindex eindex)) 	///
	(spmeasure c.(pindex eindex)##c.(pindex eindex)) 	
sjlog close, replace

sjlog using ex9, replace
poparms  , coeflegend
contrast ar.w, nowald
sjlog close, replace

sjlog using ex10, replace
contrast r.w, nowald
sjlog close, replace

sjlog using ex11, replace
margins i.w, pwcompare
marginsplot, unique plotopts(connect(none))
sjlog close, replace
graph export mpw.eps , replace

// use spmdata
sjlog using ex12, replace
set seed 12345671
poparms (w c.(pindex eindex)##c.(pindex eindex)) 	        ///
	(spmeasure c.(pindex eindex)##c.(pindex eindex)) ,      ///
	quantile(.25 .5 .75) 
sjlog close, replace

sjlog using ex13, replace
margins i.w , pwcompare predict(equation(#2))
sjlog close, replace

sjlog using ex13a, replace
margins i.w , pwcompare predict(equation(#3))
sjlog close, replace

sjlog using ex14, replace
margins i.w , pwcompare predict(equation(#4))
sjlog close, replace


// New stuff Joint inference

sjlog using ex15, replace
test (_b[mean:0.w] = _b[mean:1.w] = _b[mean:2.w])	///
     (_b[q25:0.w] = _b[q25:1.w] = _b[q25:2.w])		///
     (_b[q50:0.w] = _b[q50:1.w] = _b[q50:2.w])		///
     (_b[q75:0.w] = _b[q75:1.w] = _b[q75:2.w])
sjlog close, replace


sjlog using ex16, replace
test (_b[mean:1.w] - _b[mean:0.w] = _b[mean:2.w] - _b[mean:1.w] )	///
     (_b[q25:1.w] - _b[q25:0.w] = _b[q25:2.w] - _b[q25:1.w] )		///
     (_b[q50:1.w] - _b[q50:0.w] = _b[q50:2.w] - _b[q50:1.w] )		///
     (_b[q75:1.w] - _b[q75:0.w] = _b[q75:2.w] - _b[q75:1.w] )	
sjlog close, replace

sjlog using ex16a, replace
test _b[mean:1.w] - _b[mean:0.w] = _b[mean:2.w] - _b[mean:1.w] 
sjlog close, replace

sjlog using ex16b, replace
test _b[q50:1.w] - _b[q50:0.w] = _b[q50:2.w] - _b[q50:1.w] 	
sjlog close, replace
     

sjlog using ex17, replace
test (_b[mean:0.w] = _b[q50:0.w])	///
     (_b[mean:1.w] = _b[q50:1.w])	///
     (_b[mean:2.w] = _b[q50:2.w])
sjlog close, replace

generate ivalues  = .
label variable ivalues "Treatment level"

local eqn mean q25 q50 q75
local lpl solid dash dot dash_dot
local eqi 1
foreach eq of local eqn {
	generate bvalues_`eq'  = .
	generate ci_low_`eq'   = .
	generate ci_high_`eq'  = .

	local lpat : word `eqi' of `lpl'
	forvalues i=0/2 {
		local ip1 = `i' + 1
		if "`eq'" == "mean" {
			replace ivalues = `i' in `ip1'
		}

		replace bvalues_`eq' = _b[`eq':`i'.w] in `ip1'
		replace ci_low_`eq'  = bvalues_`eq' - invnormal(.975)*_se[`eq':`i'.w] in `ip1'
		replace ci_high_`eq' = bvalues_`eq' + invnormal(.975)*_se[`eq':`i'.w] in `ip1'
	}
	local eq1 " (connected bvalues_`eq' ivalues in 1/3, lpattern(`lpat') msize(*.5) ) "
	local eq2 " (rcap ci_low_`eq' ci_high_`eq' ivalues in 1/3, lcolor(black)) "
	local gcmd "`gcmd' `eq1' `eq2' "
	local ++eqi
}
display `"graph twoway `gcmd' "'
graph twoway `gcmd' , legend(order(1 3 5 7 8))	///
	legend(label( 1 "Mean"))		///
	legend(label( 3 "25th quantile"))	///
	legend(label( 5 "Median"))		///
	legend(label( 7 "75th quantile"))	///
	legend(label( 8 "95% confidence intervals"))	///
	xlabel(0(1)2) ytitle("Parameter estimate")

graph export effects.eps, replace

