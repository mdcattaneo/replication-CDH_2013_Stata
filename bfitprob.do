
capture log close
log using bfitp, text replace
clear all
adopath ++.

query born 
query compilenumber 
which bfit

sysuse auto, clear
generate binary=round(runiform())   // You can generate more categorical
				    //  variables if you want, but the
				    //  result is the same
bfit logit foreign binary rep78

capture log close
