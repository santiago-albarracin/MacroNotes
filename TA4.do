*-------------------------------------------------------------------------------
* Estimation of income process - Valerio Pieroni
*-------------------------------------------------------------------------------

clear 
set more off
cls

* Set up the data --------------------------------------------------------------

* load data
use panel_reg, replace

* describe 
des 

* declare panel 
xtset persnr year

* restrict household head age
drop if age_head < 25 | age_head > 60
* real earnings using consumer price index 
gen y5 = yhh5/cpi
gen y6 = yhh6/cpi
* restrict earnings > 0
gen ln_inc5 = log(y5)
gen ln_inc6 = log(y6)

* controls 
gen age2 = age_head^2
gen age3 = age_head^3
* education dummies 
gen coll = 0
replace coll = 1 if edu == 4 
gen hs = 0
replace hs = 1 if  edu == 2 | edu == 3
* interactions 
gen agecoll = age_head*coll
gen agecoll2 = age2*coll
gen agecoll3 = age3*coll
gen agehs = age_head*hs
gen agehs2 = age2*hs
gen agehs3 = age3*hs

* first stage estimates ----------------------------------------------------------------------------

xtreg ln_inc5 age_head age2 age3 coll hs famsize married /*

               */ agecoll agecoll2 agecoll3 agehs agehs2 agehs3 i.year, fe cluster(persnr)

estimates store model1

xtreg ln_inc6 age_head age2 age3 coll hs famsize married /*

               */ agecoll agecoll2 agecoll3 agehs agehs2 agehs3 i.year, fe cluster(persnr)
			   
estimates store model2

esttab model1 model2 using table1.tex, cells(b(star fmt(3)) se(par fmt(3)))


* predicted income
predict hat_ln_inc5
predict hat_ln_inc6
scatter hat_ln_inc5 age_head

* smoothed age earnings profile
reg hat_ln_inc5 age_head age2 age3 married
predict sm_inc
scatter sm_inc age_head
* correlation between hours and age
reg Hrs_hh age_head
predict hat_hours
scatter hat_hours age_head

* residual income
gen res_ln_inc5 = ln_inc5 - hat_ln_inc5
gen res_ln_inc6 = ln_inc6 - hat_ln_inc6
summarize res_ln_inc5 res_ln_inc6, d

* save data 
save panel_sel, replace
		   
* second stage estimates -----------------------------------------------------------------

clear 
cls

* GMM estimation of AR(1) model
do "gmm_reg.do" 










