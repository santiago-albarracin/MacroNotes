*------------------------------------------------------------------------------
* Project 2 - Valerio Pieroni
*------------------------------------------------------------------------------

clear 
set more off
cls

* load data
use cexdata, replace
* describe 
des 

* full year 
forvalues i = 6/9 {
replace year = 198`i' if year == 8`i'
}
forvalues i = 0/9 {
replace year = 199`i' if year == 9`i'
}
replace year = 2000 if year == 00

* sort by year quarter id
sort year quarter id

* quarterly average consumption 
by year quarter (id) : egen avgcons = mean(ndcons1)
by year quarter (id) : egen avgcons_f = mean(food)
by year quarter (id) : egen avgcons_t = mean(totcons)
* quarterly average income
by year quarter (id) : egen avginc = mean(grinc)

* date
gen date = yq(year, quarter)
format date %tqCCYY/!Qq

* descriptive stats ------------------------------------------------------------

* consumption and income stats
forvalues i = 1/4 {
sum netinc totcons ndcons1 if quarter == `i' & year == 1986
sum netinc totcons ndcons1 if quarter == `i' & year == 2000
}

* business cycles --------------------------------------------------------------

* plot the time series 
tsset id date 
tsline avgcons 
tsline avginc

* complete markets test --------------------------------------------------------

* consider only second and fifth interviews in eahc year
drop if inumb == 3 | inumb == 4

* log variables 
gen licons = log(ndcons1)
gen licons_f = log(food)
gen licons_t = log(totcons)
gen lcons = log(avgcons)
gen lcons_f = log(avgcons_f)
gen lcons_t = log(avgcons_t)
gen linc = log(netinc)

* first difference for each household and year between 2nd and 5th interview (3 quartes)
sort year id inumb
by year id (inumb) : gen dlc = lcons - lcons[_n-1]
by year id (inumb) : gen dlc_f = lcons_f - lcons_f[_n-1]
by year id (inumb) : gen dlc_t = lcons_t - lcons_t[_n-1]
by year id (inumb) : gen dlci = licons - licons[_n-1]
by year id (inumb) : gen dlci_f = licons_f - licons_f[_n-1]
by year id (inumb) : gen dlci_t = licons_t - licons_t[_n-1]
by year id (inumb) : gen dlyi = linc - linc[_n-1]

sum dlci, d

* food
regress dlci_f dlc_f dlyi, cluster(id) 
test (dlc_f = 1) (dlyi = 0) 
* nondurables
regress dlci dlc dlyi, cluster(id) 
test (dlc = 1) (dlyi = 0) 
* total consumption
regress dlci_t dlc_t dlyi, cluster(id) 
test (dlc_t = 1) (dlyi = 0) 



