# delimit ;


*-------------------------------------------------------------------------------------------;
* get time temperatures are within each 1C interval for exposure histogram;
use fips year month time*
	if (month >= 4) & (month <= 10) & (year >= 1950) & (year <= 2020) 
	using ../dataSTATA/weather_cotton, clear;
collapse (sum) time*, by(fips year);
collapse (mean) time*;
keep time0C-time40C;
* note: category above 39 is open-ended, so set time40C to zero;
qui replace time40C = 0;
* save as matrix;
mkmat time*, matrix(exposure);

*-------------------------------------------------------------------------------------------;
* initialize x-axis: 161 points in 0.25 temperature steps between 0 and 40;
clear;
local Nobs = 161;
clear;
set obs `Nobs';
gen temperature = (_n-1)/4;


*-------------------------------------------------------------------------------------------;
* load regression results - Chebyshev polynomials;
estimates use ../resultSTATA/panelRegression/cotton_chebyTemperature6_39_year1950_2020_month4_10_clusterState;

* transform temperature to [-1,1]- note: Chebyshev polynomials was between 6 and 39C;
gen tCheby = (temperature-6)/33*2 - 1 if (temperature >= 6) & (temperature <= 39);
* define 8-th order Chebyshev polynomials recursively;
gen cheby1 = tCheby;
gen cheby2 = 2*tCheby*cheby1 - 1;
forvalues p = 3/8 {;
	gen cheby`p' = 2*tCheby*cheby`=`p'-1' - cheby`=`p'-2';
};

* get the estimation results and confidence interval;
gen b = .;
gen lb = .;
gen ub = .;	
local i = 6*4+1;
qui lincom    `=cheby1[`i']' * _b[cheby1] + `=cheby2[`i']' * _b[cheby2] + `=cheby3[`i']' * _b[cheby3] + `=cheby4[`i']' * _b[cheby4]
			+ `=cheby5[`i']' * _b[cheby5] + `=cheby6[`i']' * _b[cheby6] + `=cheby7[`i']' * _b[cheby7] + `=cheby8[`i']' * _b[cheby8]; 
qui replace b = r(estimate) in `i';
qui replace lb = r(estimate)-1.96*r(se) in `i';
qui replace ub = r(estimate)+1.96*r(se) in `i';
* add constant part below lower bound of 6C - results were truncated there as very little temperature mass below that threshold;
qui replace b = r(estimate) if (temperature < 6);
qui replace lb = r(estimate)-1.96*r(se) if (temperature < 6);
qui replace ub = r(estimate)+1.96*r(se) if (temperature < 6);

* range 6-39C is where Chebyshev polynomial is defined;
forvalues i = `=6*4+2'/`=39*4+1' {;
	qui lincom    `=cheby1[`i']' * _b[cheby1] + `=cheby2[`i']' * _b[cheby2] + `=cheby3[`i']' * _b[cheby3] + `=cheby4[`i']' * _b[cheby4]
				+ `=cheby5[`i']' * _b[cheby5] + `=cheby6[`i']' * _b[cheby6] + `=cheby7[`i']' * _b[cheby7] + `=cheby8[`i']' * _b[cheby8]; 
	qui replace b = r(estimate) in `i';
	qui replace lb = r(estimate)-1.96*r(se) in `i';
	qui replace ub = r(estimate)+1.96*r(se) in `i';
};
* add constant part above upper bound of 39C - results were truncated there as very little temperature mass above that threshold;
qui replace b = r(estimate) if (temperature > 39);
qui replace lb = r(estimate)-1.96*r(se) if (temperature > 39);
qui replace ub = r(estimate)+1.96*r(se) if (temperature > 39);


*-------------------------------------------------------------------------------------------;
* load regression results - piecewise linear;
estimates use ../resultSTATA/panelRegression/cotton_piecewiseLinear15_31_39_year1950_2020_month4_10_clusterState;
gen bLin = 0 if (temperature <= 15);
qui replace  bLin = _b[ddayMod]*(temperature-15) if (temperature > 15) & (temperature <= 31);
qui replace  bLin = _b[ddayMod]*16 + _b[ddayHot]*(temperature-31) if (temperature > 31) & (temperature <= 39);
qui replace  bLin = _b[ddayMod]*16 + _b[ddayHot]*8 if (temperature > 39);


*-------------------------------------------------------------------------------------------;
* load regression results - dummy variables for 3C intervals;
estimates use ../resultSTATA/panelRegression/cotton_dummyTemperature6_39_year1950_2020_month4_10_clusterState;
gen tempDummy = .;
gen bDummy = .;
* lower truncated amount;
qui replace tempDummy = 0 in 1;
qui replace tempDummy = 3 in 2/3;
qui replace tempDummy = 6 in 4;
qui replace bDummy = _b[timeInt6_9] in 1/4;
* start dummies at 6C;
local i = 5;
forvalues b = 6(3)36 {;
	qui replace tempDummy = `b' in `i';
	qui replace tempDummy = `b'+3 in `=`i'+1';
	qui replace bDummy = _b[timeInt`b'_`=`b'+3'] in `i'/`=`i'+1';
	local i = `i' + 2;
};
* last interval is again open-ended;
qui replace bDummy = _b[timeInt39] in `i'/`=`i'+1';
qui replace tempDummy = 39 in `i';
qui replace tempDummy = 40 in `=`i'+1';


*-------------------------------------------------------------------------------------------;
* center graphs so exposure weighted average is zero (note: county fixed effects shift this up and down, all that matters is relative height of curve);
gen tempHist = .;
gen bHist = .;
local i = 1;
local avgLin = 0;
local avgCheby = 0;
local avgDummy = 0;
forvalues b = 0/39 {;
	qui replace tempHist = `b' in `=2*(`i'-1)+1';
	qui replace tempHist = `b'+1 in `=2*`i'';
	qui replace bHist = exposure[1,`i'] - exposure[1,`=`i'+1'] in `=2*(`i'-1)+1'/`=2*`i'';
	qui summ bLin if (temperature == `b'+0.5);
	local avgLin = `avgLin' + r(mean)*(exposure[1,`i'] - exposure[1,`=`i'+1'])/exposure[1,1];
	qui summ b if (temperature == `b'+0.5);
	local avgCheby = `avgCheby' + r(mean)*(exposure[1,`i'] - exposure[1,`=`i'+1'])/exposure[1,1];
	local nDummy = 2*floor((`b'+0.5)/3)+1;
	qui summ bDummy if (_n == `nDummy');
	local avgDummy = `avgDummy' + r(mean)*(exposure[1,`i'] - exposure[1,`=`i'+1'])/exposure[1,1];
	local i = `i' + 1;
};
* renormalize so exposure-weighted average is zero;
qui replace b = b - `avgCheby';
qui replace lb = lb - `avgCheby';
qui replace ub = ub - `avgCheby';
qui replace bLin = bLin - `avgLin';
qui replace bDummy = bDummy - `avgDummy';
* add zero as lower bound for histogram;
gen zero = 0 if (tempHist < .);


*-------------------------------------------------------------------------------------------;
* plot Figure 1 for cotton;
twoway 	(rarea lb ub temperature, color(gs8%30) yaxis(1))
		(line b temperature, lcolor(black) lpattern(solid) yaxis(1))
		(line bLin temperature, lcolor(red) lpattern(solid) yaxis(1))
		(line bDummy tempDummy, lcolor(blue) lpattern(solid) yaxis(1))
		(rarea bHist zero tempHist, color(green) yaxis(2))
	   ,
	   graphregion(fcolor(white)) legend(off)
	   xlabel(0(5)40) xtitle("Temperature (Celsius)")
	   ytitle("Log Cotton Yields (Pounds/Acre)", axis(1) color(black)) ylabel(-0.06(0.02)0.02, axis(1)) yscale(range(-0.08 0.02) axis(1))
	   ytitle("Exposure (Days)", axis(2) color(green)) ylabel(0(5)10, axis(2) labcolor(green) tlcolor(green)) yscale(range(0 50) axis(2))
	   yline(0, lpattern(dash) lcolor(gs6));


