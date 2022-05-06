# delimit ;


*-------------------------------------------------------------------------------------------;
* get time temperatures are within each 1C interval for exposure histogram;
use fips year month time*
	if (month >= 3) & (month <= 8) & (year >= 1950) & (year <= 2020) 
	using ../dataSTATA/weather_soybeans, clear;
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
estimates use ../resultSTATA/panelRegression/soybeans_chebyTemperature-3_36_year1950_2020_month3_8_clusterState;

* transform temperature to [-1,1]- note: Chebyshev polynomials was between -3 and 36C;
gen tCheby = (temperature+3)/39*2 - 1 if (temperature <= 36);
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
forvalues i = 1/`=36*4+1' {;
	qui lincom    `=cheby1[`i']' * _b[cheby1] + `=cheby2[`i']' * _b[cheby2] + `=cheby3[`i']' * _b[cheby3] + `=cheby4[`i']' * _b[cheby4]
				+ `=cheby5[`i']' * _b[cheby5] + `=cheby6[`i']' * _b[cheby6] + `=cheby7[`i']' * _b[cheby7] + `=cheby8[`i']' * _b[cheby8]; 
	qui replace b = r(estimate) in `i';
	qui replace lb = r(estimate)-1.96*r(se) in `i';
	qui replace ub = r(estimate)+1.96*r(se) in `i';
};
* add constant part above upper bound of 36C - results were truncated there as very little temperature mass above that threshold;
qui replace b = r(estimate) if (temperature > 36);
qui replace lb = r(estimate)-1.96*r(se) if (temperature > 36);
qui replace ub = r(estimate)+1.96*r(se) if (temperature > 36);


*-------------------------------------------------------------------------------------------;
* load regression results - piecewise linear;
estimates use ../resultSTATA/panelRegression/soybeans_piecewiseLinear10_30_36_year1950_2020_month3_8_clusterState;
gen bLin = 0 if (temperature <= 10);
qui replace  bLin = _b[ddayMod]*(temperature-10) if (temperature > 10) & (temperature <= 30);
qui replace  bLin = _b[ddayMod]*20 + _b[ddayHot]*(temperature-30) if (temperature > 30) & (temperature <= 36);
qui replace  bLin = _b[ddayMod]*20 + _b[ddayHot]*6 if (temperature > 36);

*-------------------------------------------------------------------------------------------;
* load regression results - dummy variables for 3C intervals;
estimates use ../resultSTATA/panelRegression/soybeans_dummyTemperature-3_36_year1950_2020_month3_8_clusterState;
gen tempDummy = .;
gen bDummy = .;
local i = 1;
forvalues b = 0(3)33 {;
	qui replace tempDummy = `b' in `i';
	qui replace tempDummy = `b'+3 in `=`i'+1';
	qui replace bDummy = _b[timeInt`b'_`=`b'+3'] in `i'/`=`i'+1';
	local i = `i' + 2;
};
* last interval is again open-ended;
qui replace bDummy = _b[timeInt36] in `i'/`=`i'+2';
qui replace tempDummy = 36 in `i';
qui replace tempDummy = 39 in `=`i'+1';
qui replace tempDummy = 40 in `=`i'+2';

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
* plot Figure 1 for soybeans;
twoway 	(rarea lb ub temperature, color(gs8%30) yaxis(1))
		(line b temperature, lcolor(black) lpattern(solid) yaxis(1))
		(line bLin temperature, lcolor(red) lpattern(solid) yaxis(1))
		(line bDummy tempDummy, lcolor(blue) lpattern(solid) yaxis(1))
		(rarea bHist zero tempHist, color(green) yaxis(2))
	   ,
	   graphregion(fcolor(white)) legend(off)
	   xlabel(0(5)40) xtitle("Temperature (Celsius)")
	   ytitle("Log Soybean Yields (Bushel/Acre)", axis(1) color(black)) ylabel(-0.06(0.02)0.02, axis(1)) yscale(range(-0.08 0.02) axis(1))
	   ytitle("Exposure (Days)", axis(2) color(green)) ylabel(0(5)10, axis(2) labcolor(green) tlcolor(green)) yscale(range(0 50) axis(2))
	   yline(0, lpattern(dash) lcolor(gs6));


