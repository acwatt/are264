# delimit ;


capture program drop panelRegression;
program define panelRegression;
	syntax, crop(string) pieceLinLB(integer) pieceLinBreak(integer) tempLB(integer) tempUB(integer) 
			[monthMin(integer 3) monthMax(integer 8) yearMin(integer 1950) yearMax(integer 2020) conleyCutoff(integer 1000)];

	*-------------------------------------------------------------------------------------------;
	* load weather data and aggregate over season;
	use fips year month dday`pieceLinLB'C dday`pieceLinBreak'C dday`tempUB'C prec time*
		if (month >= `monthMin') & (month <= `monthMax') & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/weather_`crop', clear;
	* save variable labels;
	foreach v of varlist _all {;
		local l_`v': variable label `v';
	};

	* summ over season - degree days, time, and precipitation;
	collapse (sum) prec dday* time*, by(fips year);
	foreach v of varlist _all {;
		label var `v' "`l_`v''";
	};
	tempfile weather;
	save `weather', replace;
	
	*-------------------------------------------------------------------------------------------;
	* load yield data;
	use crop fips latitude longitude year yield 
		if (crop == "`crop'") & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/yieldData, clear;
	merge 1:1 fips year using `weather', assert(3);
	drop _merge;

	* get variables;
	gen logYield = log(yield);
	gen state = floor(fips/1000);
	qui summ year;
	gen t = year-`yearMin';
	gen t2 = t^2;
	gen prec2 = prec^2;

	
	*-------------------------------------------------------------------------------------------;
	* piecewise linear;
	gen ddayHot = dday`pieceLinBreak'C - dday`tempUB'C;
	gen ddayMod = dday`pieceLinLB'C - dday`pieceLinBreak'C;

	* part a: cluster by state;
	xtset fips year;
	reghdfe logYield ddayMod ddayHot prec prec2 i.state#c.t i.state#c.t2, a(fips) cluster(state) resid(reghdfe_resid); 
	estimates save ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState, replace;
	* save the residuals, needed for climate impacts;
	preserve;
	keep fips year reghdfe_resid;
	save ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, replace;
	restore;
	drop reghdfe_resid;

	* part b: conley;
	reghdfe logYield ddayMod ddayHot prec prec2 i.state#c.t i.state#c.t2, a(fips) resid(reghdfe_resid);
	estimates save ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_conley`conleyCutoff', replace;
	drop reghdfe_resid;


	*-------------------------------------------------------------------------------------------;
	* dummies for 3C range;
	forvalues t = `tempLB'(3)`=`tempUB'-3' {;
		if (`t' < 0) {;
			local tString0 Minus`=abs(`t')';
		}; else {;
			local tString0 `t';
		};
		if (`t' + 3 < 0) {;
			local tString1 Minus`=abs(`t'+3)';
		}; else {;
			local tString1 `=`t'+3';
		};
		gen timeInt`tString0'_`tString1' = time`tString0'C - time`tString1'C;
	};
	* open-ended top interval;
	gen timeInt`tempUB' = time`tempUB'C;

	* part a: cluster by state;
	xtset fips year;
	reghdfe logYield timeInt* prec prec2 i.state#c.t i.state#c.t2, a(fips) cluster(state) resid(reghdfe_resid); 
	estimates save ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState, replace;
	* save the residuals, needed for climate impacts;
	preserve;
	keep fips year reghdfe_resid;
	save ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, replace;
	restore;
	drop reghdfe_resid;

	* part b: conley;
	reghdfe logYield timeInt* prec prec2 i.state#c.t i.state#c.t2, a(fips) resid(reghdfe_resid); 
	conleyPanel, cutoff(`conleyCutoff');
	estimates save ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_conley`conleyCutoff', replace;
	drop reghdfe_resid;


	*-------------------------------------------------------------------------------------------;
	* chebychev polynomials: recursive Tn(x) = 2x*Tn-1(x) - Tn-2(x);
	matrix tCheby = J(`tempUB'-`tempLB',8,.);
	forvalues t = `tempLB'/`=`tempUB'-1' {;
		matrix tCheby[`t'-`tempLB'+1,1] = (`t'+0.5-`tempLB')/(`tempUB'-`tempLB')*2 - 1;
		matrix tCheby[`t'-`tempLB'+1,2] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,1] - 1;
		matrix tCheby[`t'-`tempLB'+1,3] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,2] - tCheby[`t'-`tempLB'+1,1];
		matrix tCheby[`t'-`tempLB'+1,4] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,3] - tCheby[`t'-`tempLB'+1,2];
		matrix tCheby[`t'-`tempLB'+1,5] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,4] - tCheby[`t'-`tempLB'+1,3];
		matrix tCheby[`t'-`tempLB'+1,6] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,5] - tCheby[`t'-`tempLB'+1,4];
		matrix tCheby[`t'-`tempLB'+1,7] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,6] - tCheby[`t'-`tempLB'+1,5];
		matrix tCheby[`t'-`tempLB'+1,8] = 2 * tCheby[`t'-`tempLB'+1,1] * tCheby[`t'-`tempLB'+1,7] - tCheby[`t'-`tempLB'+1,6];
	};
	forvalues p = 1/8 {;
		* open-ended top interval;
		gen cheby`p' = time`tempUB'C * tCheby[`=`tempUB'-1'-`tempLB'+1,`p'];
		* value mid-point by time in each 1C interval;
		forvalues t = `tempLB'/`=`tempUB'-2' {;
			if (`t' < 0) {;
				local tString0 Minus`=abs(`t')';
			}; else {;
				local tString0 `t';
			};
			if (`t' + 1 < 0) {;
				local tString1 Minus`=abs(`t'+1)';
			}; else {;
				local tString1 `=`t'+1';
			};
			qui replace cheby`p' = cheby`p' + (time`tString0'C - time`tString1'C)*tCheby[`t'-`tempLB'+1,`p'];
		};
	};
	
	* part a: cluster by state;
	xtset fips year;
	reghdfe logYield cheby1-cheby8 prec prec2 i.state#c.t i.state#c.t2, a(fips) cluster(state) resid(reghdfe_resid); 
	estimates save ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState, replace;
	* save the residuals, needed for climate impacts;
	preserve;
	keep fips year reghdfe_resid;
	save ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, replace;
	restore;
	drop reghdfe_resid;

	* part b: conley;
	reghdfe logYield cheby1-cheby8 prec prec2 i.state#c.t i.state#c.t2, a(fips) resid(reghdfe_resid); 
	conleyPanel, cutoff(`conleyCutoff');
	estimates save ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_conley`conleyCutoff', replace;
	drop reghdfe_resid;
	

end;


