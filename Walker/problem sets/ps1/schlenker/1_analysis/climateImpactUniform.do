# delimit ;

capture program drop timeIntervalLabelDummy;
program define timeIntervalLabelDummy, rclass;
	syntax, tempLowerBound(integer);

	* note: this script adds Minus to temperature label if bound is negative;
	if (`tempLowerBound' < 0) {;
		* note: data does not extend below -5C, so need to truncate there;
		local tString0 Minus`=abs(max(`tempLowerBound',-5))';
	}; else {;
		local tString0 `tempLowerBound';
	};
	if (`tempLowerBound' + 3 < 0) {;
		local tString1 Minus`=abs(max(`tempLowerBound'+3,-5))';
	}; else {;
		local tString1 `=`tempLowerBound'+3';
	};
	return local intLB `tString0';
	return local intUB `tString1';
end;


capture program drop timeIntervalLabelCheby;
program define timeIntervalLabelCheby, rclass;
	syntax, tempLowerBound(integer);

	* note: this script adds Minus to temperature label if bound is negative;
	if (`tempLowerBound' < 0) {;
		* note: data does not extend below -5C, so need to truncate there;
		local tString0 Minus`=abs(max(`tempLowerBound',-5))';
	}; else {;
		local tString0 `tempLowerBound';
	};
	if (`tempLowerBound' + 1 < 0) {;
		local tString1 Minus`=abs(max(`tempLowerBound'+1,-5))';
	}; else {;
		local tString1 `=`tempLowerBound'+1';
	};
	return local intLB `tString0';
	return local intUB `tString1';
end;


capture program drop climateImpactUniform;
program define climateImpactUniform;
	syntax, crop(string) pieceLinLB(integer) pieceLinBreak(integer) tempLB(integer) tempUB(integer) 
			[monthMin(integer 3) monthMax(integer 8) yearMin(integer 1950) yearMax(integer 2020) nDraw(integer 1000)];

	*---------------------------------------------------------------------------------------------------------------------;
	* get 1000 draws of variance-covariance matrix usedto approximate distribution of impacts;
	*---------------------------------------------------------------------------------------------------------------------;
	* reload results - piecewise linear;
	estimates use ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	matrix b = e(b);
	matrix V = e(V);
	matrix V2 = V[`=rownumb(V,"ddayMod")'..`=rownumb(V,"ddayHot")',`=colnumb(V,"ddayMod")'..`=colnumb(V,"ddayHot")'];
	matrix b2 = b[1,`=colnumb(b,"ddayMod")'..`=colnumb(V,"ddayHot")'];
	drawnorm bHot bMod, n(`nDraw') cov(V2) means(b2);
	mkmat bHot bMod, matrix(bSampleLin);
	matrix drop b b2 V V2;


	*---------------------------------------------------------------------------------------------------------------------;
	* reload results - dummies for 3C intervals;
	estimates use ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	matrix b = e(b);
	matrix V = e(V);
	if (`tempLB' < 0) {;
		local tString0 Minus`=abs(`tempLB')';
	}; else {;
		local tString0 `tempLB';
	};
	if (`tempLB' + 3 < 0) {;
		local tString1 Minus`=abs(`tempLB'+3)';
	}; else {;
		local tString1 `=`tempLB'+3';
	};
	local n1 = rownumb(V,"timeInt`tString0'_`tString1'");
	local n2 = rownumb(V,"timeInt`tempUB'");
	matrix V2 = V[`n1'..`n2',`n1'..`n2'];
	matrix b2 = b[1,`n1'..`n2'];
	local tIntList;
	forvalues i = 1/`=`n2'-`n1'+1' {;
		local tIntList `tIntList' tInt`i';
	};
	drawnorm `tIntList', n(`nDraw') cov(V2) means(b2);
	mkmat `tIntList', matrix(bSampleDummy);
	matrix drop b b2 V V2;


	*---------------------------------------------------------------------------------------------------------------------;
	* reload results - Chebyshev polynomials;
	estimates use ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	matrix b = e(b);
	matrix V = e(V);
	matrix V2 = V[`=rownumb(V,"cheby1")'..`=rownumb(V,"cheby8")',`=colnumb(V,"cheby1")'..`=colnumb(V,"cheby8")'];
	matrix b2 = b[1,`=colnumb(b,"cheby1")'..`=colnumb(V,"cheby8")'];
	drawnorm bCheby1 bCheby2 bCheby3 bCheby4 bCheby5 bCheby6 bCheby7 bCheby8, n(`nDraw') cov(V2) means(b2);
	mkmat bCheby1 bCheby2 bCheby3 bCheby4 bCheby5 bCheby6 bCheby7 bCheby8, matrix(bSampleCheby);
	matrix drop b b2 V V2;


	*---------------------------------------------------------------------------------------------------------------------;
	* initializ matrix where results are stored;
	*   note: first row is point estimate, the remaining nDraw rows are to get uncertainty;
	matrix impactLin = J(`=`nDraw'+1',8,0);
	matrix impactDummy = J(`=`nDraw'+1',8,0);
	matrix impactCheby = J(`=`nDraw'+1',8,0);


	*---------------------------------------------------------------------------------------------------------------------;
	* get data;
	*---------------------------------------------------------------------------------------------------------------------;
	* load weather data and aggregate over season;
	use fips year month dday*C prec time*
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
	tempfile panelData;
	save `panelData', replace;
	
	*---------------------------------------------------------------------------------------------------------------------;
	* load yield data;
	use crop fips latitude longitude year yield area 
		if (crop == "`crop'") & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/yieldData, clear;
	merge 1:1 fips year using `panelData', assert(3);
	drop _merge;

	* get variables;
	gen logYield = log(yield);
	gen state = floor(fips/1000);
	qui summ year;
	gen t = year-`yearMin';
	gen t2 = t^2;
	gen prec2 = prec^2;
	save `panelData', replace;

	
	*---------------------------------------------------------------------------------------------------------------------;
	* climate impact piecewise linear;
	*---------------------------------------------------------------------------------------------------------------------;
	gen ddayHot = dday`pieceLinBreak'C - dday`tempUB'C;
	gen ddayMod = dday`pieceLinLB'C - dday`pieceLinBreak'C;
	
	* reload regression estimates;
	estimates use ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	merge 1:1 fips year using ../resultSTATA/panelRegression/`crop'_piecewiseLinear`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, assert(3);
	drop _merge;

	* note: climate mpacts are relative to 1960-1989 baseline - keep those years;
	keep if (year >= 1960) & (year < 1989);
	* note: some observations are not included as they are singletons and hence the county fixed effect completly absorbs them;
	*       they are not counted here;
	predict yHat0, xbd;
	qui gen prod0 = exp(yHat0)*area;
	preserve;
		
	* step 1: point estimate;	
	* loop over uniform temperature increases;
	forvalues t = 1/8 {;
		* predicted new log yield;
		qui gen yHat`t' = yHat0 + _b[ddayHot]*(dday`=`pieceLinBreak'-`t''C - dday`=`tempUB'-`t''C - dday`pieceLinBreak'C + dday`tempUB'C)
							+ _b[ddayMod]*(dday`=`pieceLinLB'-`t''C - dday`=`pieceLinBreak'-`t''C - dday`pieceLinLB'C + dday`pieceLinBreak'C);
		* predicted new total production (use area-weighting for impacts);
		qui gen prod`t' = exp(yHat`t')*area;
	};
	* sum over all counties in a year and then average years;
	collapse (sum) prod*, by(year);
	collapse (mean) prod*;
	* calculate yield declines (area-weighted yied declines are just production declines as area remains fixed);
	forvalues t = 1/8 {;
		matrix impactLin[1,`t'] = prod`t'[1]/prod0[1]-1;
	};
	
	* step 2: loop over random draws and calculate impact again for confidence interval;
	forvalues i = 1/`nDraw' {;
		restore; preserve;
		* loop over uniform temperature increases;
		forvalues t = 1/8 {;
			* predicted new log yield;
			qui gen yHat`t' = yHat0 + bSampleLin[`i',2]*(dday`=`pieceLinBreak'-`t''C - dday`=`tempUB'-`t''C - dday`pieceLinBreak'C + dday`tempUB'C)
								+ bSampleLin[`i',1]*(dday`=`pieceLinLB'-`t''C - dday`=`pieceLinBreak'-`t''C - dday`pieceLinLB'C + dday`pieceLinBreak'C);
			* predicted new total production (use area-weighting for impacts);
			qui gen prod`t' = exp(yHat`t')*area;
		};
		* sum over all counties in a year and then average years;
		collapse (sum) prod*, by(year);
		collapse (mean) prod*;
		* calculate yield declines (area-weighted yield declines are just production declines as area remains fixed);
		forvalues t = 1/8 {;
			matrix impactLin[`=`i'+1',`t'] = prod`t'[1]/prod0[1]-1;
		};
		display("Piecewise linear model - finished loop `i'");		
	};
	restore;


	*---------------------------------------------------------------------------------------------------------------------;
	* climate impact for 3c Dummy specification;
	*---------------------------------------------------------------------------------------------------------------------;
	use `panelData', clear;
	forvalues t = `tempLB'(3)`=`tempUB'-3' {;
		timeIntervalLabelDummy, tempLowerBound(`t');
		gen timeInt`=r(intLB)'_`=r(intUB)' = time`=r(intLB)'C - time`=r(intUB)'C;
	};
	* open-ended top interval;
	gen timeInt`tempUB' = time`tempUB'C;

	* reload regression estimates;
	estimates use ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	merge 1:1 fips year using ../resultSTATA/panelRegression/`crop'_dummyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, assert(3);
	drop _merge;

	* note: climate mpacts are relative to 1960-1989 baseline - keep those years;
	keep if (year >= 1960) & (year < 1989);
	* note: some observations are not included as they are singletons and hence the county fixed effect completly absorbs them;
	*       they are not counted here;
	predict yHat0, xbd;
	qui gen prod0 = exp(yHat0)*area;
	preserve;

	* step 1: point estimate;	
	* loop over uniform temperature increases;
	forvalues t = 1/8 {;
		* predicted new log yield;
		qui gen yHat`t' = yHat0;
		forvalues tBound = `tempLB'(3)`=`tempUB'-3' {;
			timeIntervalLabelDummy, tempLowerBound(`=`tBound'-`t'');
			gen timeIntNew = time`=r(intLB)'C - time`=r(intUB)'C;
			timeIntervalLabelDummy, tempLowerBound(`tBound');
			gen timeInt = time`=r(intLB)'C - time`=r(intUB)'C;
			qui replace yHat`t' = yHat`t' + _b[timeInt`=r(intLB)'_`=r(intUB)'] * (timeIntNew - timeInt);
			drop timeIntNew timeInt;
		};
		qui replace yHat`t' = yHat`t' + _b[timeInt`tempUB'] * (time`=`tempUB'-`t''C - time`tempUB'C);
		qui gen prod`t' = exp(yHat`t')*area;
	};
	* sum over all counties in a year and then average years;
	collapse (sum) prod*, by(year);
	collapse (mean) prod*;
	* calculate yield declines (area-weighted yield declines are just production declines as area remains fixed);
	forvalues t = 1/8 {;
		matrix impactDummy[1,`t'] = prod`t'[1]/prod0[1]-1;
	};

	* step 2: loop over random draws and calculate impact again for confidence interval;
	forvalues i = 1/`nDraw' {;
		restore; preserve;
		* loop over uniform temperature increases;
		forvalues t = 1/8 {;
			* predicted new log yield;
			qui gen yHat`t' = yHat0;
			local c = 1;
			forvalues tBound = `tempLB'(3)`=`tempUB'-3' {;
				timeIntervalLabelDummy, tempLowerBound(`=`tBound'-`t'');
				gen timeIntNew = time`=r(intLB)'C - time`=r(intUB)'C;
				timeIntervalLabelDummy, tempLowerBound(`tBound');
				gen timeInt = time`=r(intLB)'C - time`=r(intUB)'C;
				qui replace yHat`t' = yHat`t' + bSampleDummy[`i',`c'] * (timeIntNew - timeInt);
				drop timeIntNew timeInt;
				local c = `c'+1;
			};
			qui replace yHat`t' = yHat`t' + bSampleDummy[`i',`c'] * (time`=`tempUB'-`t''C - time`tempUB'C);
			qui gen prod`t' = exp(yHat`t')*area;
		};
		* sum over all counties in a year and then average years;
		collapse (sum) prod*, by(year);
		collapse (mean) prod*;
		* calculate yield declines (area-weighted yield declines are just production declines as area remains fixed);
		forvalues t = 1/8 {;
			matrix impactDummy[`=`i'+1',`t'] = prod`t'[1]/prod0[1]-1;
		};
		display("Dummy variable model - finished loop `i'");		
	};
	restore;
	*/;


	*---------------------------------------------------------------------------------------------------------------------;
	* climate impact for chebyshev Polynomials;
	*---------------------------------------------------------------------------------------------------------------------;
	use `panelData', clear;
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
			timeIntervalLabelCheby, tempLowerBound(`t');
			qui replace cheby`p' = cheby`p' + (time`=r(intLB)'C - time`=r(intUB)'C)*tCheby[`t'-`tempLB'+1,`p'];
		};
	};

	* reload regression estimates;
	estimates use ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState;
	merge 1:1 fips year using ../resultSTATA/panelRegression/`crop'_chebyTemperature`tempLB'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_residuals, assert(3);
	drop _merge;

	* note: climate mpacts are relative to 1960-1989 baseline - keep those years;
	keep if (year >= 1960) & (year < 1989);
	* note: some observations are not included as they are singletons and hence the county fixed effect completly absorbs them;
	*       they are not counted here;
	predict yHat0, xbd;
	qui gen prod0 = exp(yHat0)*area;
	preserve;

	* step 1: point estimate;	
	* loop over uniform temperature increases;
	forvalues t = 1/8 {;
		* predicted new log yield;
		qui gen yHat`t' = yHat0;
		forvalues p = 1/8 {;
			* open-ended top interval;
			gen chebyNew`p' = time`=`tempUB'-`t''C * tCheby[`=`tempUB'-1'-`tempLB'+1,`p'];
			* value mid-point by time in each 1C interval;
			forvalues tBound = `tempLB'/`=`tempUB'-2' {;
				timeIntervalLabelCheby, tempLowerBound(`=`tBound'-`t'');
				qui replace chebyNew`p' = chebyNew`p' + (time`=r(intLB)'C - time`=r(intUB)'C)*tCheby[`tBound'-`tempLB'+1,`p'];
			};
			qui replace yHat`t' = yHat`t' + _b[cheby`p'] * (chebyNew`p' - cheby`p');
			drop chebyNew`p';
		};
		qui gen prod`t' = exp(yHat`t')*area;
	};
	* sum over all counties in a year and then average years;
	collapse (sum) prod*, by(year);
	collapse (mean) prod*;
	* calculate yield declines (area-weighted yield declines are just production declines as area remains fixed);
	forvalues t = 1/8 {;
		matrix impactCheby[1,`t'] = prod`t'[1]/prod0[1]-1;
	};

	* step 2: loop over random draws and calculate impact again for confidence interval;
	forvalues i = 1/`nDraw' {;
		restore; preserve;
		* loop over uniform temperature increases;
		forvalues t = 1/8 {;
			* predicted new log yield;
			qui gen yHat`t' = yHat0;
			forvalues p = 1/8 {;
				* open-ended top interval;
				gen chebyNew`p' = time`=`tempUB'-`t''C * tCheby[`=`tempUB'-1'-`tempLB'+1,`p'];
				* value mid-point by time in each 1C interval;
				forvalues tBound = `tempLB'/`=`tempUB'-2' {;
					timeIntervalLabelCheby, tempLowerBound(`=`tBound'-`t'');
					qui replace chebyNew`p' = chebyNew`p' + (time`=r(intLB)'C - time`=r(intUB)'C)*tCheby[`tBound'-`tempLB'+1,`p'];
				};
				qui replace yHat`t' = yHat`t' + bSampleCheby[`i',`p'] * (chebyNew`p' - cheby`p');
				drop chebyNew`p';
			};
			qui gen prod`t' = exp(yHat`t')*area;
		};
		* sum over all counties in a year and then average years;
		collapse (sum) prod*, by(year);
		collapse (mean) prod*;
		* calculate yield declines (area-weighted yield declines are just production declines as area remains fixed);
		forvalues t = 1/8 {;
			matrix impactCheby[`=`i'+1',`t'] = prod`t'[1]/prod0[1]-1;
		};
		display("Chebyshev model - finished loop `i'");		
	};
	restore;

	
	*---------------------------------------------------------------------------------------------------------------------;
	* load results and save;
	*---------------------------------------------------------------------------------------------------------------------;
	clear;
	svmat impactLin; 
	svmat impactDummy; 
	svmat impactCheby; 

	save ../resultSTATA/climateImpactUniform/`crop'_temperatureBounds`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_clusterState_draws`nDraw';

end;


