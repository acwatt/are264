# delimit ;


capture program drop outOfSamplePrediction;
program define outOfSamplePrediction;
	syntax, crop(string) pieceLinLB(integer) pieceLinBreak(integer) tempLB(integer) tempUB(integer) 
			[monthMin(integer 3) monthMax(integer 8) yearMin(integer 1950) yearMax(integer 2020) nDraw(integer 1000)];

	*#############################################################################################;
	* load weather data and aggregate over season;
	* replicates MNS - quadratic in average temperature for January, April, July, October;
	use fips year month tMin tMax
		if inlist(month,1,4,7,10) & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/weather_`crop', clear;
	gen tAvgMonth = (tMin + tMax)/2;
	drop tMin tMax;
	reshape wide tAvg, i(fips year) j(month);
	tempfile weather;
	save `weather', replace;
	
	* merge with our new data;
	* load weather data and aggregate over season;
	use fips year month dday`pieceLinLB'C dday`pieceLinBreak'C dday`tempUB'C prec time* ddayThom8C ddayThom32C ddayThom34C ddayTavg8C ddayTavg32C 
		if (month >= `monthMin') & (month <= `monthMax') & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/weather_`crop', clear;
	* summ over season - degree days, time, and precipitation;
	collapse (sum) prec dday* time*, by(fips year);
	merge 1:1 fips year using `weather', assert(3);
	drop _merge;
	save `weather', replace;
	
	* merge in yield data;
	use crop fips latitude longitude year yield 
		if (crop == "`crop'") & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/yieldData, clear;
	merge 1:1 fips year using `weather', assert(3);

	* get additional variables;
	gen logYield = log(yield);
	gen state = floor(fips/1000);
	qui summ year;
	gen t = year-`yearMin';
	gen t2 = t^2;
	gen prec2 = prec^2;
	forvalues m = 1(3)10 {;
		gen tAvgMonth`m'_2 = tAvgMonth`m'^2;
	};
	
	*#############################################################################################;
	* prepare data - piecewise linear;
	gen ddayHot = dday`pieceLinBreak'C - dday`tempUB'C;
	gen ddayMod = dday`pieceLinLB'C - dday`pieceLinBreak'C;
	* replication using Thom's forumla (monthly instead of daily data to construct degree days);
	gen ddayHotThom = sqrt(ddayThom34C);
	gen ddayModThom = ddayThom8C - ddayThom32C;
	* replication using average daily temperature;
	gen ddayHotDailyAvg = ddayTavg8C - ddayTavg32C;
	gen ddayModDailyAvg = ddayTavg32C;

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
	

	*#############################################################################################;
	keep fips year logYield ddayMod ddayHot prec prec2 timeInt* cheby* tAvgMonth* ddayModThom ddayHotThom ddayModDailyAvg ddayHotDailyAvg state t t2;
	tempfile panel;
	save `panel', replace;
	local nObs =_N;

	* initialize matrix to store results;
	matrix rmse = J(`nDraw',7,0);

	* set random seed so results are reproducible;
	set seed 193432;

	* loop over random draws of years;
	forvalues n = 1/`nDraw' {;
		* randomly pick 85% of the years;
		clear;
		set obs `=`yearMax'-`yearMin'+1';
		gen year = `yearMin' -1 +_n;
		gen draw = runiform();
		sort draw;
		gen sampleYears = (_n <= 0.85*(`=`yearMax'-`yearMin'+1'));
		
		* merge back in panel data;
		merge 1:m year using `panel', assert(3);
		drop _merge;

		* --------------------------------------------------------------------------;
		* note: it is possible that for a county with few annual observations, none are in the 85% estimation sample;
		*       for those, take out county fixed effect in out-of-sample prediction as it is not identifed in estmation sample; 
		bysort fips: egen inSampleCheck = max(sampleYears);
		levelsof fips if (inSampleCheck == 0), local(fipsListDemean);
		drop inSampleCheck;
		

		* --------------------------------------------------------------------------;
		* baseline model - no weather;
		reghdfe logYield i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 		
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',1] = sqrt(r(mean));
		drop yHat d fe r rSquare;
		

		* --------------------------------------------------------------------------;
		* piecewise-linear model;
		reghdfe logYield ddayMod ddayHot prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',2] = sqrt(r(mean));
		drop yHat d fe r rSquare;

		
		* --------------------------------------------------------------------------;
		* dummy intervals;
		reghdfe logYield timeInt* prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',3] = sqrt(r(mean));
		drop yHat d fe r rSquare;


		* --------------------------------------------------------------------------;
		* chebychecv polynomials;
		reghdfe logYield cheby* prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',4] = sqrt(r(mean));
		drop yHat d fe r rSquare;


		* --------------------------------------------------------------------------;
		* four months;
		reghdfe logYield tAvgMonth* prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',5] = sqrt(r(mean));
		drop yHat d fe r rSquare;


		* --------------------------------------------------------------------------;
		* Thom's formula;
		reghdfe logYield ddayModThom ddayHotThom prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',6] = sqrt(r(mean));
		drop yHat d fe r rSquare;
		
		
		* --------------------------------------------------------------------------;
		* Degree days using monthly average;
		reghdfe logYield ddayModDailyAvg ddayHotDailyAvg prec prec2 i.state#c.t i.state#c.t2 if (sampleYears == 1), a(fips) resid; 
		* construct error out of sample;
		predict yHat if (sampleYears == 0), xb;
		predict d if (sampleYears == 1), d;
		bysort fips: egen fe = mean(d);
		gen r = logYield - yHat - fe if (sampleYears == 0);
		foreach f of local fipsListDemean {;
			qui summ r if (fips == `f') & (sampleYears == 0);
			qui replace r = r - r(mean) if (fips == `f') & (sampleYears == 0);
		};
		gen rSquare = r^2;
		summ rSquare if (sampleYears == 0);
		matrix rmse[`n',7] = sqrt(r(mean));
		drop yHat d fe r rSquare;
		
	};	
	
	* clear and reload rmse;
	clear;
	svmat rmse;
	
	label var rmse1 "baseline: no weather, just fixed effects and time trend";
	label var rmse2 "piecewise-linear in temperature / quadratic in precipitation";
	label var rmse3 "dummy for 3C ranges of temperature / quadratic in precipitation";
	label var rmse4 "Chebyshev polynomial in temperature / quadratic in precipitation";
	label var rmse5 "quadratic in avg. temperature for 4 months / quadratic in precipitation";
	label var rmse6 "degree days using Thom's formula on monthly temperature / quadratic in precip.";
	label var rmse7 "degree days using daily average temperature / quadratic in precipitation";

	save ../resultSTATA/outOfSamplePrediction/`crop'_temperatureBounds`pieceLinLB'_`pieceLinBreak'_`tempUB'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax'_draws`nDraw', replace;
		
end;


