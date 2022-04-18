# delimit ;


capture program drop bestLinearModel;
program define bestLinearModel;
	syntax, crop(string) lowerBound(integer) [monthMin(integer 4) monthMax(integer 8) yearMin(integer 1950) yearMax(integer 2020)];

	*-------------------------------------------------------------------------------------------;
	* load weather data and aggregate over season;
	use if (month >= `monthMin') & (month <= `monthMax') & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/weather_`crop';
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
	use if (crop == "`crop'") & (year >= `yearMin') & (year <= `yearMax') 
		using ../dataSTATA/yieldData, clear;
	merge 1:1 fips year using `weather', assert(3);
	drop _merge;
	
	* get additional variables;
	gen logYield = log(yield);
	gen state = floor(fips/1000);
	gen t = year-`yearMin';
	gen t2 = t^2;
	gen prec2 = prec^2;
	

	*-------------------------------------------------------------------------------------------;
	* loop over bounds;
	xtset fips year;
	* initialize matrix to save results;
	matrix r2 = J(`=40-`lowerBound'',2,0);
	forvalues b = `lowerBound'/39 {;
		gen ddayMod = dday`lowerBound'C - dday`b'C;
		gen ddayHot = dday`b'C - dday39C;
		reghdfe logYield ddayMod ddayHot prec prec2, a(fips i.state#c.t i.state#c.t2); 
		matrix r2[`=`b'-`lowerBound'+1',1] = `b';
		matrix r2[`=`b'-`lowerBound'+1',2] = e(r2);
		drop ddayMod ddayHot;
	};

	* save results;
	clear;
	svmat r2;
	rename r21 breakpoint;
	label var breakpoint "breakpoint b in piecewise linear function lower bound to b and above b";
	rename r22 r2; 
	label var r2 "R-square of model";
	save ../resultSTATA/bestLinearModel/`crop'_year`yearMin'_`yearMax'_month`monthMin'_`monthMax', replace;
	

end;


