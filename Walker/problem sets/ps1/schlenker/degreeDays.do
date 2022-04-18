# delimit;
clear all;
set more off;
set tracedepth 2;
* set trace on;


* #############################################################################;
* set parameters in following section;
* #############################################################################;
* part 1: degree days bound list;
local boundList 0 5 8 10 12 15 20 25 29 30 31 32 33 34;

* part 2: range of years in analysis (can be between 1950 and 2005);
local yearMin = 1957;
local yearMax = 1957;

* part 3: define growing season;
* first month and day of month used in a year;
local monthBegin = 4;  local dayMonthBegin =  1;
* last month and day of month used in a year;
local monthEnd   = 9; local dayMonthEnd   = 30;

* part 4: states that are included in analysis;
local stateList;
* NOTE: delete star at beginning of line to include a state;
* local stateList `stateList'  1; * Arkansas;
* local stateList `stateList'  4; * Arizona;
* local stateList `stateList'  5; * Arkansas;
* local stateList `stateList'  6; * California;
* local stateList `stateList'  8; * Colorado;
* local stateList `stateList'  9; * Connecticut;
* local stateList `stateList' 10; * Delaware;
* local stateList `stateList' 11; * District of Columbia;
* local stateList `stateList' 12; * Folorida;
* local stateList `stateList' 13; * Georgia;
* local stateList `stateList' 16; * Idaho;
* local stateList `stateList' 17; * Illinois;
* local stateList `stateList' 18; * Indiana;
* local stateList `stateList' 19; * Iowa;
* local stateList `stateList' 20; * Kansas;
* local stateList `stateList' 21; * Kentucky;
* local stateList `stateList' 22; * Louisiana;
* local stateList `stateList' 23; * Maine;
* local stateList `stateList' 24; * Maryland;
* local stateList `stateList' 25; * Massachusetts;
* local stateList `stateList' 26; * Michigan;
* local stateList `stateList' 27; * Minnesota;
* local stateList `stateList' 28; * Mississippi;
* local stateList `stateList' 29; * Missouri;
* local stateList `stateList' 30; * Montana;
* local stateList `stateList' 31; * Nebraska;
* local stateList `stateList' 32; * Nevada;
* local stateList `stateList' 33; * New Hampshire;
* local stateList `stateList' 34; * New Jersey;
* local stateList `stateList' 35; * New Mexico;
* local stateList `stateList' 36; * New York;
* local stateList `stateList' 37; * North Carolina;
* local stateList `stateList' 38; * North Dakota;
* local stateList `stateList' 39; * Ohio;
* local stateList `stateList' 40; * Oklahoma;
* local stateList `stateList' 41; * Oregon;
* local stateList `stateList' 42; * Pennsylvania;
* local stateList `stateList' 44; * Rhose Island;
* local stateList `stateList' 45; * South Carolina;
* local stateList `stateList' 46; * South Dakota;
* local stateList `stateList' 47; * Tennessee;
* local stateList `stateList' 48; * Texas;
* local stateList `stateList' 49; * Utah;
* local stateList `stateList' 50; * Vermont;
* local stateList `stateList' 51; * Virginia;
* local stateList `stateList' 53; * Washington;
* local stateList `stateList' 54; * West Virginia;
* local stateList `stateList' 55; * Wisconsin;
* local stateList `stateList' 56; * Wyoming;
* #############################################################################;


* #############################################################################;
* #############################################################################;
* #############################################################################;
* program code from here on, no need to change anything;
* #############################################################################;
* load grid information;
use ../metaData/cropArea, clear;
gen longitude = -125 + mod(gridNumber-1,1405)/24;
label var longitude "longitude of grid centroid (decimal degrees)";
gen latitude  = 49.9375+1/48 - ceil(gridNumber/1405)/24;
label var latitude  "latitude of grid centroid (decimal degrees)";
merge 1:1 gridNumber using ../metaData/linkGridnumberFIPS;
drop _merge;
compress;
sort gridNumber;
save gridInfo, replace;

******************************************************************************;
* loop over states;
local appendID = 0;
foreach s of local stateList {;
    * loop over years;
    forvalues t = `yearMin'/`yearMax' {;
        ******************************************************************************;
        display("Working on state `s' in year `t'");
        * load data;
    	use if (dateNum >= mdy(`monthBegin',`dayMonthBegin',year(dateNum)))
    			& (dateNum <= mdy(`monthEnd',`dayMonthEnd',year(dateNum)))
    			using ../rawDataByYear/year`t'/state`s', clear;
        gen tAvg = (tMin+tMax)/2;
        label var tAvg "average temperature - average of minimum and maximum temperature";

        * merge in grid weights;
        * make sure there is no (_merge == 1) as gridInfo should have data for each gridNumber;
        qui merge m:1 gridNumber using gridInfo, assert(2 3);

        * keep only were weights are nonzero;
        qui keep if (_merge == 3) & (cropArea > 0);
        drop _merge;

        if (_N > 0) {;
            ******************************************************************************;
            * generate degree days;
            foreach b of local boundList {;
                * create label for negative bounds by adding Minus;
                if (`b' < 0) {; 
                	local b2 = abs(`b'); 
                	local bLabel Minus`b2'; 
                }; else {; 
                	local bLabel `b'; 
                };

                * default case 1: tMax <= bound;
                qui gen dday`bLabel'C = 0;

                * case 2: bound <= tMin;
                qui replace dday`bLabel'C = tAvg - `b' if (`b' <= tMin);

                * case 3: tMin < bound < tMax;
                qui gen tempSave = acos( (2*`b'-tMax-tMin)/(tMax-tMin) );
                qui replace dday`bLabel'C = ( (tAvg-`b')*tempSave + (tMax-tMin)*sin(tempSave)/2 )/_pi 
                        if ( (tMin < `b') & (`b' < tMax) );
                drop tempSave;
            };

            ******************************************************************************;
			* save variable labels;
			foreach v in tMin tMax tAvg prec {;
				local l_`v': variable label `v';
			};
             
            * collapse by aggregation level;
            collapse tMin tMax tAvg prec dday* [aweight=cropArea], by(fips dateNum);

            * collapse by year;
            gen year = year(dateNum);
            collapse (sum) prec dday* (mean) tMin tMax tAvg, by(fips year);

			* label variables;
			foreach v in tMin tMax tAvg prec {;
				label var `v' "`l_`v''";
			};
			foreach b of local boundList {;
				* create label for negative bounds by adding Minus;
				if (`b' < 0) {; 
					local b2 = abs(`b'); 
					local bLabel Minus`b2'; 
				}; else {; 
					local bLabel `b'; 
				};
				label var dday`bLabel'C "total degree days above `b'C (Celcius and days)";
			};
			label var year "year";

            ******************************************************************************;
            * save data;
            if (`appendID' > 0) {;
	            append using ddayByYearandFips_cropAreaWeighted;
	        };
	        local appendID = 1;
            sort fips year;
            qui save ddayByYearandFips_cropAreaWeighted, replace;
        };
    };
};

* clear up files;
! rm gridInfo.dta;
