# delimit ;

* save results in matrix;
*    columns 1:3: crops;
*    rows 1-7: model specification;
matrix rmse = J(6,3,0);

*-------------------------------------------------------------------------------------------;
* load out-of-sample predictions error - corn;
use ../resultSTATA/outOfSamplePrediction/corn_temperatureBounds10_29_36_year1950_2020_month3_8_draws1000, clear;
qui summ rmse1;
* rearrange order to graph (from piecewise linear, dummy, cheby to dummy, cheby, piecewise linear);
rename rmse2 rmse2a;
rename rmse3 rmse2;
rename rmse4 rmse3;
rename rmse2a rmse4;
local baselineModel = r(mean);
forvalues m = 1/6 {;
	qui summ rmse`=`m'+1';
	local newModel = r(mean);
	matrix rmse[`m',1] = (1-`newModel'/`baselineModel')*100;
};

*-------------------------------------------------------------------------------------------;
* load out-of-sample predictions error - soybeans;
use ../resultSTATA/outOfSamplePrediction/soybeans_temperatureBounds10_30_36_year1950_2020_month3_8_draws1000, clear;
* rearrange order to graph (from piecewise linear, dummy, cheby to dummy, cheby, piecewise linear);
rename rmse2 rmse2a;
rename rmse3 rmse2;
rename rmse4 rmse3;
rename rmse2a rmse4;
qui summ rmse1;
local baselineModel = r(mean);
forvalues m = 1/6 {;
	qui summ rmse`=`m'+1';
	local newModel = r(mean);
	matrix rmse[`m',2] = (1-`newModel'/`baselineModel')*100;
};

*-------------------------------------------------------------------------------------------;
* load out-of-sample predictions error - cotton;
use ../resultSTATA/outOfSamplePrediction/cotton_temperatureBounds15_31_39_year1950_2020_month4_10_draws1000, clear;
* rearrange order to graph (from piecewise linear, dummy, cheby to dummy, cheby, piecewise linear);
rename rmse2 rmse2a;
rename rmse3 rmse2;
rename rmse4 rmse3;
rename rmse2a rmse4;
qui summ rmse1;
local baselineModel = r(mean);
forvalues m = 1/6 {;
	qui summ rmse`=`m'+1';
	local newModel = r(mean);
	matrix rmse[`m',3] = (1-`newModel'/`baselineModel')*100;
};


clear;
set obs 24; 
gen x = 0.625 +(_n-1)*0.25;
gen rmseCorn = .;
gen rmseSoybeans = .;
gen rmseCotton = .;

forvalues m = 1/6 {;
	qui replace rmseCorn = rmse[`m',1] in `=(`m'-1)*4+1'/`=(`m'-1)*4+2';
	qui replace rmseSoybeans = rmse[`m',2] in `=(`m'-1)*4+2'/`=(`m'-1)*4+3';
	qui replace rmseCotton = rmse[`m',3] in `=(`m'-1)*4+3'/`=(`m'-1)*4+4';
};

local xLine0 -0.5;
local xLine1 -1;
twoway	(area rmseCorn x, color(red) cmissing(n))
		(area rmseSoybeans x, color(blue) cmissing(n))
		(area rmseCotton x, color(black) cmissing(n))
		,
		graphregion(fcolor(white)) legend(pos(12) cols(3) order(1 "Corn" 2 "Soybeans" 3 "Cotton"))
		ytitle("Percent Reduction in RMS") ylabel(0(2)16) yscale(range(0 16))
		xlabel(none) xtitle(" ") xscale(range(0.5 6.5))
		text(`xLine0' 1 "Step", placement(s) color(black) size(small))
		text(`xLine1' 1 "Function", placement(s) color(black) size(small))
		text(`xLine0' 2 "Polynomial", placement(s) color(black) size(small))
		text(`xLine1' 2 "(8th order)", placement(s) color(black) size(small))
		text(`xLine0' 3 "Piecewise", placement(s) color(black) size(small))
		text(`xLine1' 3 "Linear", placement(s) color(black) size(small))
		text(`xLine0' 4 "Monthly", placement(s) color(black) size(small))
		text(`xLine1' 4 "Averages", placement(s) color(black) size(small))
		text(`xLine0' 5 "Degree Days", placement(s) color(black) size(small))
		text(`xLine1' 5 "(Thom)", placement(s) color(black) size(small))
		text(`xLine0' 6 "Degree Days", placement(s) color(black) size(small))
		text(`xLine1' 6 "(Daily Mean)", placement(s) color(black) size(small))
		;


