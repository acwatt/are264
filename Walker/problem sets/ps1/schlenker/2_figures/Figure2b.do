# delimit ;

* save results in matrix;
*    columns 1:3: crops;
*    rows 1-7: model specification;
matrix impact = J(12,3,0);
matrix lb = J(12,3,0);
matrix ub = J(12,3,0);

*-------------------------------------------------------------------------------------------;
* load climate impact - corn;
use ../resultSTATA/climateImpactUniform/corn_temperatureBounds10_29_36_year1950_2020_month3_8_clusterState_draws1000, clear;
local c = 1;
forvalues t = 5/8 {;
	matrix impact[`=`t'-4',1] = impactDummy`t'[1];
	centile impactDummy`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'-4',1] = r(c_1);
	matrix ub[`=`t'-4',1] = r(c_2);

	matrix impact[`=`t'-4',2] = impactCheby`t'[1];
	centile impactCheby`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'-4',2] = r(c_1);
	matrix ub[`=`t'-4',2] = r(c_2);

	matrix impact[`=`t'-4',3] = impactLin`t'[1];
	centile impactLin`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'-4',3] = r(c_1);
	matrix ub[`=`t'-4',3] = r(c_2);
};


*-------------------------------------------------------------------------------------------;
* load climate impact - soybeans;
use ../resultSTATA/climateImpactUniform/soybeans_temperatureBounds10_30_36_year1950_2020_month3_8_clusterState_draws1000, clear;
local c = 1;
forvalues t = 5/8 {;
	matrix impact[`t',1] = impactDummy`t'[1];
	centile impactDummy`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`t',1] = r(c_1);
	matrix ub[`t',1] = r(c_2);

	matrix impact[`t',2] = impactCheby`t'[1];
	centile impactCheby`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`t',2] = r(c_1);
	matrix ub[`t',2] = r(c_2);

	matrix impact[`t',3] = impactLin`t'[1];
	centile impactLin`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`t',3] = r(c_1);
	matrix ub[`t',3] = r(c_2);
};

*-------------------------------------------------------------------------------------------;
* load climate impact - cotton;
use ../resultSTATA/climateImpactUniform/cotton_temperatureBounds15_31_39_year1950_2020_month4_10_clusterState_draws1000, clear;
local c = 1;
forvalues t = 5/8 {;
	matrix impact[`=`t'+4',1] = impactDummy`t'[1];
	centile impactDummy`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'+4',1] = r(c_1);
	matrix ub[`=`t'+4',1] = r(c_2);

	matrix impact[`=`t'+4',2] = impactCheby`t'[1];
	centile impactCheby`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'+4',2] = r(c_1);
	matrix ub[`=`t'+4',2] = r(c_2);

	matrix impact[`=`t'+4',3] = impactLin`t'[1];
	centile impactLin`t' if (_n > 1), centile(2.5 97.5);
	matrix lb[`=`t'+4',3] = r(c_1);
	matrix ub[`=`t'+4',3] = r(c_2);
};


*-------------------------------------------------------------------------------------------;

clear;
set obs 36; 
gen x = ceil(_n/3) + mod(_n-1,3)*0.15 - 0.15;
qui replace x = x+1 if (x > 4.5);
qui replace x = x+1 if (x > 9.5);

qui gen b = .;
qui gen lb = .;
qui gen ub = .;
qui gen m = .;

forvalues m = 1/12 {;
	forvalues j = 1/3 {;
		qui replace b = impact[`m',`j']*100 in `=(`m'-1)*3+`j'';
		qui replace lb = lb[`m',`j']*100 in `=(`m'-1)*3+`j'';
		qui replace ub = ub[`m',`j']*100 in `=(`m'-1)*3+`j'';
		qui replace m = `j' in `=(`m'-1)*3+`j'';
	};
};

local xLine0 -104;
local xLine1 -109;
twoway	(rcap lb ub x if (m==1), color(blue))
		(scatter  b x if (m==1), color(blue) msym(x))
		(rcap lb ub x if (m==2), color(black))
		(scatter  b x if (m==2), color(black) msym(x))
		(rcap lb ub x if (m==3), color(red))
		(scatter  b x if (m==3), color(red) msym(x))
		,
		graphregion(fcolor(white)) legend(pos(12) cols(3) order(1 "Step Function" 3 "Chebyshev Polynomial" 5 "Piecewise Linear"))
		ytitle("Impact (Percent)") ylabel(-100(20)0) yscale(range(-100 14))
		yline(0, lcolor(black) lpattern(solid))
		xlabel(none) xtitle(" ") xscale(range(0.5 6.5)) xsize(7)
		text(`xLine0'  1 "+5C", placement(s) color(black) size(medium))
		text(`xLine0'  2 "+6C", placement(s) color(black) size(medium))
		text(`xLine0'  3 "+7C", placement(s) color(black) size(medium))
		text(`xLine0'  4 "+8C", placement(s) color(black) size(medium))
		text(`xLine0'  6 "+5C", placement(s) color(black) size(medium))
		text(`xLine0'  7 "+6C", placement(s) color(black) size(medium))
		text(`xLine0'  8 "+7C", placement(s) color(black) size(medium))
		text(`xLine0'  9 "+8C", placement(s) color(black) size(medium))
		text(`xLine0' 11 "+5C", placement(s) color(black) size(medium))
		text(`xLine0' 12 "+6C", placement(s) color(black) size(medium))
		text(`xLine0' 13 "+7C", placement(s) color(black) size(medium))
		text(`xLine0' 14 "+8C", placement(s) color(black) size(medium))
		text(`xLine1'  2.5 "Corn", placement(s) color(black) size(medium))
		text(`xLine1'  7.5 "Soybeans", placement(s) color(black) size(medium))
		text(`xLine1' 12.5 "Cotton", placement(s) color(black) size(medium))
		;


