# delimit ;


capture program drop r2_piecewiseLinear;
program define r2_piecewiseLinear;
	syntax, [yearMin(integer 1950) yearMax(integer 2020)];

	*-------------------------------------------------------------------------------------------;
	* load results of r2 as function of breakpoint - corn;
	use ../resultSTATA/bestLinearModel/corn_year`yearMin'_`yearMax'_month3_8, clear;
	
	qui summ r2;
	qui summ breakpoint if (r2 == r(max));
	local b_max = r(mean);
	
	twoway (line r2 breakpoint, lcolor(blue) lpattern(solid))
		   ,
		   graphregion(fcolor(white))
		   xlabel(10(5)40) xscale(range(10 40)) xtitle("Temperature Breakpoint (Celsius)")
		   xline(`b_max', lcolor(red) lpattern(dash))
		   text(0.861 `b_max' "Optimal breakpoint: `b_max'C", color(red) placement(ne))
		   ylabel(0.83(0.005)0.86) yscale(range(0.829 0.861)) ytitle("R-square");


	*-------------------------------------------------------------------------------------------;
	* load results of r2 as function of breakpoint - soybeans;
	use ../resultSTATA/bestLinearModel/soybeans_year`yearMin'_`yearMax'_month3_8, clear;
	
	qui summ r2;
	qui summ breakpoint if (r2 == r(max));
	local b_max = r(mean);
	
	twoway (line r2 breakpoint, lcolor(blue) lpattern(solid))
		   ,
		   graphregion(fcolor(white))
		   xlabel(10(5)40) xscale(range(10 40)) xtitle("Temperature Breakpoint (Celsius)")
		   xline(`b_max', lcolor(red) lpattern(dash))
		   text(0.791 `b_max' "Optimal breakpoint: `b_max'C", color(red) placement(ne))
		   ylabel(0.745(0.005)0.79) yscale(range(0.745 0.791)) ytitle("R-square");

	*-------------------------------------------------------------------------------------------;
	* load results of r2 as function of breakpoint - cotton;
	use ../resultSTATA/bestLinearModel/cotton_year`yearMin'_`yearMax'_month4_10, clear;
	
	qui summ r2;
	qui summ breakpoint if (r2 == r(max));
	local b_max = r(mean);
	
	twoway (line r2 breakpoint, lcolor(blue) lpattern(solid))
		   ,
		   graphregion(fcolor(white))
		   xlabel(15(5)40) xscale(range(15 40)) xtitle("Temperature Breakpoint (Celsius)")
		   xline(`b_max', lcolor(red) lpattern(dash))
		   text(0.566 `b_max' "Optimal breakpoint: `b_max'C", color(red) placement(ne))
		   ylabel(0.555(0.005)0.565) yscale(range(0.555 0.566)) ytitle("R-square");
	
end;


