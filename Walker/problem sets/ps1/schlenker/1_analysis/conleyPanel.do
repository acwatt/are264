# delimit ;


capture program drop conleyPanel;
program define conleyPanel, eclass;
	syntax, cutoff(real);

	* save data frame;
	preserve;
	keep if e(sample);
	
	* get number of variables in regression;
	local K = rowsof(e(V));
	local df = e(df_r);

	* get independent variables and error term in matrix form;
	predict e if e(sample), residuals;
	
	* get list of exogenous variables;
	local varList = e(indepvars);
	* replace constant in variable list with contantTerm;
	local varList = subinstr("`varList'","_cons","constantTerm", .);
	gen constantTerm = 1;
	

	* -----------------------------------------------------------------------------------;
	* pass data set on to MATA;
	mata: e = st_data(.,"e");
	mata: X = st_data(.,"`varList'");
	mata: latLong = st_data(.,"latitude longitude");
	mata: year = st_data(.,"year");
	mata: XeeX = J(`K',`K',0);

	* -----------------------------------------------------------------------------------;
	* loop over different years - and derive within year spatial correlation;
	qui levelsof year, local(yearList);
	foreach y of local yearList {;
		mata: rowsYear = selectindex(year :== `y');
		mata: eYear = e[rowsYear,1];
		mata: Xyear = X[rowsYear,.];
		mata: obsYear = rows(Xyear);
		mata: latLongYear = latLong[rowsYear,.]*pi()/180;
		mata: mata drop rowsYear;
		
		* get number of observations for year;
		mata: st_numscalar("obsYear",obsYear);
		local obsYear = obsYear;

		* loop over observations;
		forvalues i = 1/`obsYear' {;
			* step a: get non-parametric weight, Bartlett window of Newey-West;
			mata: weight = J(`obsYear',1,1) - 20038/pi()/`cutoff' * 
								acos( (sin(latLongYear[.,1]):*sin(J(`obsYear',1,latLongYear[`i',1]))) 
									+ (cos(latLongYear[.,1]):*cos(J(`obsYear',1,latLongYear[`i',1])):*cos(J(`obsYear',1,latLongYear[`i',2])-latLongYear[.,2])) );
			mata: zeroWeight = selectindex(weight :< 0);
			mata: weight[zeroWeight,1] = J(rows(zeroWeight),1,0);
			* acos of same location sometimes comes back at missing because of precicion issues, weight should be 1;
			mata: weight[`i',1] = 1;
			
			* step b: construct X'e'eX for given observation;
			mata: XeeXh = (J(1,`obsYear',Xyear[`i',.]') :* J(`K', 1, eYear[`i',1]*eYear' :* weight')) * Xyear;
			mata: XeeX = XeeX + (XeeXh + XeeXh')/2;
			mata: mata drop weight zeroWeight XeeXh;
		};
		mata: mata drop eYear Xyear obsYear latLongYear;
		display("Finished year `y'");
	};
	
	* normalize for degrees of freedom;
	mata: XeeX = XeeX / `df';

	* now get corrected variance-covariance matrix;
	mata: invXX = invsym(X'*X) * `df';
	mata: V = invXX * XeeX * invXX / `df';

	* reload the adjusted variance-covariance matrix;
	mata: st_matrix("V",V);
	mata: mata drop e X latLong year XeeX invXX V;

	* restore old data frame;
	restore;

	* return revised variance-covariance matrix;
    ereturn repost V = V;

end;






