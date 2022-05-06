# delimit ;
version 16.1;
clear all;
set more off;
set tracedepth 2;
* set trace on;


* ------------------------------------------------------------------------------------------;
* call subscripts for figures;
* ------------------------------------------------------------------------------------------;
/* best piecewise linear fit;
do 2_figures/r2_piecewiseLinear;
r2_piecewiseLinear;
*/;


* ------------------------------------------------------------------------------------------;
/* Figure 1 in paper: panel regression;
do 2_figures/Figure1a;
do 2_figures/Figure1b;
do 2_figures/Figure1c;
*/;


* ------------------------------------------------------------------------------------------;
/* Figure 2 in paper: climate impact;
do 2_figures/Figure2a;
do 2_figures/Figure2b;
*/;


* ------------------------------------------------------------------------------------------;
/* Figure 3 in paper: out-of-sample prediction error;
do 2_figures/Figure3;
*/;



