# delimit ;
version 16.1;
clear all;
set more off;
set tracedepth 2;
* set trace on;


* ------------------------------------------------------------------------------------------;
* call subscripts for analysis;
* ------------------------------------------------------------------------------------------;
* best piecewise linear fit;
do 1_analysis/bestLinearModel;
bestLinearModel, crop(corn) lowerBound(10) monthMin(3) monthMax(8);
bestLinearModel, crop(soybeans) lowerBound(10) monthMin(3) monthMax(8);
bestLinearModel, crop(cotton) lowerBound(15) monthMin(4) monthMax(10);
*/;


* ------------------------------------------------------------------------------------------;
/* panel regression for Figure 1;
do 1_analysis/panelRegression;
do 1_analysis/conleyPanel.do;
panelRegression, crop(corn) pieceLinLB(10) pieceLinBreak(29) tempLB(-3) tempUB(36) monthMin(3) monthMax(8) conleyCutoff(1000);
panelRegression, crop(soybeans) pieceLinLB(10) pieceLinBreak(30) tempLB(-3) tempUB(36) monthMin(3) monthMax(8) conleyCutoff(1000);
panelRegression, crop(cotton) pieceLinLB(15) pieceLinBreak(31) tempLB(6) tempUB(39) monthMin(4) monthMax(10) conleyCutoff(1000);
*/;


* ------------------------------------------------------------------------------------------;
/* climate impact for Figure 2;
do 1_analysis/climateImpactUniform.do;
climateImpactUniform, crop(corn) pieceLinLB(10) pieceLinBreak(29) tempLB(-3) tempUB(36) monthMin(3) monthMax(8);
climateImpactUniform, crop(soybeans) pieceLinLB(10) pieceLinBreak(30) tempLB(-3) tempUB(36) monthMin(3) monthMax(8);
climateImpactUniform, crop(cotton) pieceLinLB(15) pieceLinBreak(31) tempLB(6) tempUB(39) monthMin(4) monthMax(10);
*/;


* ------------------------------------------------------------------------------------------;
/* out of sample predictions for Figure 3;
do 1_analysis/outOfSamplePrediction.do;
outOfSamplePrediction, crop(corn) pieceLinLB(10) pieceLinBreak(29) tempLB(-3) tempUB(36) monthMin(3) monthMax(8);
outOfSamplePrediction, crop(soybeans) pieceLinLB(10) pieceLinBreak(30) tempLB(-3) tempUB(36) monthMin(3) monthMax(8);
outOfSamplePrediction, crop(cotton) pieceLinLB(15) pieceLinBreak(31) tempLB(6) tempUB(39) monthMin(4) monthMax(10);
*/;


