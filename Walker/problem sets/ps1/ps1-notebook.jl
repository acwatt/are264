### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9bf60494-b58e-495e-adc0-8c9d0585880b
# ╠═╡ show_logs = false
begin
	using ZipFile  # unzip compressed .zip folders
	using StatFiles  # read Stata files
	using DataFrames, Chain, DataFramesMeta
	using CSV
	# using Econometrics
	using CovarianceMatrices
	using CategoricalArrays
	using Statistics
	using StatsModels
	using StatsBase:sample
	using CUDA  # must be run before FixedEffectModels
	using FixedEffectModels
	using RegressionTables
	using Dates
	using Plots
	using GLM
	using Latexify
	using LaTeXStrings
	using RDatasets
	using StatsPlots
	using LinearAlgebra
	using Distributions
	using TimeSeries
	using Loess: loess, predict
	using PlutoUI
	Latexify.set_default(fmt = "%.3f")
end

# ╔═╡ 220c7746-b558-4b42-9e47-e227730a1aad
html"<span style='font-size:4em;'>ARE 264-b Problem Set 1</span>"

# ╔═╡ 8a325a10-4c17-43a1-a75f-9bc68f72fbdb
md"# Problem 1"

# ╔═╡ 00edd77c-be64-44fc-b510-28234152ca56
md"See the code appendix at the end for function definitions."

# ╔═╡ ac7d2e34-53a5-4c0e-b413-c1129378d91a
md"## 1.1 Temperature Aggregation"

# ╔═╡ 22681583-0918-49da-87b8-3a86c027cf81
md"### 1.1.1 Construct 4 temperature response variables
See next section."

# ╔═╡ b0a24d4d-016f-45a4-9124-90d463d1db93
md"### 1.1.2 Aggregate to year-county
Construction and aggregation are both handled in the `create_temperature_vars()` function."

# ╔═╡ 28d5ff08-5fbd-45c8-9816-dacbe9d20493
html"<br><br><br><br><br><br><br><br>"

# ╔═╡ f9fbde9e-b9fd-4407-a63d-992a28408245
md"## 1.2 US Climate Impacts: County-Year Damages"

# ╔═╡ 65cebd89-47ba-4773-9b7f-98ce0ef42d4b
md"### 1.2.1 log-transformed farm employment vs binned temperature"

# ╔═╡ c94aac41-7071-4cbd-acd6-9f1cd53176cd
md"
> Explore the relationship between log tranformed emp_farm and the vector of binned temperature controls. Include additional controls for county FE, year FE. Provide an interpretation of the coefficient of the 32+ bin.
"

# ╔═╡ 225f183d-aa3b-46ce-8edf-4336bba8da30
md"Below, I produce a plot of the coefficients on the temperature bin variables and their 95% confidence intervals, where the points are plotted at the midpoint of their respective bins. I used the 20-24°C bin as the reference bin in the regression, so all coefficients represent percent of annual Farm Employment in the county changed by moving a single day's average temperature from the 20-24°C bin to each respective bin. For example, the second point from the left should be interpreted as: *On average, if we move a single day that has a daily average of 22°C to a daily average of 2°C, we expect to see a decrease in that county's farm employment by roughly 0.1 percentage points, all other days held constant*."

# ╔═╡ 7214407b-7850-4a69-8798-59c463394cf2
md"Note that the left-most and right-most bin represent *below 0°C* and *above 32°C* average temperature, respectively. Thus the right most point in the plot represents: *% change in farm employees for the average county for each additional day above 32°C that would have been between 20°C and 24°C*"

# ╔═╡ 9cf39312-272b-4b73-867c-faf98da039f7
md"### 1.2.2 per capita farm property income vs cubic spline temperature"

# ╔═╡ fe21bed1-dcd3-47fa-997c-9827db649e3f
md"
> Explore the relationship between log tranformed `inc_farm_prop_income`/ population (i.e. log(per capita farm prop income) using the restricted cubic spline. Include additional controls for county FE and year FE. Plot the predicted marginal effects with associated confidence intervals, and compare to the binned temperature response function above.
"

# ╔═╡ 84fbd2e7-82cd-45ea-bb89-d75a621f79ca
md"
Let's compare the shift from 20°C days to cooler days between the plot above and the spline plot below. In the previous plot, we saw that farm employment decreases when shifting to cooler days, but in the below plot, the farm's income (normalized to the county population) seems to be fairly steady when moving moderate days to cooler days. Perhaps this is because it takes fewer hands to harvest when the season is cooler because the product ripens slower. 

In the plot above, farm employment incrases when shifting 20°C days to hotter days, and in the plot below, the farm's income seems to be signficantly decreasing when more moderate days are made hot. Similarly, perhaps hotter seasons require more intense (faster) harvesting of the produce. So more people could be employed in the county, but the farm may not make as much due to more product over-ripening.

The confidence intervals were created using the delta method. Bootstrapping was also explored -- clustered bootstrapping (at the county level,w ith 10,000 iterations) resulted in confidence intervals that were nearly identical to the delta method confidence intervals. The bootstrapping results are omitted here for brevity."

# ╔═╡ aa60b5ef-3aeb-42ff-b679-e00b28205267


# ╔═╡ 6a978c31-51c4-47a8-b626-14ca00de0370
md"### 1.2.3 Treatment heterogeneity"

# ╔═╡ 952ba66d-3264-405e-ad52-9f0a951e244c
md"
> Use the binned temperature estimator to design a test for whether we observe treatment effect heterogeneity. One possibility for this test is to interact the temperature bins with the average # of days in a county for which the temperature falls into each respective bin (i.e. are counties that experience more 32+ days more/less responsive to temperature).

To estimate treatment heterogeniety, I first calculated the 10-year moving average (backward looking) for each year-county-temperature bin. This is functionally a medium-term historical average of days of the year with average temperatures in each bin. So the first 10 years of data are now unused for these plots. For each year-county-bin, I then calculated the difference of the bin from the 10-year average (the deviance from the average). I then regressed both Farm Employment and Farm Income per capita on the new deviance bins. "

# ╔═╡ 9a4d0f65-7c09-4c92-99fb-ff31b0453054
md"To compare, I again plot the basic temperature bins' effects, but for both Farm Employment and Farm Income"

# ╔═╡ 68e3f5e2-382a-431e-a19e-3547e70c3a94


# ╔═╡ b460634b-72c1-4f56-ad00-2a550f1c660f


# ╔═╡ 288e1ee5-809f-4108-80d1-15b5c8e67452


# ╔═╡ 4fe46d77-791f-4e80-b715-2ccb12de37f5


# ╔═╡ 714b920e-a30a-4930-bd05-fb070dc4f641
md"# Problem 2"

# ╔═╡ 8255f716-16ad-44c2-9c87-1710f95ff1f6
md"## 2.2 Hedonic Air Quality Analysis"

# ╔═╡ 478897a9-13e7-453a-ad0a-4b4031e7932d
# ╠═╡ disabled = true
#=╠═╡
md"Load pollution data"
  ╠═╡ =#

# ╔═╡ a613f88c-ea92-4e66-8a62-6f9012437cb9
md"### 2.2.1  Housing prices vs air pollution"

# ╔═╡ 71b0a94f-67ae-4240-ba7a-5e9ebd33af45
md"
> Estimate the relationship between changes in air pollution and housing prices, both not adjusting and adjusting for other housing price determinants (use pop7080 as weights). What do your estimates imply and do they make sense? Describe the potential omitted variables biases. What is the likely relationship between economic shocks and pollution and housing price changes? Using the observable measures of economic shocks (dincome, dunemp, dmnfcg), provide evidence on this."

# ╔═╡ cac486e4-3bd6-45ab-b6d0-24dbafc0da71
md"
Both regressions (1) and (2) seem to indicate that increased pollution increases
housing values! This doesn't makes sense since air pollution
should be a disamenity. This result could be caused by OVB,
if the omitted variable has correlation with housing price in
the same direction as the correlation with pollution (ie,
if x2 is the OV, cov(x2, pollution) and cov(x2, houseprice) are both the same sign).
For example, if there is a recession, and the economy slows, 
this generally decreases the amount of pollution being produced
by firms supplying goods to the economy. So an indicator of positive
economic activity would be positively correlated with pollution. At the same
time, a recession would likely cause the average price of housing to decrease.

So an indicator of positive
economic activity would be positively correlated with housing price.
Then bias from an omitted economic variable that is positively correlated with
both pollution and housing prices would cause an upward bias in our estimate
of the partial derivative of house price changes w.r.t. pollution changes.
If the bias is large enough, it could flip the sign of the coefficient,
and since our coeffience on pollution is nearly statistically indistinguishable from 0
after controlling for other observed housing determinants, upward OVB could
be likely here.

There was a large recession in the 70s, so variables describing the impact of 
the recession on the supply and demand for housing seem like they
would be important to control for. For example, many manufacturing jobs were lost
during the recession (more relative to other sectors)
which affects the housing market. We can test three of our economic
indicators to see how correlated with pollution they are in order to see if
we should be concerned about OVB.

We could further
control for unemployment more generally and average county
income changes to ensure we are comparing willingness to pay across
similar counties.  After controlling for three economic 
indicators, regression (3) above shows that the hedonic estimate of willingness to pay for
changes in pollution become statistically indistinguishable from zero.
"

# ╔═╡ af154adc-56f4-45d8-9c5c-7647fc7ed12e
md"
Above is the regression of changes in pollution levels on changes in various eocnomic indicators. Economic shocks can change pollution because pollution
is generated by economic activity (i.e., if there is 
a recession, there will be less pollution) and economic shocks
can change what people are able/willing to pay for a house 
(ie, recession -> house prices decrease), which creates
endogeniety (house price is correlated with the error term).
"

# ╔═╡ b27cc951-9bd5-45f6-b30f-51010f2c9537


# ╔═╡ fb3ad8e9-08f9-4ffa-891b-e6bd926897a7


# ╔═╡ bd0603ce-3a15-4428-af4d-c2ee2127e537
md"### 2.2.2  Attainment status as IV for pollution changes"

# ╔═╡ d084de41-4f40-4e51-8dca-54ec0608e72b
md"""
> Suppose that federal EPA pollution regulation is
> a potential instrumental variable for pollution
> changes during the 1970s. What are the
> assumptions required for mid-decade regulatory
> status (tsp7576) to be a valid instrument for
> pollution changes when the outcome of interest is
> housing price changes?

*Assumptions*:

1. Relevance: tsp7576 needs to be correlated with pollution changes.

2. Exogeneity: tsp7576 only affects housing prices through it's effect on changes in pollution.
"""

# ╔═╡ 1185d096-15c8-42e1-bbff-8047722dbabc
md"""
> Provide evidence on the
> relationship between the regulatory status
> indicator and the observable economic shock
> measures. Interpret your findings.
"""

# ╔═╡ 0df4996b-c5e2-4ae1-9909-09b442cc9ff2
md"
The above regression of the regulator status indicator on our measures of economic shock shows that there is vary little statistical relationship between the regulatory status indicator and economic shocks. This supports our assumption of exogeneity of the regulatory status as an IV for changes in pollution.
"

# ╔═╡ adbf35de-93e8-492d-a89d-c61697da3391


# ╔═╡ bc11f842-1b19-4250-a969-81b85a9f7dc4


# ╔═╡ d9face0e-651d-4256-8aab-3c257a18619d
md"### 2.2.3 First stage, reduced form"

# ╔═╡ 6bbe502f-3750-4321-87c1-a01d67e228a1
md"""
> Document the first-stage relationship between
> regulation and air pollution changes and the
> reduced form relationship between regulation and
> housing price changes, both not adjusting and
> adjusting for other covariates. Interpret your
> findings. How does two-stage least squares use
> these two equations? Now estimate the effect of
> air quality changes on housing price changes
> using two-stage least squares and the tsp7576
> indicator as an instrument (not conditioning and
> conditioning on other observables). Interpret the
> results.
"""

# ╔═╡ 95412c61-3a1b-40bf-9b7c-ddd54ce17ba7
md"*First-stage* (relevance) relationship between regulation and air pollution changes; and **Reduced form* relationship between regulation and housing price changes. (SE in parentheses)"

# ╔═╡ 0cc1daaa-9358-4074-a576-3f4263c18e95
md"**How does two-stage least squares use these two equations?**

The 2SLS regression coefficient is the ratio of the reduced form relationship to the first-stage relationship. Formally, if the first stage relationship is given by:

$x = \alpha z + u$

and the reduced form relationship is given by

$y = \gamma z + \varepsilon$

then the IV (2SLS) coefficient $\beta$ that we want to estimate from 

$y = \beta x + e$

will be 

$\beta = \frac{\gamma}{\alpha}$

So, without using statistical controls, our IV estimate of the effect of pollution on house prices can be calculated by dividing the coefficient from regression (3) above by the coefficient from regression (1). With statistical controls, we would divide the coefficient from (4) by that from (2).

I run the 2SLS model below, and then test if the ratio of the above coefficients are equal to the 2SLS coeficients below. They do indeed match, both with and without controls. However, the 2SLS model below corrects the estimates the standard errors for the ratio.

Looking at regression (2) in the table below, the coefficient of -0.007 indicates that a 1-unit decrease in *the change in mean TSPs* from 1970 to 1980 results in an average
housing price change from 1970 to 1980 *being 0.7 percentage points larger* (more positive).

To make the interpretation more concrete, let's look at an example county. First note that the dependent variable (`dlhouse`) is the natural log of the ratio of housing prices in the county (1970 / 1980). So, for a given value of `dlhouse`,

$dlhouse = log\left( \frac{hprice_{1980}}{hprice_{1970}}\right)$

a percentage increase would be

$\text{\% increase from 1970 to 1980} = 100\cdot(e^{dlhouse}-1)$

"

# ╔═╡ 6c00497f-91bf-4122-a366-e2bb119bbcdc
md"*2SLS Regression using an indicator for regulation as an IV for pollution changes* (SE in parentheses)"

# ╔═╡ 76650858-6af6-4060-9487-4d6424dc7cb7
md"*Testing the 2SLS regression coefficients (regressions `223e` and `223f`) for equality with the ratio of coefficients from the reduced-form and first-stage regressions:*"

# ╔═╡ c6e1abd9-0490-40cf-a8db-9dfe3bd8249a


# ╔═╡ 26a2778a-fe46-4772-8e78-e09b63ae4152


# ╔═╡ f89066c8-39d0-4e85-b924-3e493aabda88
md"### 2.2.4 Mean of TSPs as IV for pollution changes"

# ╔═╡ 22802ab4-ed98-40fe-9c20-4207843c2cb2
md"
> In principle, the regulation indicator variable
> should be a discrete function of pollution levels
> in 1974. Specifically, the EPA is supposed to
> regulate those counties in 1975 who had either an
> annual geometric mean of TSPs above 75 units
> (mg/m3) or a 2nd highest daily concentration
> above 260 units in 1974. Based on this notion,
> redo part (3) using mtspgm74 as an instrument for
> pollution changes. Interpret your findings."

# ╔═╡ 4aaa6d2c-1ea7-400b-a8e1-9e0993c9a5f7
md"*First stage and reduced form regressions* (SE in parentheses)"

# ╔═╡ b989539a-09f4-4ecc-8f24-d3efb5962035
md"*2SLS Regression using an indicator for above 75 pollution units as an IV for pollution changes* (SE in parentheses)"

# ╔═╡ ae43ae55-ad11-43a9-9f25-9123ee28a885


# ╔═╡ 90903d07-605a-4b51-91e7-b2f5bac56c29


# ╔═╡ 2f203087-0dab-49a9-9115-e3b8477d3089
md"### 2.2.5 Discontinuity in treatment assignment"

# ╔═╡ 78ac5bab-75c5-4ccf-aa63-1681bce62b9c
md"
> Describe how one could use this discontinuity in treatment assignment to derive an alternative estimate of the capitalization of pollution changes. Under what conditions will this estimate be valid? Based on this concept, estimate the nonparametric bivariate relationships between pollution changes and 1974 TSPs levels and housing price changes and 1974 TSPs levels. Do this using the lowess STATA command or a similar command in R (experiment with relatively small bandwidths between 2-4). Plot the two conditional mean functions on either side of the discontinutiy (changes in air quality for discrete 1974 values on either side of threshold). Interpret your findings and relate them to the results in (4).
"

# ╔═╡ 573af536-2a9e-4f7a-8a68-0b5e3157bddf
md"
I can imagine using either regression discontinuity (RD) or differences-in-differences (DiD) around the 1974 TSP threshold to analyze the capitalization of pollution changes. Because RD takes a relatively small group of observations around the threshold, it depends on having enough data to estimate the effects. Because the regulation is not only dependend on the 1974 threshold for the mean TSPs (it also depends on the 2nd highest daily concentration), this would be a fuzzy RD. DiD on the other hand could use the entire dataset above and below the threshold and would almost surely have more smaller standard errors on the treatment effect estimates of the regulation.
"

# ╔═╡ 2145f8ab-2e3a-4b88-803b-f73ffdd0b4e6
md"*Assumptions of the RD:*
- No manipulation: there shouldn't be signs of units manipulating their values very close to the cutoff to avoid being regulated.
- Units (counties) just above and just below the threshold should be very similar along other observed and unobserved characteristics that are important in determining the outcomes.
"

# ╔═╡ 85e03d2c-af57-4feb-9107-dedd0d59ecd0
Bandwidth = @bind bandwidth Slider(0.02:0.005:1.5, default=0.1)

# ╔═╡ f7aaaee5-262c-40cc-8ff6-793a904253fa
PredictionPoints = @bind npoints Select(200:100:2000, default=1000)

# ╔═╡ 3ab934a9-4dc0-473b-be1d-f290761dddd8
PolynomialDegree = @bind degrees Select([2, 3], default=3)

# ╔═╡ 4df0e4c6-55f4-401d-b18e-c4cba40d0de4
begin
	function nomissing(df, cols)
	    df2 = dropmissing(df, cols)
	    return (df2[!, c] for c in cols)
	end
	±(x, y) = x+y, x-y
	
	myloess(x::Vector{Float32}, y::Vector{Float32}; bandwidth) = loess(Float64.(x), Float64.(y); span=bandwidth, degree=degrees)
	myloess(x::Vector{Float32}, y::Vector{Float64}; bandwidth) = loess(Float64.(x), y; span=bandwidth, degree=degrees)
	myloess(x::Vector{Float64}, y::Vector{Float32}; bandwidth) = loess(x, Float64.(y); span=bandwidth, degree=degrees)
	myloess(x::Vector{Float64}, y::Vector{Float64}; bandwidth) = loess(x, y; span=bandwidth, degree=degrees)


	function evenly_spaced(a, b, n)
	    h = (b-a)/(n-1)
	    Float64.(collect(a:h:b))
	end

	function get_loess_prediction(df, xcol, ycol; bandwidth, npoints=200)
		x = df[!, xcol]; y = df[!, ycol];
		model = myloess(x, y; bandwidth=bandwidth)
		xs = evenly_spaced(extrema(x)..., npoints)
		ys = predict(model, xs)
		return xs, ys
	end
	end;

# ╔═╡ 2d5817b1-a1cd-4882-9c3d-ac53eb9362f2
PlotWidth=@bind plotwidth Slider(0.2:0.2:100, default=10)

# ╔═╡ 0c8a7d25-7233-45d5-b3c1-f46bd423cb94
md"But when we zoom out, the nonparametric estimation further away from the threshold is even more noisey. So how do we know the difference in the estimates near the threshold aren't just due to sampling variation? One way to examine the variance due to sampling this is to use permutation inference: hold the bandwidth constant, and resample the data many times, plotting the mean of the LOESS estimate and the 2.5``^{th}`` and 97.5``^{th}`` percentiles of the prediction at each point."

# ╔═╡ eb368d82-03b9-47a0-839f-b90f8eb10e09
md"When we look below at changes in house prices based on the 1974 TSP levels, and using a wider bandwidth, it appears that the house price changes were smaller just below the threshold."

# ╔═╡ b59f26c1-9792-4958-86ae-3712f19854fa
BandwidthHousePrice = @bind bwdlhouse Slider(0.02:0.005:1.5, default=0.325)

# ╔═╡ 9036f9b8-e157-495b-b2a5-b583935170f0
PlotWidthHousePrice = @bind pwdlhouse Slider(0.2:0.2:100, default=3)

# ╔═╡ 44aeda75-7419-40cb-9234-f01945e4d9b9
md"Same bandwidth, just zoomed out to see more of the data."

# ╔═╡ f0a30424-670b-4200-98bd-e9f4ecd38eb2
md"In part (4), we saw that less improvement in pollution levels, as instrumented by being below this threshold, results in smaller increases in housing prices. If we were to take the plot above at face value, we can see that counties just below the threshold experienced a smaller increase in housing prices compared to counties just above the threshold. From the first two plots above of the pollution changes, we saw that counties just below the threshold (presumably less regulated) saw less of a pollution improvement (the decrease in pollution was smaller in magnitude, on average, just below the threshold).

So, here we see that less improvement in pollution levels (just below the threshold) is correlated with a smaller increase in housing prices (though this result is highly dependent on the choice of bandwidth)."

# ╔═╡ 834fe856-c053-4076-901c-754eae2c8846


# ╔═╡ d4dee95d-f991-44ea-9c52-95d51f8c1ab2
md"### 2.2.6 The Smoothness Criteria of RD"

# ╔═╡ bd95a3c4-553e-418f-9831-dda93068a9be
md"
> Now use linear regression to estimate the predicted change in housing prices based on the control variables (excluding the change in TSPs). This provides a single-index measure of the housing price changes predicted to occur due to other variables changing. Estimate the nonparametric bivariate relation between this index and 1974 TSPs levels and plot it against the smoothed housing price changes from (5). Explain how this provides an indirect test of the smoothness condition required by the regression discontinuity design and interpret your findings. "

# ╔═╡ 202a5a35-7633-44cc-887e-22c73d3019d0
md"The second assumption about RD from part (5) describes covariate balance between observations just above and just below the cutoff.
This is roughly equivalent to a smoothness condition across the cutoff: the predicted change in housing prices due to all other variables
should be smooth across the cutoff, so that we can attribute any observed discontinuity in the change in housing prices
across the cutoff to the treatment (regulation)."

# ╔═╡ 4b46b382-5349-4bc7-bb03-229c303f37be
md"
Ideally, to visually test the smoothness condition of the RD design, we would
want to see a noticable discontinuity at the threshold for the kernel regression
of observed changes in housing price on TSP levels, while also seeing a smooth
transition across the threshold for kernel regression of predicted
changes in housing price on TSP levels. Since the predicted house price changes
are based on control only (and not pollution changes), if those changes
show a smooth relationship across the TSP threshold, that is evidence
of good covariate balance across the threshold (ie, it's less likely that
there is a large difference in the covariates across the threshold that could
also create a discontinuity in house price changes). 

In the plots above, we see some evidence of this, but it is not conclusive. With
this bandwidth, the blue kernel regression (observed data) shows what could be a noticable
discontinuity across the threshold. Due to the lack of more data very close
to the threshold, it is hard to tell what could happen in the middle. But
the green kernel regression (predicted data) looks a little more like it could 
meet at the threshold. However, even the green lines show some difference on
either side of the threshold. But this might provide some indirect evidence that
there is some discontinuity happening at the threshold that is not accounted for
by the control variables.
"

# ╔═╡ 98857139-6cd1-4e39-9e1e-8a35f2629e28


# ╔═╡ 8dc1d528-a9a1-4f99-8b83-a864e5beac6b


# ╔═╡ bdf52a46-4f47-4bc5-9e07-a20fb7b76dbf
md"### 2.2.7 Regulated vs Unregulated"

# ╔═╡ ce91dbb9-70ec-43e5-978a-68cd7817a6fd
md"
> A number of counties with annual geometric mean TSPs below 75 units in 1974 were regulated due to having as few as 2 bad days. This implies that we can compare regulated and unregulated counties with identical average TSPs levels in the regulation selection year (below 75 units). For those counties with 1974 mean TSPs between 50 and 75 units, estimate the bivariate relation between TSPs changes and 1974 TSPs levels separately for regulated and unregulated counties and plot the results. Now do the same for the relation between log-housing price changes and 1974 TSPs levels. Interpret your findings. Since there are fewer observations, you may need to use bigger bandwidths than those in part (5) (e.g., bandwidths between 6-9). "

# ╔═╡ 948f96b3-6f50-4687-8d4f-db309c8e29a2
BandwidthPollChange = @bind bwdgtsp Slider(0.04:0.005:1.0, default=0.325)

# ╔═╡ 16d94c66-4a4e-443a-9f0b-ea2af6e4555f
BandwidthHousePrice2 = @bind bwdlhouse2 Slider(0.04:0.005:1.0, default=0.325)

# ╔═╡ e5981e73-8029-4e7b-aedb-76b1080d882a
md"With a fairly large bandwidth (compared to earlier plots):
For changes in pollution, we can see that as we get closer to the 75-unit
threshold, the regulated counties show pollution improvements 
(larger negative pollution changes), while the unregulated counties
seem to have fairly uniform average pollution improvements.

For changes in house prices, we also see that house prices for 
regulated counties that are higher in average TSPs in 1974
(closer to the 75-unit threshold) saw larger house price increases
than unregulated counties around the same 1974 TSP levels. This provides
more evidence that regulation, conditional on average TSPs, seems
to have increased house prices (thus, less pollution is valuable
to house buyers)."

# ╔═╡ 4321d60b-fc24-4acf-9184-f7a68e7909d8


# ╔═╡ 9d9db365-0799-4699-8c75-c16b1820d322


# ╔═╡ 40b4be47-8c52-466a-8130-cfc9b256452d
md"### 2.2.8 2SLS and the ATE"

# ╔═╡ da4a328f-ef3f-44d4-90cb-b22b1519fad0
md"
> Under what assumptions will two-stage least squares identify the average treatment effect (ATE)? If ATE is not identified, describe what may be identifiable with two-stage least squares estimation. Under what conditions is this effect identified? Give some intuition on what this effect may represent when one uses EPA regulation as an instrument variable.
"

# ╔═╡ 57754134-0b34-41b5-9366-fe7e469ee1af
md"

"

# ╔═╡ 872ffa35-9f86-461c-bf77-29f4a8ad5e26


# ╔═╡ b5cf0297-5489-488f-a75f-586adca9f222


# ╔═╡ 730c1994-cbee-4ec9-95a8-51e0923ef4ed
md"### 2.2.9 Summary of TSP regulation effects on House Price Changes"

# ╔═╡ 26d3b991-cc5c-455e-b38d-5925069cd713
md"
> Provide a concise synthesis/summary of your results. Discuss the credibility of the various research designs underlying the results.
"

# ╔═╡ 5f31d109-b75e-4cc5-9cae-a06493009221
md"

"

# ╔═╡ 96ba58aa-8028-4db7-bd76-ef3fc7f24637


# ╔═╡ 01ce9f00-209a-4f4b-bd63-d7b02767f5a0


# ╔═╡ 8b3041fe-3b6d-4ef2-9223-5f81a6366019
md"# Code Appendix: Functions"

# ╔═╡ 7d7b39f2-1db6-4c0b-a41b-710030b7537d
md"## File Functions"

# ╔═╡ 8893d6e2-2e3a-440a-ba2b-e1b60a5b9d73
"""Extracts file_name from zip folder zip_name in directory root_dir."""
function extract_file_from_zip(root_dir, zip_name, file_name)
    zip_path = joinpath(root_dir, zip_name)
    save_path = joinpath(root_dir, file_name)
    # If the file is already extracted, done.
    if isfile(save_path)
        return
    end

    println("Extracting $file_name from $root_dir")
    # Open zip file
    zarchive = ZipFile.Reader(zip_path)
    # Find file_name in zip archive
    files = [f.name for f in zarchive.files]
    idx = findfirst(x -> x == file_name, files)
    # Save file_name to disk
    write(save_path, read(zarchive.files[idx]))
end

# ╔═╡ 892333c3-0c5e-4d58-bc0d-84502afe4c96
"""Return dataframe, after unzipping and saving the file."""
function df_from_zip(root_dir, zip_name, filename)
    extract_file_from_zip(root_dir, zip_name, filename)
    println("Loading $filename from file.")
    return DataFrame(load(joinpath(root_dir, filename)))
end

# ╔═╡ af4f8032-c318-4035-913f-351c8e0130f7
"""Return a dataframe, with types Float64 and no missing for given columns."""
function convert_float32(df, cols)
    df1 = dropmissing(df, cols=cols; disallowmissing=true)
    df1 = transform!(df1, cols .=> ByRow(Float64), renamecols=false)
    return df1
end

# ╔═╡ 93eae3b1-187b-4026-81e7-a587d1477af9
md"## Temperature Functions"

# ╔═╡ a72bfb4b-0a93-48c5-81f7-6156c88c4ef4
"""Degree day calculation from Snyder 1985, using integral of sine."""
function degree_day_calc(MAX, MIN, THR)
    M = (MAX + MIN) / 2  # average temp
    W = (MAX - MIN) / 2  # half temp range
    θ = asin((THR - M) / W)  # radian between -π/2 and π/2 where sine crosses THR
    DD = ((M - THR)*(π / 2 - θ) + W * cos(θ)) / π
    return DD
end

# ╔═╡ cdc1dd7f-683c-40c8-9410-3fcdb4d74322
"""Return degree days above threshold temperature for single day based on max and min temperatures."""
function degree_day_single_day(Tmax, Tmin, Threshold)
    if Threshold ≥ Tmax
        return 0
    elseif Threshold ≤ Tmin
        return (Tmax + Tmin) / 2 - Threshold
    else
        return degree_day_calc(Tmax, Tmin, Threshold)
    end

end

# ╔═╡ b0bf6cbc-e1e4-42f1-8d0c-055f5a249faa
"""Return dataframe with degree day columns, degree days over each threshold in t."""
function add_degree_days(df::DataFrame, t::Vector)
    println("Adding degree days above thresholds in t.")
    for k ∈ 1:length(t)
        df[!, "dday$(t[k])C"] = degree_day_single_day.(df[!, :tMax], df[!, :tMin], t[k])
    end
    return df
end

# ╔═╡ 79bf9c02-b2a1-40a2-b87e-d436d587fe13
"""Return dataframe with binned temperature variables columns (indicator for tAvg in bin)."""
function add_bin_indicator_days(df::DataFrame, edges::Vector)
    println("Adding inidicator for average daily temperature in each bin")
    for k ∈ 1:length(edges)
        if k == 1
            # tAvg below lowest threshold
            df[!, "tempB$(edges[k])"] = (df[!, :tAvg] .< edges[k])
        else
            #  tAvg between this threshold and previous threshold
            df[!, "temp$(edges[k-1])to$(edges[k])"] = (edges[k-1] .≤ df[!, :tAvg] .< edges[k])

            if k == length(edges)
                # tAvg above threshold
                df[!, "tempA$(edges[k])"] = (df[!, :tAvg] .≥ edges[k])
            end
        end
    end
    return df
end

# ╔═╡ 790ce5d9-6119-4c5c-b057-6e0a52ecdaf2
"""Return temperature threshold after converting to -1,1 interval (for sine function)"""
function convert_temp_threshold(Tmax, Tmin, Thresh)
    proportion = (Thresh - Tmin) / (Tmax - Tmin)
    return (2 * proportion - 1)
end

# ╔═╡ 0e4a2b18-9a0f-4acd-a4e8-920a1bac0b87
"""Return portion (0 to 1) of day's temperature above temperature Threshold,
    given max and min temperature, 
    based on sinusoidal curve to model diurnal temperature cycle.
"""
function fraction_above_threshold(Tmax, Tmin, Threshold)
    α = convert_temp_threshold(Tmax, Tmin, Threshold)
    if α ≥ 1
        return 0
    elseif α ≤ -1
        return 1
    else
        return (π - 2 * asin(α)) / (2π)
    end
end

# ╔═╡ 24c043c5-dc9b-4a7a-a815-bdd4bb1997bb
"""Return portion (0 to 1) of day's temperature between temperature Thresholds,
    given max and min temperature, 
    based on sinusoidal curve to model diurnal temperature cycle.
"""
function fraction_between_thresholds(Tmax, Tmin, Threshold_lower, Threshold_upper)
    # Fraction between = (fraction above lower) - (fraction above upper)
    lower = fraction_above_threshold(Tmax, Tmin, Threshold_lower)
    upper = fraction_above_threshold(Tmax, Tmin, Threshold_upper)
    return lower - upper
end

# ╔═╡ ac20872f-f7e9-480f-8155-b1da7c616909
"""Return dataframe with binned temperature variables columns (portion of day in bin)."""
function add_bin_portion_days(df::DataFrame, bins::Vector)
    println("Adding portion of days between bin edges in bins.")
    for k ∈ 1:length(bins)
        if k == 1
            # Portion of day below lowest threshold
            df[!, "tempPropB$(bins[k])"] = 1 .- fraction_above_threshold.(df[!, :tMax], df[!, :tMin], bins[k])
        else
            # Portion of day between this threshold and previous threshold
            df[!, "tempProp$(bins[k-1])to$(bins[k])"] = fraction_between_thresholds.(df[!, :tMax], df[!, :tMin], bins[k-1], bins[k])

            if k == length(bins)
                # Portion of day above threshold
                df[!, "tempPropA$(bins[k])"] = fraction_above_threshold.(df[!, :tMax], df[!, :tMin], bins[k])
            end
        end
    end
    return df
end

# ╔═╡ 07e6f171-e8bb-4d96-acc8-0d0694692141
"""Return vector of element-wise minimum between a vector and constant."""
min_vecs(vec, c) = minimum.(zip(vec, repeat([c], length(vec))))

# ╔═╡ 78d8e678-4148-451d-a2de-9d84986d522b
"""Return vector of element-wise maximum between a vector and constant."""
max_vecs(vec, c) = maximum.(zip(vec, repeat([c], length(vec))))

# ╔═╡ 779e117a-3e39-4e34-9e09-d6d10c8a4f22
"""(⋅)₊ function from Harrell 2001"""
m(a) = maximum([0, a])

# ╔═╡ b6218121-c5cb-4d53-bfab-69a10a8aa241
"""Return jth restricted cubic spline variable using the list of knots t.
    From Frank Harrell's Stats Textbook: Regression Modeling Strategies, 2001, sect. 2.4.4
"""
function cubic_spline_var_harrell(x, j, t)
    k = length(t)
    if (j < 1) | (j > k - 1)
        println("Out of Bounds (j=$j): j is < 1 or > k-1. Not a valid variable.")
    elseif j == 1
        return x
    else
        xⱼ = m.(x .- t[j-1]) .^ 3 .-
             m.(x .- t[k-1]) .^ 3 * (t[k] - t[j-1]) / (t[k] - t[k-1]) +
             m.(x .- t[k]) .^ 3 * (t[k-1] - t[j-1]) / (t[k] - t[k-1])
        # Divide by the scale factor from the stata formula (range of the knots)
        return xⱼ / (t[k] - t[1])^2
    end
end

# ╔═╡ b9a9e3eb-f707-4e25-87de-9d6d53a9db60
"""Add restricted cubic spline basis variables to dataframe df.
    Using continuous variable varname from dataframe and list of knots.
    After trying both the stata and Harrell formulas, neither matched the example file given. After comparing the formulas, I recognized the only difference was a scaling factor of 
        `julia 1/(knots[last] - knots[first])^2`
    so I added that scaling factor to the Harrell formula and it produces results that match the example file. So I probably just messed up the stata formula somehow. I assume the example file was created in stata with mksplines.
"""
function add_cubic_spline_vars(df, varname, knots; newbasename="splineC")
    for i in 1:(length(knots)-1)
        newname = "$newbasename$i"
        println("Making $newname")
        df[!, newname] = cubic_spline_var_harrell(df[!, varname], i, knots)
    end
    return df
end

# ╔═╡ fe94d889-551e-47ea-8363-6c6064a3df1b
"""Return ith restricted cubic spline variable using the list of knots k.
    From stata's mkspline documentation
"""
function cubic_spline_var_stata(V, i, k)
    n = length(k)
    if (i < 1) | (i > n - 1)
        println("Out of Bounds (i=$i, n=$n): i < 1 or i > n-1. Not a valid variable.")
    elseif i == 1
        return V
    else
        Vᵢ = (m.(V .- k[i-1]) .^ 3 .- (k[n] - k[n-1])^(-1) *
                                      (m.(V .- k[n-1]) .^ 3 * (k[n] - k[i-1]) .- (V .- k[n]) .^ 3 * (k[n-1] - k[i-1]))) ./
             (k[n] - k[1])^2
        return Vᵢ
    end
end

# ╔═╡ 0cd2f8f9-a72d-4600-89c5-8ee99f3009f7
"""Return ith linear piecewise variable using the list of knots k.
    From stata's mkspline documentation
"""
function piecewise_linear_var_stata(V, i, k)
    n = length(k) + 1
    if i == 1
        return min_vecs(V, k[i])
    elseif i > 1 && i < n
        return max_vecs(min_vecs(V, k[i]), k[i-1]) .- k[i-1]
    elseif i == n
        return max_vecs(V, k[n-1]) .- k[n-1]
    else
        println("Out of Bounds (i=$i, n=$n): i < 1 or i > n-1. Not a valid variable.")
    end
end

# ╔═╡ d18c8a79-d9c5-4efb-aeae-c2e58e07251c
"""Return dataframe with piecewise linear basis variable columns using given knots."""
function add_piecewise_linear_vars(df, varname, knots; newbasename="piece")
    println("Adding piecewise linear basis variable columns with breakpoints.")
    n = length(knots)+1
    for i in 1:n
        e = i == 1 ? "B$(knots[i])" :
            i == n ? "A$(knots[n-1])" :
            #= o.w =# "$(knots[i-1])to$(knots[i])"
        newname = "$newbasename$e"
        println("Making $newname")
        df[!, newname] = piecewise_linear_var_stata(df[!, varname], i, knots)
    end
    return df
end

# ╔═╡ bcfdbe23-5290-40fa-9965-a1f254c21bee
"""Return grid-day data aggregated to county-year level"""
function aggregate_df(df, tags)
    # List of columns to average over
    avg_list = [:tMin, :tMax, :tAvg, :prec]
    # List of columns to sum over the year (all cols with a substring from tags)
    sum_list = [n for n in names(df) if any(occursin.(tags, n))]
    # sum_list = [:splineC1, :splineC2, :splineC3, :splineC4]

    # Sum over all days in each year, for each grid point
    df[!, "date"] = Date(1960, 1, 1) + Day.(df.dateNum)
    df[!, "year"] = Year.(df.date)
    gd_gy = groupby(df, [:gridNumber, :year])
    df_gy = combine(gd_gy,
                    avg_list .=> mean .=> avg_list,
                    sum_list .=> sum .=> sum_list)
    # Take the average over all grid points in the county
    gd_y = groupby(df_gy, :year)
    df_y = combine(gd_y, valuecols(gd_y) .=> mean .=> valuecols(gd_y))

    return df_y
end

# ╔═╡ fec1fd19-b0e0-48c4-9ade-140c0c068379
md"## Analysis Functions"

# ╔═╡ de2d137d-a9cc-43b2-ba9f-caa0502c6210
md"### 1.2.2 Functions"

# ╔═╡ 45a1bad3-2d2c-409b-a742-8401ff15f166
"""Return DF of avg temp, cubic spline, and predicted outcome variables needed for regressions and plotting.
    
    @param: reg_output: output from a regression on 4 cubic spline variables
    
    - create an x-axis for plotting: define an array of avg temps from 0 to 40°C, 0.25 steps
    - define the spline knots/breakpoints
    - Create spline basis variables. After regressing, the absolute level of the outcome (y-axis)
        is not identified because the county & year fixed effects shift the intercept around.
        So these basis will be compared to the spline basis variables created for X°C.
    - Create spline basis variables for X°C
    - Create relative basis variables (base - X°C variables) which is equivalent to shifting
        the base prediction by the X°C prediction, but allows for adjusted standard errors
        in the delta method calculation.
"""
function plotdf_init(reg_output; X=20)
    df = DataFrame(tAvg=0:0.25:40)      # Initialize x-axis
    knots = [0 8 16 24 32]              # breakpoints in spline
    # Create base basis spline variables
    df = add_cubic_spline_vars(df, :tAvg, knots, newbasename="splinebaseC")
    # Create X°C basis spline variables
    df[!, "XC"] = repeat([X], nrow(df))
    df = add_cubic_spline_vars(df, "XC", knots, newbasename="XC")
    # Create relative basis spline variables (relative to X°C)
    for k in 1:4
        df[!, "splinerelativeC$k"] = df[!, "splinebaseC$k"] - df[!, "XC$k"]
    end
    # Construct the predicted spline outcome
    df[!, :yhat] = Matrix(df[!,r"splinerelativeC"]) * reg_output.coef
    # Return the x-axis, predicted outcome, splinerelCnd relative spline variables (for delta method SE estimation)
    return df[!, [["tAvg", "yhat"]; ["splinerelativeC$k" for k ∈ 1:4]]]
end

# ╔═╡ 17ec6db5-aa8c-4488-acfb-595c245f2c23
"""Return vector of diagonal elements of matrix"""
function get_diagonal(M::Matrix)
    n,_ = size(M)
    return [M[i,i] for i∈1:n]
end

# ╔═╡ 4ac46ab5-1615-4feb-bfb6-bc9cef26c353
"""Return DF with upper and lower bounds created from delta method.

    Delta method: SE =  √(∂f(β)/∂β ⋅ VCOV ⋅ ∂f(β)/∂β)
    Calculate the β-derivative of the function of parameters
    f(β) = Σxₖ⋅βₖ ⟹ ∂f∂βₖ = xₖ
    ∂fβ_∂β = Matrix(df_plot[!,r"splineC"])
    SEs are the diagonal of .√(∂fβ_∂β * VCOV * ∂fβ_∂β'))
    where VCOV is the variance-covariance matrix of parameter estimates
    from the regression results
"""
function deltamethod_CIbounds!(df, reg_output; column_stem="splinerelativeC", α=0.05)
    VCOV = reg_output.vcov
    ∂fβ_∂β = Matrix(df[!, Regex(column_stem)])
    n = reg_output.nobs
    df[!, :SE] = .√(get_diagonal(∂fβ_∂β * VCOV * ∂fβ_∂β'))
    # Use a t-distribution just in case the sample is small
    tcrit = quantile(TDist(reg_output.dof_residual), 1-α/2)
    # Use √(1 + 1/n) correction because both the mean (predicted value) and the variance are unknown estimated parameters
    width = tcrit * df[!, :SE] * √(1 + 1/n)
    df[!, :y_lb] = df[!, :yhat] - width
    df[!, :y_ub] = df[!, :yhat] + width
    #! If the sample size is small, consider using the second-order delta method
    return df
end

# ╔═╡ aa2420c2-8681-49a1-9ad7-87ac6df4466e
"""Return DF of confidence intervals for each temperature for bootstraped simulation results"""
function confidence_interval(df, α)
    df_bs = @chain df begin
        groupby(:tAvg)
        combine(:yhat => mean => :y_mean,
                :yhat => (x -> quantile!(x, α/2)) => :y_lb,
                :yhat => (x -> quantile!(x, 1-α/2)) => :y_ub,
                :yhat => minimum => :y_min, :yhat => maximum => :y_max)
    end
    return df_bs
end

# ╔═╡ ca4b5b53-7eeb-4ee1-b6fd-3306405b78b3
md"### 1.2.3 Functions"

# ╔═╡ 450033bf-3e98-422e-a5c4-a4b934155c12
"""Return dataframe with backward-looking `windowsize` moving averages of all columns"""
function add_moving_avgs!(df; windowsize=10)
    # group each county
    cols = [["fips", "Date"]; names(df, r"temp[^P]")]
    gdf = groupby(df[!, cols], :fips)
    dfs = []
    # for each county, calculate a 10-year backward-looking moving average
    for group in gdf
        # Create TimeArray
        temp = TimeArray(group, timestamp = :Date)
        # create moving average of last 10 years
        temp = moving(mean, temp, windowsize)
        # convert back to dataframe
        temp = DataFrame(temp)
        append!(dfs, [temp])
    end
    # concatenate county dataframes and leftjoin onto df using fips and Date
    df_avg = reduce(vcat, dfs)
    renamingdict = [names(df, r"temp[^P]") .=> names(df, r"temp[^P]") .* "_avg"; ["timestamp" => "Date"]]
    df_avg = DataFrames.rename(df_avg, renamingdict)
    df = leftjoin(df, df_avg, on=[:fips, :Date])
    return df
end

# ╔═╡ e2702c60-89f5-45fd-8903-4c27ea878e28
"""Regress the outcomes on the deviance from the running average bins"""
function reg123(df, y; cluster=:fips, regex=r"temp(?!20to24).*_dev")
    formula_ = term(y) ~ sum([term.(names(df, regex)); [term(fe(:fips)), term(fe(:year))]])
    reg_ = @time reg(df, formula_, Vcov.cluster(cluster))
    se = .√[reg_.vcov[i,i] for i in 1:length(reg_.coef)]
    ub = reg_.coef + 1.96*se
    lb = reg_.coef - 1.96*se
    df2 = DataFrame(tAvg=[-2,2,6,10,14,18,26,30,34],
                    yhat=reg_.coef,
                    yerror=1.96*se,
                    y_ub=ub,
                    y_lb=lb)
    df2 = reduce(vcat, [[df2]; [DataFrame(tAvg=22,yhat=0,yerror=0,y_lb=0,y_ub=0)]])
    df2 = sort!(df2, :tAvg)
    return df2
end

# ╔═╡ 4ecef41b-ace5-49fe-9b67-aee256d4884c
# ╠═╡ show_logs = false
""" Plot deviance bins: above avg # of days this year that were in this bin (days above the 10-year average # of days)"""
function plot123(y, df)
	labels = Dict(:emp_farm_ln => "Farm Employment", :inc_farm_prop_inc_lpc => "Farm Income")
	ylabels = Dict(:emp_farm_ln => "log(Farm Employment)", :inc_farm_prop_inc_lpc => "log(Farm Income per capita)")
	temp = reg123(df, y)
	plot(temp[!,:tAvg], repeat([0], nrow(temp)), c="gray", s=:dash, label="")
	@df temp plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
	fillrange=:y_ub, fillalpha=0.35,
	label="95% CI")
	xlabel!("Average Temperature Bins")
	ylabel!(ylabels[y])
	title!("Deviance from 10-yr Avg Binned Days: $(labels[y])          ")
	@df temp plot!(:tAvg, :yhat, 
		label="Predicted Response", 
		legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
end

# ╔═╡ 02b1800f-6748-49a1-826e-710f4d1f32d6
"""Plot basic bins: # of days in a county-year that day's avg temp was in this bin"""
function plot123b(y, df)
	labels = Dict(:emp_farm_ln => "Farm Employment", :inc_farm_prop_inc_lpc => "Farm Income")
	ylabels = Dict(:emp_farm_ln => "log(Farm Employment)", :inc_farm_prop_inc_lpc => "log(Farm Income per capita)")
	temp = reg123(df, y; regex=r"^temp(?!20to24)[^P].*[^vg]$")
    plot(temp[!,:tAvg], repeat([0], nrow(temp)), c="gray", s=:dash, label="")
    @df temp plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
    xlabel!("Average Temperature Bins")
    ylabel!(ylabels[y])
    title!("Binned Temperature Response: $(labels[y])         ")
    @df temp plot!(:tAvg, :yhat, 
        label="Predicted Response", 
        legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
end

# ╔═╡ 05a1bce4-fe91-4e89-82cd-09a7457331ad
md"### Problem 2 Functions (in order of the problems)"

# ╔═╡ 36f14f2c-330b-473a-a03b-05289d360983
"""2.2.5: Plot separate kernel regressions on either side of the 75-unit threshold."""
function left_right_plot(df_in, xcol, ycol; bandwidth, plotwidth, threshold=75, npoints=200)
	# Plot all observations as circles
	df = dropmissing(df_in, [xcol, ycol])
	ylabels = Dict(:dlhouse => "House Price", :dgtsp => "Pollution")
	plot(df[!, xcol], df[!, ycol], 
		seriestype=:scatter, 
		markeralpha=0.3,
		label="observations",
		xlabel="1974 TSPs levels",
		ylabel="$(ylabels[ycol]) changes",
		title="LOESS regressions: $(ylabels[ycol]) Changes",
		xlims = (75-plotwidth, 75+plotwidth), #75±width
		ylims = extrema(df[!, ycol]),
		legend=:bottomleft
	)
	thickness = 3
	vline!([threshold], label="$xcol threshold", w=thickness)
	# Plot Loess below threshold
	df1 = subset(df, xcol => ByRow(<(threshold)))
	xs1, ys1 = get_loess_prediction(df1, xcol, ycol; bandwidth=bandwidth,  npoints=npoints)
	plot!(xs1, ys1, label="loess regression (below)", w=thickness)
	# Plot Loess above threshold
	df2 = subset(df, xcol => ByRow(≥(threshold)))
	xs2, ys2 = get_loess_prediction(df2, xcol, ycol; bandwidth=bandwidth, npoints=npoints)
	plot!(xs2, ys2, label="loess regression (above)", w=thickness)
end

# ╔═╡ ecbd00f1-ca73-44d8-96ba-56e01dcae8fc
"""2.2.6: Plot unconditional kernal regression and control-conditioned-prediction kernal regression on either side of the 75-unit threshold."""
function controls_smoothed_plot(df, xcol, ycol1, ycol2; bandwidth, plotwidth, threshold=75, npoints=200)
	thickness = 3
    alpha = 0.4
	c3 = 3
	# Plot all observations as circles
	ylabels = Dict(:dlhouse => "House Price", :dgtsp => "Pollution")
	plot(df[!, xcol], df[!, ycol1], seriestype=:scatter, markeralpha=alpha, label="observations", c=1)
	plot!(df[!, xcol], df[!, ycol2], seriestype=:scatter, label="controls' prediction", c=c3, shape=:+, markerstrokewidth=5, markersize=5)
	# Plot Loess below threshold
	df1 = subset(df, xcol => ByRow(<(threshold)))
	x1, y1 = get_loess_prediction(df1, xcol, ycol1; bandwidth=bandwidth,  npoints=npoints)
	x2, y2 = get_loess_prediction(df1, xcol, ycol2; bandwidth=bandwidth,  npoints=npoints)
	plot!(x1, y1,
        label="obseservations conditioned on TSP", 
        xlabel="1974 TSPs levels",
		ylabel="$(ylabels[ycol1]) changes",
		title="LOESS regressions: $(ylabels[ycol1]) Changes",
		xlims = (75-plotwidth, 75+plotwidth), #75±width
		ylims = extrema([df[!, ycol1]; df[!, ycol2]]),
		legend=:topleft,
        w=thickness, c=1)
    plot!(x2, y2, label="predictions conditioned on TSP", w=thickness, c=c3)
	# Plot Loess above threshold
	df2 = subset(df, xcol => ByRow(≥(threshold)))
	x3, y3 = get_loess_prediction(df2, xcol, ycol1; bandwidth=bandwidth, npoints=npoints)
	x4, y4 = get_loess_prediction(df2, xcol, ycol2; bandwidth=bandwidth, npoints=npoints)
	plot!(x3, y3, label="", w=thickness, c=1)
	plot!(x4, y4, label="", w=thickness, c=c3)
    # Vertical line at the threshold
	vline!([threshold], label="$xcol threshold", w=thickness, c=2)
end

# ╔═╡ 3a721fea-a72f-4645-b45a-48f5d03d469f
"""2.2.7: Plot Regulated vs Unregulated kernel regressions for 50-75 1974 TSP."""
function un_regulated_plot(df, xcol, ycol; bandwidth, npoints=2000)
	thickness = 3
    alpha = 0.2
    c3=3
	# Plot all observations as circles
    df1 = @subset(df, :tsp7576 .== 0)
	df2 = @subset(df, :tsp7576 .== 1)
	ylabels = Dict(:dlhouse => "House Price", :dgtsp => "Pollution")
	plot(df1[!, xcol], df1[!, ycol], 
        xlabel="1974 TSPs levels",
        ylabel="$(ylabels[ycol]) changes",
        title="LOESS regressions: $(ylabels[ycol]) Changes",
        xlims = (50, 75),
        ylims = extrema(df[!, ycol]),
        seriestype=:scatter,
        markeralpha=alpha, 
        label="Unregulated observations",
        c=1)
	plot!(df2[!, xcol], df2[!, ycol], seriestype=:scatter, 
        label="Regulated observations", c=c3, shape=:+, markerstrokewidth=5, markersize=5)
	# Plot unregulated
	x1, y1 = get_loess_prediction(df1, xcol, ycol; bandwidth=bandwidth,  npoints=npoints)
	plot!(x1, y1,
        label="Unregulated Kernel Regression", 
		legend=:bottomright,
        w=thickness, c=1)
	# Plot regulated
	x2, y2 = get_loess_prediction(df2, xcol, ycol; bandwidth=bandwidth, npoints=npoints)
	plot!(x2, y2, label="Regulated Kernel Regression", w=thickness, c=c3)
end

# ╔═╡ d9fdb9cf-dbb0-4656-8ded-ad7f7c4d3a1f
md"# Notebook Settings"

# ╔═╡ 9fae08d1-ad8b-47a6-b284-389d7e56e8be
md"Packages"

# ╔═╡ d0b71169-8720-44da-93cb-03b4fd6566cd
md"File Name variables and notebook printing macros"

# ╔═╡ a0fcb7a7-e4b4-4f55-826f-63074f2686b9
# Set local directory paths and filenames
begin
root = dirname(@__FILE__)
zip_fn = "Walker-ProblemSet1-Data.zip"
county_fn = "fips1001.dta"
employment_fn = "reis_combine.dta"
example_fn = "CountyAnnualTemperature1950to2012.dta"
pollution_fn = "poll7080.dta"
example_url = "https://www.dropbox.com/s/fnl1u0ix4e493vv/CountyAnnualTemperature1950to2012.dta?dl=0"

macro seeprints(expr)
	quote
		stdout_bk = stdout
		rd, wr = redirect_stdout()
		$expr
		redirect_stdout(stdout_bk)
		close(wr)
		read(rd, String) |> Text
	end
end
macro with_stdout(expr)
        escaped_expr = esc(expr)
	return quote
		stdout_bk = stdout
		rd, wr = redirect_stdout()
		result = ($escaped_expr)
		redirect_stdout(stdout_bk)
		close(wr)
		print_result = read(rd, String) |> Text
		print_result
	end
end
macro with_stdout_string(expr)
        escaped_expr = esc(expr)
	return quote
		stdout_bk = stdout
		rd, wr = redirect_stdout()
		result = ($escaped_expr)
		redirect_stdout(stdout_bk)
		close(wr)
		print_result = read(rd, String)
		print_result
	end
end
end;

# ╔═╡ 4d0870b4-12bd-4000-a8c4-2caabb5242cc
md"This creates the new file of county-year temperature response variables to regress on (degree days, temperature bins, cubic splines, and piecewise linear). I have compared the values of these county-year variables with the $example_fn file for Autauga County and they match."

# ╔═╡ e786bc92-cb7f-4ab5-8a60-bfe77d9a70df
begin
	extract_file_from_zip(root, zip_fn, pollution_fn)
	pollution = DataFrame(load(joinpath(root, pollution_fn)));
end;

# ╔═╡ 259d5729-692b-4bb2-a89a-baefb962a253
left_right_plot(pollution, :mtspgm74, :dgtsp; bandwidth=bandwidth, plotwidth, npoints=npoints)

# ╔═╡ 7ecb6794-39b9-40b4-b77c-0c15400eea51
left_right_plot(pollution, :mtspgm74, :dgtsp; bandwidth=bandwidth, plotwidth=15, npoints=npoints)

# ╔═╡ c124c041-8bd8-45fc-a32e-7c2b4195c7cd
left_right_plot(pollution, :mtspgm74, :dlhouse; bandwidth=bwdlhouse, plotwidth=pwdlhouse, npoints=npoints)

# ╔═╡ a6dca9d2-e63f-4211-bef7-3699e3fb2f9d
left_right_plot(pollution, :mtspgm74, :dlhouse; bandwidth=bwdlhouse, plotwidth=15, npoints=npoints)

# ╔═╡ 9dd7da67-0f16-4b08-9847-40e20bd0ffac
begin
#! Predict dlhouse using all controls (and not pollution changes)
cols226 = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
pollution2 = dropmissing(pollution, [cols226; [:dlhouse, :mtspgm74, :dgtsp]])
formula226 = (term(:dlhouse) ~ sum(term.(Symbol.(cols226))))
reg226 =  reg(pollution2, formula226, Vcov.robust(); weights=:pop7080)
pollution2[!, :dlhouse_hat] = convert(Vector{Float64}, StatsModels.predict(reg226, pollution2))
end;

# ╔═╡ 2fed03d7-c01a-4710-b550-bb6470ce90eb
controls_smoothed_plot(pollution2, :mtspgm74, :dlhouse, :dlhouse_hat; bandwidth=bwdlhouse, plotwidth=pwdlhouse, npoints=npoints)

# ╔═╡ 76e84935-d551-47b9-97f9-440e67682b95
controls_smoothed_plot(pollution2, :mtspgm74, :dlhouse, :dlhouse_hat; bandwidth=bwdlhouse, plotwidth=15, npoints=npoints)

# ╔═╡ 16e765de-6a0e-4b23-a47d-f8bc12f3eff8
pollution3 = @subset(pollution2, 50 .< :mtspgm74 .< 75);

# ╔═╡ 1d1e7307-c4d6-4c6e-9093-56583596b9ea
un_regulated_plot(pollution3, :mtspgm74, :dgtsp; bandwidth);

# ╔═╡ 35ec6b0f-74a5-472e-963d-24ddf82ca996
un_regulated_plot(pollution3, :mtspgm74, :dgtsp; bandwidth=bwdgtsp)

# ╔═╡ 325ae653-88a6-45e9-91f2-0e679b2593f3
un_regulated_plot(pollution3, :mtspgm74, :dlhouse; bandwidth=bwdlhouse2)

# ╔═╡ 4a83b89f-1ad4-477c-8508-c2a15bf18663
# ╠═╡ show_logs = false
begin
# 
reg211 = reg(pollution, @formula(dlhouse ~ dgtsp), Vcov.robust(); weights=:pop7080)
#
cols221 = [:dgtsp, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
formula221 = (term(:dlhouse) ~ sum(term.(Symbol.(cols221))))
reg221 = reg(pollution, formula221, Vcov.robust(); weights=:pop7080)
#
cols221c = [:dgtsp, :ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :dincome, :dunemp, :dmnfcg]
formula221c = (term(:dlhouse) ~ sum(term.(Symbol.(cols221c))))
reg221c = reg(pollution, formula221c, Vcov.robust(); weights=:pop7080)
means_ = [mean(pollution.dlhouse[reg211.esample]),
		mean(pollution.dlhouse[reg221.esample]),
		mean(pollution.dlhouse[reg221c.esample])]
controls_ = ["N", "Y", "Y"]
shocks_ = ["N", "N", "Y"]
mystats = NamedTuple{(:means, :controls, :shocks)}((means_, controls_, shocks_))
@with_stdout regtable(reg211, reg221, reg221c,
	regressors = ["dgtsp"], print_estimator_section=false,
	custom_statistics=mystats,
	labels = Dict("__LABEL_CUSTOM_STATISTIC_controls__" => "Controls", "__LABEL_CUSTOM_STATISTIC_shocks__" => "Econ. Shocks",
	"__LABEL_CUSTOM_STATISTIC_means__" => "dlhouse Mean"),
	renderSettings = asciiOutput())
end

# ╔═╡ 3ae3a3a0-ce47-49c7-b027-a665129d1807
# ╠═╡ show_logs = false
begin
reg221b = @time reg(pollution, @formula(dgtsp ~ dincome + dunemp + dmnfcg), Vcov.robust(); weights=:pop7080)
@with_stdout regtable(reg221b, renderSettings = asciiOutput())
end

# ╔═╡ bbb320bc-04f6-48f9-a0cb-9e05c360e793
# ╠═╡ show_logs = false
begin
	reg222 = @time reg(pollution, @formula(tsp7576 ~ dincome + dunemp + dmnfcg), Vcov.robust(); weights=:pop7080)
	@with_stdout regtable(reg222, renderSettings = asciiOutput())
end

# ╔═╡ 46146a3f-8c17-4774-8f22-8eeeaa12451f
begin
	# Reg: dgtsp on tsp7576, weight with pop7080
	reg223a = reg(pollution, @formula(dgtsp ~ tsp7576), Vcov.robust(); weights=:pop7080)
	
	# Reg: dgtsp on tsp7576, weight with pop7080, with controls
	cols223b = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :tsp7576]
	formula223b = (term(:dgtsp) ~ sum(term.(Symbol.(cols223b))))
	reg223b = reg(pollution, formula223b, Vcov.robust(); weights=:pop7080)
	
	# Reg: dlhouse on tsp7576, weight with pop7080
	reg223c = reg(pollution, @formula(dlhouse ~ tsp7576), Vcov.robust(); weights=:pop7080)
	
	# Reg: dlhouse on tsp7576, weight with pop7080, with controls
	cols223d = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :tsp7576]
	formula223d = (term(:dlhouse) ~ sum(term.(Symbol.(cols223d))))
	reg223d = reg(pollution, formula223d, Vcov.robust(); weights=:pop7080)
	regtable(reg223c, reg223d, renderSettings = latexOutput("reg223-reducedform.tex"))
	
	# Make Reg Table
	mystats223a = NamedTuple{(:means, :controls)}(
		([mean(pollution.dgtsp[reg223a.esample]),
			mean(pollution.dgtsp[reg223b.esample]),
			mean(pollution.dlhouse[reg223c.esample]),
			mean(pollution.dlhouse[reg223d.esample])], 
		["N", "Y", "N", "Y"]))
	@with_stdout regtable(reg223a, reg223b, reg223c, reg223d,
		regressors = ["tsp7576"], print_estimator_section=false,
		custom_statistics=mystats223a,
		labels = Dict("__LABEL_CUSTOM_STATISTIC_controls__" => "Controls", "__LABEL_CUSTOM_STATISTIC_shocks__" => "Econ. Shocks",
		"__LABEL_CUSTOM_STATISTIC_means__" => "dlhouse Mean"),
		renderSettings = asciiOutput())
end

# ╔═╡ ba1525cf-44a4-4704-94c0-9dae18cb2b67
begin
	# 2SLS Reg: dlhouse on dgtsp, weight with pop7080
	reg223e = reg(pollution, @formula(dlhouse ~ (dgtsp ~ tsp7576)), Vcov.robust(); weights=:pop7080)
	
	# 2SLS Reg: dlhouse on dgtsp, weight with pop7080, with controls
	cols223f = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
	formula223f = (term(:dlhouse) ~ sum(term.(Symbol.(cols223f))) + (term(:dgtsp) ~ term(:tsp7576)))
	reg223f = reg(pollution, formula223f, Vcov.robust(); weights=:pop7080)
	regtable(reg223e, reg223f, renderSettings = latexOutput("reg223-IV.tex"))
	
	# Make Reg Table
	mystats223b = NamedTuple{(:ymeans, :xmeans, :controls)}(
		([mean(pollution.dlhouse[reg223e.esample]),
			mean(pollution.dlhouse[reg223f.esample])],
		[mean(pollution.dgtsp[reg223e.esample]),
			mean(pollution.dgtsp[reg223f.esample])],
		["N", "Y"]))
	@with_stdout regtable(reg223e, reg223f,
		regressors = ["dgtsp"], print_estimator_section=false,
		custom_statistics=mystats223b,
		labels = Dict("__LABEL_CUSTOM_STATISTIC_controls__" => "Controls", "__LABEL_CUSTOM_STATISTIC_shocks__" => "Econ. Shocks",
		"__LABEL_CUSTOM_STATISTIC_ymeans__" => "dlhouse Mean",
		"__LABEL_CUSTOM_STATISTIC_xmeans__" => "dgtsp Mean"),
		renderSettings = asciiOutput(),
		regression_statistics = [:nobs, :r2, :adjr2, :f, :f_kp])
end

# ╔═╡ a33ac919-c4f3-4ce8-89cc-d7adbbb8d044
begin
	Δx = -round(mean(pollution.dgtsp[reg223f.esample]), digits=1)
	Δy = round(mean(pollution.dlhouse[reg223f.esample]), digits=3)
	Δlny = round((ℯ^mean(pollution.dlhouse[reg223f.esample])-1)*100,digits=1)
	β100 = -round(100*last(reg223f.coef), digits=1)
	eβ100 = round(100*(ℯ^last(reg223f.coef)-1), digits=4)
end;

# ╔═╡ ebce0eca-ff0b-4e07-bfcd-b0e5f9691ae3
Markdown.parse("

Then, take an average county (County X) that saw a decrease of $Δx units of pollution (an improvement in air quality) and an increase in housing prices by about  \$100(e^{$Δy}-1)\$  = $Δlny% between 1970 and 1980. If instead they saw one unit more of pollution improvement (a $(Δx+1)-unit improvement instead), then we would expect that County X housing prices would have instead increased by

$Δlny% - (100\$\\beta_{IV}\$)% \$\\approx\$ $Δlny% + $β100% = $(round(Δlny+β100, digits=1))%

showing evidence that, on average, pollution decreases caused by the regulation have increased house prices more than they would have otherwise increased.
")

# ╔═╡ d97e9d6b-9050-487b-a6aa-b947c607a38e
last(reg223e.coef) ≈ last(reg223c.coef) / last(reg223a.coef)

# ╔═╡ 62d53ad4-5a70-4311-80da-0e5e11269453
last(reg223f.coef) ≈ last(reg223d.coef) / last(reg223b.coef)

# ╔═╡ ac0f080d-3068-48ca-9995-dfbd7fb122df
begin
	# Create an indicator for mtspgm74<75
	pollution[!, :above75] = pollution[!, :mtspgm74] .> 75	
	
	#! first-stage relationship between regulation and air pollution changes
	#! Relevance
	# Reg: dgtsp on above75, weight with pop7080
	reg224a = reg(pollution, @formula(dgtsp ~ above75), Vcov.robust(); weights=:pop7080)
	
	# Reg: dgtsp on above75, weight with pop7080, with controls
	cols224b = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :above75]
	formula224b = (term(:dgtsp) ~ sum(term.(Symbol.(cols224b))))
	reg224b = reg(pollution, formula224b, Vcov.robust(); weights=:pop7080)
	regtable(reg224a, reg224b, renderSettings = latexOutput("reg224-firststage.tex"))
	
	#! reduced form relationship between regulation and housing price changes
	#! Reduced form estimate
	# Reg: dlhouse on above75, weight with pop7080
	reg224c = reg(pollution, @formula(dlhouse ~ above75), Vcov.robust(); weights=:pop7080)
	
	# Reg: dlhouse on above75, weight with pop7080, with controls
	cols224d = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop, :above75]
	formula224d = (term(:dlhouse) ~ sum(term.(Symbol.(cols224d))))
	reg224d = reg(pollution, formula224d, Vcov.robust(); weights=:pop7080)
	
	# Make Reg Table
	mystats224a = NamedTuple{(:means, :controls)}(
		([mean(pollution.dgtsp[reg224a.esample]),
			mean(pollution.dgtsp[reg224b.esample]),
			mean(pollution.dlhouse[reg224c.esample]),
			mean(pollution.dlhouse[reg224d.esample])], 
		["N", "Y", "N", "Y"]))
	@with_stdout regtable(reg224a, reg224b, reg224c, reg224d,
		regressors = ["above75"], print_estimator_section=false,
		custom_statistics=mystats224a,
		labels = Dict("__LABEL_CUSTOM_STATISTIC_controls__" => "Controls", "__LABEL_CUSTOM_STATISTIC_shocks__" => "Econ. Shocks",
		"__LABEL_CUSTOM_STATISTIC_means__" => "dlhouse Mean"),
		renderSettings = asciiOutput())
end

# ╔═╡ d279d5d2-c82b-4084-897a-7ecf8d969229
begin
	#! air quality changes on housing price changes using two-stage least squares and the above75 indicator as an instrument
	#! IV estimate
	# 2SLS Reg: dlhouse on dgtsp (above75 as IV), weight with pop7080
	reg224e = reg(pollution, @formula(dlhouse ~ (dgtsp ~ above75)), Vcov.robust(); weights=:pop7080)
	
	# 2SLS Reg: dlhouse on dgtsp (above75 as IV), weight with pop7080, with controls
	cols224f = [:ddens, :dwhite, :dfeml, :dage65, :dhs, :dcoll, :durban, :dpoverty, :dvacant, :downer, :dplumb, :dtaxprop]
	formula224f = (term(:dlhouse) ~ sum(term.(Symbol.(cols224f))) + (term(:dgtsp) ~ term(:above75)))
	reg224f = reg(pollution, formula224f, Vcov.robust(); weights=:pop7080)
	
	# Make Reg Table
	mystats224b = NamedTuple{(:ymeans, :xmeans, :controls)}(
		([mean(pollution.dlhouse[reg224e.esample]),
			mean(pollution.dlhouse[reg224f.esample])],
		[mean(pollution.dgtsp[reg224e.esample]),
			mean(pollution.dgtsp[reg224f.esample])],
		["N", "Y"]))
	@with_stdout regtable(reg224e, reg224f,
		regressors = ["dgtsp"], print_estimator_section=false,
		custom_statistics=mystats224b,
		labels = Dict("__LABEL_CUSTOM_STATISTIC_controls__" => "Controls", "__LABEL_CUSTOM_STATISTIC_shocks__" => "Econ. Shocks",
		"__LABEL_CUSTOM_STATISTIC_ymeans__" => "dlhouse Mean",
		"__LABEL_CUSTOM_STATISTIC_xmeans__" => "dgtsp Mean"),
		renderSettings = asciiOutput(),
		regression_statistics = [:nobs, :r2, :adjr2, :f, :f_kp])
end

# ╔═╡ ff5cef88-d4d7-42b7-8bd4-5dd591b08def
"""Read file at save_path to dataframe, if not present, download from url."""
function df_from_url(url, save_path)
    if isfile(save_path)
        println("Loading $save_path from file.")
        return DataFrame(load(joinpath(root, example_fn)))
    else
        println("Downloading $save_path from url.")
        download(url, joinpath(root, example_fn))
        println("Loading $save_path from file.")
        return DataFrame(load(joinpath(root, example_fn)))
    end
end

# ╔═╡ 0d56492d-b6c8-4585-8454-8ffe3517d340
# ╠═╡ show_logs = false
begin
	# Extract and read employment file from zip
	df_employ = df_from_zip(root, zip_fn, employment_fn)
	df_employ[!, :fips] = parse.(Int32, df_employ[!, :fips])
	
	# Download and load example county file to compare variables
	df_temp = df_from_url(example_url, joinpath(root, example_fn))
	
	# drop first row because it's missing
	df_temp = last(df_temp, nrow(df_temp)-1)
	
	# Join onto climate data
	df = innerjoin(df_temp, df_employ, on=[:fips, :year])
	
	# Make year and fips categegorical for regressions
	df[!,:year] = categorical(df[!,:year])
	df[!,:fips] = categorical(df[!,:fips])
	
	# Create log(farm employees in county) and replace log(0) with missing
	df[!, :emp_farm_ln] = replace(log.(df[!, :emp_farm]), -Inf => missing)
	
	# Create log tranformed inc_farm_prop_income/population 
	logfn(x) = x===missing ? missing : (x ≤ 0 ? missing : log(x))
	df[!, :inc_farm_prop_inc_lpc] = logfn.(df[!,:inc_farm_prop_income] ./ df[!,:pop_population])
	# 40,385 missing
	
	# Convert the year into a Date
	df[!, :Date] = Date.(convert(Vector{Int}, df[!,:year]))

	# For each bin, add the 10-year moving average
	df = add_moving_avgs!(df)
	# Calculate the deviance from the running mean
	for col in names(df, r"temp[^P].*(?<!_avg)$")
	    df[!, "$(col)_dev"] = df[!, col] - df[!, "$(col)_avg"]
	end
end;

# ╔═╡ 309ea09a-12ae-4a5a-b099-0aa366206856
# ╠═╡ show_logs = false
begin
	formula121 = @formula(emp_farm_ln ~ tempB0 + temp0to4 + temp4to8 + temp8to12 + temp12to16 + temp16to20 + temp24to28 + temp28to32 + tempA32 +
	    fe(fips) + fe(year)
	)
	# Reference category = temp20to24 (days in county-year where tAvg ∈ (20, 24]°C)
	
	reg121 = reg(df, formula121, Vcov.cluster(:fips))
	se121 = .√[reg121.vcov[i,i] for i in 1:length(reg121.coef)]
	ub121 = reg121.coef + 1.96*se121
	lb121 = reg121.coef - 1.96*se121
	df121 = DataFrame(tAvg=[-2,2,6,10,14,18,26,30,34],
	                  yhat=reg121.coef,
	                  yerror=1.96*se121,
	                  y_ub=ub121,
	                  y_lb=lb121)
	df121 = reduce(vcat, [[df121]; [DataFrame(tAvg=22,yhat=0,yerror=0,y_lb=0,y_ub=0)]])
	df121 = sort!(df121, :tAvg)
	
	plot(df121[!,:tAvg], repeat([0], nrow(df121)), c="gray", s=:dash, label="")
	@df df121 plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
	    fillrange=:y_ub, fillalpha=0.35,
	    label="95% CI")
	xlabel!("Average Temperature Bins")
	ylabel!("log(Farm Employment)")
	title!("Binned Temperature Response: Farm Employment      ")
	p121 = @df df121 plot!(:tAvg, :yhat, 
	    label="Predicted Response", 
	    legend=:bottomright, c=2, shape=:circle, markerstrokewidth=0)
end

# ╔═╡ 75ffadd1-7c23-4a5b-9a96-aa2c86a7d831
# ╠═╡ show_logs = false
begin
	formula122 = @formula(inc_farm_prop_inc_lpc ~ splineC1 + splineC2 + splineC3 + splineC4 + fe(year) + fe(fips))
	reg122 = reg(df, formula122, Vcov.cluster(:fips))
	# Initialize dataframe used for plotting
	# - x-axis with temperatures 0-40°C at 0.25° intervals
	# - spline variables
	df_plot = plotdf_init(reg122)
	# Use delta method on recentered data to produce confidence intervals
	# - outcome centered on predicted outcome at 20°C
	df_plot = deltamethod_CIbounds!(df_plot, reg122)
end;

# ╔═╡ 3a91edc6-18ba-4f0b-9cf9-5aa3b00022e1
begin
# Plot Delta method confidence intervals and predicted outcome
plot(df_plot[!,:tAvg], repeat([0], nrow(df_plot)), c="gray", s=:dash, label="")
@df df_plot plot!(:tAvg, :y_lb, w=0, msw = 0, ms = 0, c=1,
    fillrange=:y_ub, fillalpha=0.35,
    label="95% CI")
xlabel!("Daily Average Temperature")
ylabel!("log(Farm Income per capita)")
title!("Spline Temperature Response: Farm Employment      ")
@df df_plot plot!(:tAvg, :yhat, label="Predicted Response", legend=:bottomleft, c=2)
end

# ╔═╡ a86ea794-73ad-4012-99c2-12778898c56d
# Let's try bootstrapping!
"""Given a dataframe sample and dataframe with tAvg values to predict on,
    estimate the predicted values vector (for bootstrapping)."""
function prediction_sample(df, df_plot)
    reg_ = reg(df, formula122, Vcov.cluster(:fips))
    df2 = DataFrame(tAvg=0:0.25:40)
    df2[!, :yhat] = Matrix(df_plot[!,r"splinerelativeC"]) * reg_.coef
    return df2[!, [:tAvg, :yhat]]
end

# ╔═╡ 1cc81fad-6309-4865-9576-e2fb637a02b6
"""Return dataframe of bootstrapped predictions with α confidence intervals."""
function bootstrap_predicted_outcome(df::DataFrame,
                                     prediction_function::Function;
                                     screenwidth=50,
                                     nsimulations=10000,
                                     α=0.05
                                     )
    println("\nStarting simple bootstrap")
    cols = [:year, :fips, :tAvg, :inc_farm_prop_inc_lpc, :splineC1, :splineC2, :splineC3, :splineC4]
    preds = []
    t = @timed for k ∈ 1:nsimulations
        k%ceil(nsimulations/screenwidth) == 0 ? print("⋅") : nothing
        sample_rows = sample(1:nrow(df), nrow(df), replace=true)
        df_sample = df[sample_rows, cols]
        append!(preds, [prediction_function(df_sample, df_plot)])
    end
    # Vertically concatenate the bootstrap predictions
    preds_comb = reduce(vcat, preds)
    # Get mean and 95% CI for each value of tAvg
    df_bs = confidence_interval(preds_comb, α)
    print("$(t.time) seconds\n")
    return df_bs
end

# ╔═╡ f11d7396-cde5-47fc-9c4f-71c58cd5f906
"""Return dataframe of cluster-bootstrapped predictions with α confidence intervals.

    Clustered Bootstrap (over each fips county code)
    Need to sample 3055 counties of 3055 counties (with replacement)
    Requires repeating some counties in the data and keeping all years
    of their data
"""
function bootstrap_predicted_outcome(df::DataFrame,
                                     prediction_function::Function,
                                     cluster::Union{Symbol, String};
                                     screenwidth=50,
                                     nsimulations=10000,
                                     α=0.05
                                     )
    println("\nStarting clustered bootstrap")
    cols = [:year, :fips, :tAvg, :inc_farm_prop_inc_lpc, :splineC1, :splineC2, :splineC3, :splineC4]
    preds = []
    df_grouped = groupby(df[!, cols], cluster)
    ngroups = length(df_grouped)
    t = @timed for k ∈ 1:nsimulations
        k%ceil(nsimulations/screenwidth) == 0 ? print("⋅") : nothing
        sample_idx = sample(1:ngroups, ngroups; replace = true, ordered = false)
        df_sample = reduce(vcat, [df_grouped[i] for i in sample_idx])
        append!(preds, [prediction_function(df_sample, df_plot)])
    end
    # Vertically concatenate the bootstrap predictions
    preds_comb = reduce(vcat, preds)
    # Get mean and 95% CI for each value of tAvg
    df_bs = confidence_interval(preds_comb, α)
    print("$(t.time) seconds\n")
    return df_bs
end

# ╔═╡ e38183bf-7403-49b7-987d-d4f6df7562e5
begin
"""Bootstrapping confidence intervals for Problem 1"""
function run_bootstrap(df, df_plot, prediction_sample)
	# Get mean and 95% CI for each value of tAvg using bootstraps
	nsimulations = 100
	# Simple bootstrap (over all individuals)
	df_bs = bootstrap_predicted_outcome(df, prediction_sample, nsimulations=nsimulations)
	# Clustered bootstrap
	df_bs_cluster = bootstrap_predicted_outcome(df, prediction_sample, :fips; nsimulations=nsimulations)   
	
	
	# Plot Delta method and bootstraps confidence intervals and predicted outcome
	s = 0; dpi=400;
	@df df_plot plot(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=1,
		fillrange=:y_ub, fillalpha=0.35,
		label="95% CI Delta Method")
	@df df_bs plot!(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=2,
		fillrange=:y_ub, fillalpha=0.35,
		label="95% CI BootStrap ($nsimulations)")
	@df df_bs_cluster plot!(:tAvg, :y_lb, w=0, msw = 0, ms = s, c=3,
		fillrange=:y_ub, fillalpha=0.35,
		label="95% CI BootStrap-clustered ($nsimulations, FIPS)")
	p8 = @df df_plot plot!(:tAvg, :yhat, label="Predicted Response", legend=:bottomleft)
	savefig(p8, "plot122-delta-boot-bootcluster.svg")
end
end

# ╔═╡ 37882765-26f6-4e68-af76-e64a5a5068fe
# ╠═╡ show_logs = false
plot123(:emp_farm_ln, df)

# ╔═╡ e044b2e2-35db-4b3d-9bb0-de63cb791767
# ╠═╡ show_logs = false
plot123(:inc_farm_prop_inc_lpc, df)

# ╔═╡ 8f0a8f38-5d62-47e9-a647-e7b607149c35
# ╠═╡ show_logs = false
plot123b(:emp_farm_ln, df)

# ╔═╡ c4261949-5b6a-4727-90e0-bde6ca907186
# ╠═╡ show_logs = false
plot123b(:inc_farm_prop_inc_lpc, df)

# ╔═╡ 8879e813-c350-4fe0-8c92-44a80e361588
"""Save temperature variables
    1.1.1 Construct 4 temperature response variables
    1.1.2 Aggregate to year-county
    - sum each new temperaature variable in each grid point over each year
    - average tMin, tMax, tAvg, perc in each gridpoint over each year
    - take a simple average over all grid points in the county to get county-year level
"""
function create_temperature_vars()
    # Extract single-county fips file from zip if not present, and read
    extract_file_from_zip(root, zip_fn, county_fn)

    # Load data and create average daily temp
    df_temp1 = DataFrame(load(joinpath(root, county_fn)))
    df_temp1[!, :tAvg] = (df_temp1[!, :tMax] .+ df_temp1[!, :tMin]) ./ 2

    # Add degree days for each day-gridpoint
    thresholds = [30, 32, 34]
    df_temp1 = add_degree_days(df_temp1, thresholds)

    # Add binned temperature variables
    bin_edges = [0, 4, 8, 12, 16, 20, 24, 28, 32]
    df_temp1 = add_bin_indicator_days(df_temp1, bin_edges)
    df_temp1 = add_bin_portion_days(df_temp1, bin_edges)

    # Add cubic spline basis variables for each day-gridpoint
    knots = [0 8 16 24 32]
    df_temp1 = add_cubic_spline_vars(df_temp1, :tAvg, knots)

    # Add piecewise linear basis variables for each day-gridpoint
    knots = [28, 32]
    df_temp1 = add_piecewise_linear_vars(df_temp1, :tAvg, knots)

    # sum over year, then average over all grid points in county
    column_tags = ["dday", "temp", "splineC", "piece"]
    df_temp = aggregate_df(df_temp1, column_tags)

    # Save the new file
    CSV.write(joinpath(root, "CountyAnnualTemperature_aaron.csv"), df_temp)
end

# ╔═╡ 93470fdb-b32d-4cff-9205-717b0d873447
if !isfile(joinpath(root, "CountyAnnualTemperature_aaron.csv"))
    create_temperature_vars()
end

# ╔═╡ f46dd00c-afd6-4f48-b860-b326c3d303f7
md"Table of contents sidebar"

# ╔═╡ 621d4ab7-680e-4ae2-b1d0-4a43817c81b4
begin
	# Create table of contents on right side of screen
	# TableOfContents()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
CovarianceMatrices = "60f91f6f-d783-54cb-84f9-544141854719"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FixedEffectModels = "9d5cd8c9-2029-5cab-9928-427838db53e3"
GLM = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Loess = "4345ca2d-374a-55d4-8d30-97f9976e7612"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RDatasets = "ce6b1742-4840-55fa-b093-852dadbb1d8b"
RegressionTables = "d519eb52-b820-54da-95a6-98e1306fdade"
StatFiles = "1463e38c-9381-5320-bcd4-4134955f093a"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsModels = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
TimeSeries = "9e3dc215-6440-5c97-bce1-76c03772f85e"
ZipFile = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"

[compat]
CSV = "~0.10.4"
CUDA = "~3.9.0"
CategoricalArrays = "~0.10.5"
Chain = "~0.4.10"
CovarianceMatrices = "~0.10.4"
DataFrames = "~1.3.3"
DataFramesMeta = "~0.11.0"
Distributions = "~0.25.53"
FixedEffectModels = "~1.6.5"
GLM = "~1.7.0"
LaTeXStrings = "~1.3.0"
Latexify = "~0.15.14"
Loess = "~0.5.4"
Plots = "~1.27.6"
PlutoUI = "~0.7.38"
RDatasets = "~0.7.7"
RegressionTables = "~0.5.4"
StatFiles = "~0.8.0"
StatsBase = "~0.33.16"
StatsModels = "~0.6.29"
StatsPlots = "~0.14.33"
TimeSeries = "~0.23.0"
ZipFile = "~0.9.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "a598ecb0d717092b5539dbbe890c98bac842b072"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.2.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CompilerSupportLibraries_jll", "ExprTools", "GPUArrays", "GPUCompiler", "LLVM", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "ba75320aaa092b3e17c020a2d8b9e0a572dbfa6a"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "3.9.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "109664d3a6f2202b1225478335ea8fea3cd8706b"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.5"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CovarianceMatrices]]
deps = ["CategoricalArrays", "Documenter", "LinearAlgebra", "Requires", "Statistics", "StatsBase", "StatsModels", "Tables"]
git-tree-sha1 = "f9184166d19489a365a5f816088195cae5306bef"
uuid = "60f91f6f-d783-54cb-84f9-544141854719"
version = "0.10.4"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "6c19003824cbebd804a51211fd3bbd81bf1ecad5"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.3"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "f1d89a07475dc4b03c08543d1c6b4b2945f33eca"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.11.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5a4168170ede913a2cd679e53c2123cb4b889795"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.53"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "Base64", "Dates", "DocStringExtensions", "IOCapture", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Markdown", "REPL", "Test", "Unicode"]
git-tree-sha1 = "6edbf28671b4df4f692e54ae72f1e35851cfbf38"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "0.27.16"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FixedEffectModels]]
deps = ["DataFrames", "FixedEffects", "LinearAlgebra", "Printf", "Reexport", "Statistics", "StatsBase", "StatsFuns", "StatsModels", "Tables", "Vcov"]
git-tree-sha1 = "bd741fd9b058179064ede5adaef9424d6ebbf946"
uuid = "9d5cd8c9-2029-5cab-9928-427838db53e3"
version = "1.6.5"

[[deps.FixedEffects]]
deps = ["GroupedArrays", "LinearAlgebra", "Printf", "Requires", "StatsBase"]
git-tree-sha1 = "63fdaffe29fdd3dc7c96de87ce890d4398ec7661"
uuid = "c8885935-8500-56a7-9867-7708b20db0eb"
version = "2.1.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "92b8d38886445d6d06e5f13201e57d018c4ff880"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.7.0"

[[deps.GPUArrays]]
deps = ["Adapt", "LLVM", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "c783e8883028bf26fb05ed4022c450ef44edd875"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.3.2"

[[deps.GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "556190e1e0ea3e37d83059fc9aa576f1e2104375"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.14.1"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.GroupedArrays]]
deps = ["DataAPI", "Missings"]
git-tree-sha1 = "93e21548d0a4b8ac793fea1aa1d720f5c9eaf11a"
uuid = "6407cd72-fade-4a84-8a1e-56e431fc1533"
version = "0.3.1"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "c9b86064be5ae0f63e50816a5a90b08c474507ae"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.9.1"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5558ad3c8972d602451efe9d81c78ec14ef4f5ef"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.14+2"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "6f14549f7760d84b2db7a9b10b88cd3cc3025730"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.14"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a970d55c2ad8084ca317a4658ba6ce99b7523571"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.12"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "7008a3412d823e29d370ddc77411d593bd8a3d03"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.9.1"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3b429f37de37f1fc603cc1de4a799dc7fbe4c0b6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "6f2dd1cf7a4bbf4f305a0d8750e351cb46dfbe80"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "28ef6c7ce353f0b35d0df0d5930e0d072c1f5b9b"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.RData]]
deps = ["CategoricalArrays", "CodecZlib", "DataFrames", "Dates", "FileIO", "Requires", "TimeZones", "Unicode"]
git-tree-sha1 = "19e47a495dfb7240eb44dc6971d660f7e4244a72"
uuid = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
version = "0.8.3"

[[deps.RDatasets]]
deps = ["CSV", "CodecZlib", "DataFrames", "FileIO", "Printf", "RData", "Reexport"]
git-tree-sha1 = "2720e6f6afb3e562ccb70a6b62f8f308ff810333"
uuid = "ce6b1742-4840-55fa-b093-852dadbb1d8b"
version = "0.7.7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "afeacaecf4ed1649555a19cb2cad3c141bbc9474"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.5.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.ReadStat]]
deps = ["DataValues", "Dates", "ReadStat_jll"]
git-tree-sha1 = "f8652515b68572d3362ee38e32245249413fb2d7"
uuid = "d71aba96-b539-5138-91ee-935c3ee1374c"
version = "1.1.1"

[[deps.ReadStat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "afd287b1031406b3ec5d835a60b388ceb041bb63"
uuid = "a4dc8951-f1cc-5499-9034-9ec1c3e64557"
version = "1.1.5+0"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegressionTables]]
deps = ["Compat", "Distributions", "FixedEffectModels", "Formatting", "GLM", "Statistics", "StatsBase", "StatsModels", "UnPack"]
git-tree-sha1 = "caec9d037af643dfc511a3923102a5bdbfcd704a"
uuid = "d519eb52-b820-54da-95a6-98e1306fdade"
version = "0.5.4"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StatFiles]]
deps = ["DataValues", "FileIO", "IterableTables", "IteratorInterfaceExtensions", "ReadStat", "TableShowUtils", "TableTraits", "TableTraitsUtils", "Test"]
git-tree-sha1 = "28466ea10caec61c476a262172319d2edf248187"
uuid = "1463e38c-9381-5320-bcd4-4134955f093a"
version = "0.8.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "03c99c7ef267c8526953cafe3c4239656693b8ab"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.29"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "4d9c69d65f1b270ad092de0abe13e859b8c55cad"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.14.33"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeSeries]]
deps = ["Dates", "DelimitedFiles", "DocStringExtensions", "RecipesBase", "Reexport", "Statistics", "Tables"]
git-tree-sha1 = "3c91141a9f2276c37c3b6bc2bd83e652d50fecbc"
uuid = "9e3dc215-6440-5c97-bce1-76c03772f85e"
version = "0.23.0"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "0a359b0ee27e4fbc90d9b3da1f48ddc6f98a0c9e"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.3"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "11db03dd5bbc0d2b57a570d228a0f34538c586b1"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.17"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Vcov]]
deps = ["Combinatorics", "GroupedArrays", "LinearAlgebra", "StatsBase", "Tables"]
git-tree-sha1 = "b382811c8beba117f70c07a42f3f18b3075a39db"
uuid = "ec2bfdc2-55df-4fc9-b9ae-4958c2cf2486"
version = "0.5.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "505c31f585405fc375d99d02588f6ceaba791241"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─220c7746-b558-4b42-9e47-e227730a1aad
# ╟─8a325a10-4c17-43a1-a75f-9bc68f72fbdb
# ╟─00edd77c-be64-44fc-b510-28234152ca56
# ╟─ac7d2e34-53a5-4c0e-b413-c1129378d91a
# ╟─22681583-0918-49da-87b8-3a86c027cf81
# ╟─b0a24d4d-016f-45a4-9124-90d463d1db93
# ╠═93470fdb-b32d-4cff-9205-717b0d873447
# ╟─4d0870b4-12bd-4000-a8c4-2caabb5242cc
# ╟─28d5ff08-5fbd-45c8-9816-dacbe9d20493
# ╟─f9fbde9e-b9fd-4407-a63d-992a28408245
# ╟─65cebd89-47ba-4773-9b7f-98ce0ef42d4b
# ╟─c94aac41-7071-4cbd-acd6-9f1cd53176cd
# ╟─225f183d-aa3b-46ce-8edf-4336bba8da30
# ╟─7214407b-7850-4a69-8798-59c463394cf2
# ╟─309ea09a-12ae-4a5a-b099-0aa366206856
# ╟─0d56492d-b6c8-4585-8454-8ffe3517d340
# ╟─9cf39312-272b-4b73-867c-faf98da039f7
# ╟─fe21bed1-dcd3-47fa-997c-9827db649e3f
# ╟─75ffadd1-7c23-4a5b-9a96-aa2c86a7d831
# ╟─84fbd2e7-82cd-45ea-bb89-d75a621f79ca
# ╟─3a91edc6-18ba-4f0b-9cf9-5aa3b00022e1
# ╟─aa60b5ef-3aeb-42ff-b679-e00b28205267
# ╟─6a978c31-51c4-47a8-b626-14ca00de0370
# ╟─952ba66d-3264-405e-ad52-9f0a951e244c
# ╠═37882765-26f6-4e68-af76-e64a5a5068fe
# ╠═e044b2e2-35db-4b3d-9bb0-de63cb791767
# ╟─9a4d0f65-7c09-4c92-99fb-ff31b0453054
# ╠═8f0a8f38-5d62-47e9-a647-e7b607149c35
# ╠═c4261949-5b6a-4727-90e0-bde6ca907186
# ╟─68e3f5e2-382a-431e-a19e-3547e70c3a94
# ╟─b460634b-72c1-4f56-ad00-2a550f1c660f
# ╟─288e1ee5-809f-4108-80d1-15b5c8e67452
# ╟─4fe46d77-791f-4e80-b715-2ccb12de37f5
# ╟─714b920e-a30a-4930-bd05-fb070dc4f641
# ╟─8255f716-16ad-44c2-9c87-1710f95ff1f6
# ╟─478897a9-13e7-453a-ad0a-4b4031e7932d
# ╟─e786bc92-cb7f-4ab5-8a60-bfe77d9a70df
# ╟─a613f88c-ea92-4e66-8a62-6f9012437cb9
# ╟─71b0a94f-67ae-4240-ba7a-5e9ebd33af45
# ╟─4a83b89f-1ad4-477c-8508-c2a15bf18663
# ╟─cac486e4-3bd6-45ab-b6d0-24dbafc0da71
# ╟─3ae3a3a0-ce47-49c7-b027-a665129d1807
# ╟─af154adc-56f4-45d8-9c5c-7647fc7ed12e
# ╟─b27cc951-9bd5-45f6-b30f-51010f2c9537
# ╟─fb3ad8e9-08f9-4ffa-891b-e6bd926897a7
# ╟─bd0603ce-3a15-4428-af4d-c2ee2127e537
# ╟─d084de41-4f40-4e51-8dca-54ec0608e72b
# ╟─1185d096-15c8-42e1-bbff-8047722dbabc
# ╟─bbb320bc-04f6-48f9-a0cb-9e05c360e793
# ╟─0df4996b-c5e2-4ae1-9909-09b442cc9ff2
# ╟─adbf35de-93e8-492d-a89d-c61697da3391
# ╟─bc11f842-1b19-4250-a969-81b85a9f7dc4
# ╟─d9face0e-651d-4256-8aab-3c257a18619d
# ╟─6bbe502f-3750-4321-87c1-a01d67e228a1
# ╟─95412c61-3a1b-40bf-9b7c-ddd54ce17ba7
# ╟─46146a3f-8c17-4774-8f22-8eeeaa12451f
# ╟─0cc1daaa-9358-4074-a576-3f4263c18e95
# ╟─ebce0eca-ff0b-4e07-bfcd-b0e5f9691ae3
# ╟─a33ac919-c4f3-4ce8-89cc-d7adbbb8d044
# ╟─6c00497f-91bf-4122-a366-e2bb119bbcdc
# ╟─ba1525cf-44a4-4704-94c0-9dae18cb2b67
# ╟─76650858-6af6-4060-9487-4d6424dc7cb7
# ╠═d97e9d6b-9050-487b-a6aa-b947c607a38e
# ╠═62d53ad4-5a70-4311-80da-0e5e11269453
# ╟─c6e1abd9-0490-40cf-a8db-9dfe3bd8249a
# ╟─26a2778a-fe46-4772-8e78-e09b63ae4152
# ╟─f89066c8-39d0-4e85-b924-3e493aabda88
# ╟─22802ab4-ed98-40fe-9c20-4207843c2cb2
# ╟─4aaa6d2c-1ea7-400b-a8e1-9e0993c9a5f7
# ╟─ac0f080d-3068-48ca-9995-dfbd7fb122df
# ╟─b989539a-09f4-4ecc-8f24-d3efb5962035
# ╟─d279d5d2-c82b-4084-897a-7ecf8d969229
# ╟─ae43ae55-ad11-43a9-9f25-9123ee28a885
# ╟─90903d07-605a-4b51-91e7-b2f5bac56c29
# ╟─2f203087-0dab-49a9-9115-e3b8477d3089
# ╟─78ac5bab-75c5-4ccf-aa63-1681bce62b9c
# ╟─573af536-2a9e-4f7a-8a68-0b5e3157bddf
# ╟─2145f8ab-2e3a-4b88-803b-f73ffdd0b4e6
# ╟─4df0e4c6-55f4-401d-b18e-c4cba40d0de4
# ╟─85e03d2c-af57-4feb-9107-dedd0d59ecd0
# ╟─f7aaaee5-262c-40cc-8ff6-793a904253fa
# ╟─3ab934a9-4dc0-473b-be1d-f290761dddd8
# ╟─2d5817b1-a1cd-4882-9c3d-ac53eb9362f2
# ╟─259d5729-692b-4bb2-a89a-baefb962a253
# ╟─0c8a7d25-7233-45d5-b3c1-f46bd423cb94
# ╟─7ecb6794-39b9-40b4-b77c-0c15400eea51
# ╟─eb368d82-03b9-47a0-839f-b90f8eb10e09
# ╟─b59f26c1-9792-4958-86ae-3712f19854fa
# ╟─9036f9b8-e157-495b-b2a5-b583935170f0
# ╟─c124c041-8bd8-45fc-a32e-7c2b4195c7cd
# ╟─44aeda75-7419-40cb-9234-f01945e4d9b9
# ╟─a6dca9d2-e63f-4211-bef7-3699e3fb2f9d
# ╟─f0a30424-670b-4200-98bd-e9f4ecd38eb2
# ╟─834fe856-c053-4076-901c-754eae2c8846
# ╟─d4dee95d-f991-44ea-9c52-95d51f8c1ab2
# ╟─bd95a3c4-553e-418f-9831-dda93068a9be
# ╟─202a5a35-7633-44cc-887e-22c73d3019d0
# ╟─9dd7da67-0f16-4b08-9847-40e20bd0ffac
# ╟─2fed03d7-c01a-4710-b550-bb6470ce90eb
# ╟─76e84935-d551-47b9-97f9-440e67682b95
# ╟─4b46b382-5349-4bc7-bb03-229c303f37be
# ╟─98857139-6cd1-4e39-9e1e-8a35f2629e28
# ╟─8dc1d528-a9a1-4f99-8b83-a864e5beac6b
# ╟─bdf52a46-4f47-4bc5-9e07-a20fb7b76dbf
# ╟─ce91dbb9-70ec-43e5-978a-68cd7817a6fd
# ╟─16e765de-6a0e-4b23-a47d-f8bc12f3eff8
# ╟─1d1e7307-c4d6-4c6e-9093-56583596b9ea
# ╟─948f96b3-6f50-4687-8d4f-db309c8e29a2
# ╟─35ec6b0f-74a5-472e-963d-24ddf82ca996
# ╟─16d94c66-4a4e-443a-9f0b-ea2af6e4555f
# ╟─325ae653-88a6-45e9-91f2-0e679b2593f3
# ╟─e5981e73-8029-4e7b-aedb-76b1080d882a
# ╟─4321d60b-fc24-4acf-9184-f7a68e7909d8
# ╟─9d9db365-0799-4699-8c75-c16b1820d322
# ╟─40b4be47-8c52-466a-8130-cfc9b256452d
# ╟─da4a328f-ef3f-44d4-90cb-b22b1519fad0
# ╠═57754134-0b34-41b5-9366-fe7e469ee1af
# ╟─872ffa35-9f86-461c-bf77-29f4a8ad5e26
# ╟─b5cf0297-5489-488f-a75f-586adca9f222
# ╟─730c1994-cbee-4ec9-95a8-51e0923ef4ed
# ╟─26d3b991-cc5c-455e-b38d-5925069cd713
# ╠═5f31d109-b75e-4cc5-9cae-a06493009221
# ╟─96ba58aa-8028-4db7-bd76-ef3fc7f24637
# ╟─01ce9f00-209a-4f4b-bd63-d7b02767f5a0
# ╟─8b3041fe-3b6d-4ef2-9223-5f81a6366019
# ╟─7d7b39f2-1db6-4c0b-a41b-710030b7537d
# ╠═8893d6e2-2e3a-440a-ba2b-e1b60a5b9d73
# ╠═892333c3-0c5e-4d58-bc0d-84502afe4c96
# ╠═ff5cef88-d4d7-42b7-8bd4-5dd591b08def
# ╠═af4f8032-c318-4035-913f-351c8e0130f7
# ╟─93eae3b1-187b-4026-81e7-a587d1477af9
# ╠═8879e813-c350-4fe0-8c92-44a80e361588
# ╠═b0bf6cbc-e1e4-42f1-8d0c-055f5a249faa
# ╠═cdc1dd7f-683c-40c8-9410-3fcdb4d74322
# ╠═a72bfb4b-0a93-48c5-81f7-6156c88c4ef4
# ╠═79bf9c02-b2a1-40a2-b87e-d436d587fe13
# ╠═ac20872f-f7e9-480f-8155-b1da7c616909
# ╠═0e4a2b18-9a0f-4acd-a4e8-920a1bac0b87
# ╠═24c043c5-dc9b-4a7a-a815-bdd4bb1997bb
# ╠═790ce5d9-6119-4c5c-b057-6e0a52ecdaf2
# ╠═b9a9e3eb-f707-4e25-87de-9d6d53a9db60
# ╠═b6218121-c5cb-4d53-bfab-69a10a8aa241
# ╠═fe94d889-551e-47ea-8363-6c6064a3df1b
# ╠═07e6f171-e8bb-4d96-acc8-0d0694692141
# ╠═78d8e678-4148-451d-a2de-9d84986d522b
# ╠═779e117a-3e39-4e34-9e09-d6d10c8a4f22
# ╠═d18c8a79-d9c5-4efb-aeae-c2e58e07251c
# ╠═0cd2f8f9-a72d-4600-89c5-8ee99f3009f7
# ╠═bcfdbe23-5290-40fa-9965-a1f254c21bee
# ╟─fec1fd19-b0e0-48c4-9ade-140c0c068379
# ╟─de2d137d-a9cc-43b2-ba9f-caa0502c6210
# ╠═45a1bad3-2d2c-409b-a742-8401ff15f166
# ╠═17ec6db5-aa8c-4488-acfb-595c245f2c23
# ╠═4ac46ab5-1615-4feb-bfb6-bc9cef26c353
# ╠═a86ea794-73ad-4012-99c2-12778898c56d
# ╠═aa2420c2-8681-49a1-9ad7-87ac6df4466e
# ╠═1cc81fad-6309-4865-9576-e2fb637a02b6
# ╠═f11d7396-cde5-47fc-9c4f-71c58cd5f906
# ╠═e38183bf-7403-49b7-987d-d4f6df7562e5
# ╟─ca4b5b53-7eeb-4ee1-b6fd-3306405b78b3
# ╠═450033bf-3e98-422e-a5c4-a4b934155c12
# ╠═e2702c60-89f5-45fd-8903-4c27ea878e28
# ╠═4ecef41b-ace5-49fe-9b67-aee256d4884c
# ╠═02b1800f-6748-49a1-826e-710f4d1f32d6
# ╟─05a1bce4-fe91-4e89-82cd-09a7457331ad
# ╠═36f14f2c-330b-473a-a03b-05289d360983
# ╠═ecbd00f1-ca73-44d8-96ba-56e01dcae8fc
# ╠═3a721fea-a72f-4645-b45a-48f5d03d469f
# ╟─d9fdb9cf-dbb0-4656-8ded-ad7f7c4d3a1f
# ╟─9fae08d1-ad8b-47a6-b284-389d7e56e8be
# ╟─9bf60494-b58e-495e-adc0-8c9d0585880b
# ╟─d0b71169-8720-44da-93cb-03b4fd6566cd
# ╟─a0fcb7a7-e4b4-4f55-826f-63074f2686b9
# ╟─f46dd00c-afd6-4f48-b860-b326c3d303f7
# ╟─621d4ab7-680e-4ae2-b1d0-4a43817c81b4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
