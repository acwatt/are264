# 2022-03-03
## Prompt
Take the example of OPower information treatments. If you wanted to establish the welfare effects of this treatment, how might you do it? (Here I'm getting at the idea that the nudge might induce guilt versus provide valuable information versus some other channel. The experiment reveals the behavioral response, but how should we interpret that in terms of the welfare of the people in the experiment?)

## Response
I imagine four changes in utility for the energy consumer that result from the informational nudge. Because the consumer assumes their changes in consumption probably won't affect the total amount of CO2 or other pollutants coming from energy use in their area, I assume that they have no change in utility from decreased externalities from energy use. The four effects on utility I imagine that result from a decrease in energy use or a shift to a "better" time to consume energy are: (1) decreases in utility from shifting or decreasing consumption (because they were probably at a weak optimum in terms of the timing and amount of energy consumed before the nudge), (2) decrease in utility from the guilt of not being a good environmentalist, (3) increase in utility from income savings from energy bill savings (if decreasing consumption or shifting to a time that has a lower rate), and (4) increase in utility from the pride associated with being a better environmentalist and shifting/decreasing consumption. 

If the household changes their consumption in reaction to the nudge, which we would see in the data, then I assume that the absolute value of increases in utility (3 and 4) is larger than decreases in utility (1 and 2) so that the new choice is utility-maximizing. From this, it's hard to get a very clear measure of welfare changes. But I think we could at least sign the changes: if people change their behavior after the nudge, their utility has increased (positive welfare change); if they do not change their behavior, then their utility has weakly decreased from guilt but has not changed w.r.t. points 1, 3, or 4.



# 2022-03-01
## Prompt
Upload a quick file or enter text explaining the main contribution of Myers, Puller and West. Assume that the reader you are addressing is an applied economist with no particular expertise in energy or behavioral economics.

## Response
Of interest to both policy makers and economists, this paper provides some of the first evidence on efficacy of mandatory disclosure programs. Specifically, in a program where it is mandatory to disclose the energy efficiency of a house before sale, this paper shows evidence that the program increases investments in energy efficiency and thus increases the quality of the good being sold.

Additionally, of more general interest to economists, this paper provides evidence from a structural behavior model that the market failure being corrected by the policy is that of symmetrically incomplete information where both the seller and buyer are missing information about the energy efficiency of the house relative to other houses (as opposed to asymmetric where the seller has more information than the buyer). This is useful because it helps further decompose the class of incomplete information market failures and provide evidence for a policy instrument that could alleviate this particular form of incomplete information.




# 2022-02-24
## Prompt
Goulder, Hafstead and Williams 2016 
What does this article find that is surprising?
What do you think about the nature of the evidence presented in this paper?

## Response
Goulder, Hafstead, and Williams find that there are plausible circumstances where an efficiency standard can actually out-perform a carbon tax or cap and trade policy on cost-effectiveness grounds. Using a model that incorporates the complexity of the US tax system and the energy market, and calibrates on real input-output tables of international trade, they find that efficiency standards can outperform carbon taxes due to the tax-interaction effect.

The results, both analytical and numerical, seem like existence proofs -- that there are some cases in which CES can do better than a tax. These results are not necessarily robust to most situations, but the fact that there are reasonable conditions in which CES could be more cost-effective was important to show so that economists don't walk around badmouthing standards when they actually could be more useful than a carbon tax given the complicated existing market distortions already in place in the US.



# 2022-02-17
## Prompt
see problem sets/2_writing-workshop/WorkshopPrompt.pdf

## Response
See pre-notes pdf





# 2022-02-09
## Prompt
Below is a table from Muller and Mendelsohn (AER 2009) that describes the distribution of marginal damages from a ton of each pollutant. The distribution is across points sources in the United States. The differences derive from an air transport model that maps how emissions at each point interact with the ambient environment and determine exposure to people.

Your task is to think about the implications of this type of heterogeneity in damages for the design of a corrective policies. For example, suppose that you were designing an SO2 permit program. What are the implications of this type of heterogeneity for such a program?

## Response
Generally, the economist's purpose of instituting a market-based policy solution to an environmental externality is to create incentives that arrive at an efficient solution -- where the efficient outcome is where marginal damages of the externality caused by the service or good are equal to the marginal benefits. In a permit program, the permits are generally tradable and bankable such that each permit is equal to the same amount of marginal damages.

If we assume that marginal damages are equal across sources of pollution, then a simple permit trading program can result in an efficent outcome. For example, if marginal damages from an additional ton of CO2 are identical no matter where the CO2 is generated (becuase most of the damages are caused by the CO2 rising into the upper atomosphere where it does long-run damage), then the damage done by the CO2 allowed by each permit is independent of the firm who emits the CO2 and independent of how many permits the firm has.

Assume we want to design an SO2 permit trading program (Cap and Trade, for example). If we know that marginal damages of SO2 are heterogeneous across point sources, then we know we need to adjust the permits so when they are traded across firms/point sources, the marginal damages are equal. Assume that we know the approximate marginal damages per ton of emissions for each firm. Then we can assign a marginal damage weight or exchange coefficent to each firm, such that one permit of firm j, times the exchange coefficient of firm j, is equal to the marginal damage of one ton of SO2 produced by that firm (or some multiple of marginal damage). So when firm j goes to trade permits with firm i, they need to convert between j and i permits. There would then be an exchange rate between every pair of firms, and if we wanted, we could define an n by n exchange matrix (where n is the number of firms regulated under the program).

For example, if firm j's marginal damage is 200 $/ton/year and firm i's marginal damage is 1000 $/ton/year, then one of firm i's permits would be worth 5 of firm j's permits. This is all assuming both firms have constant marginal damage.

We could make this even more realistic and assume that marginal damages are not constant -- that firm j's 10th ton of SO2 produces higher social damages than it's first ton. The only modification needed is that the i,j element of the n by n exchange matrix would need to be a function of firms i and j's total amount of permits/pollution.






# 2022-02-03
## Prompt
There are several channels through which (a) pollution or (b) policies that mitigate pollution can have a welfare impact on individuals. If we want to describe the distributional impacts of pollution (or policies that try to reduce it), we need to have some sense of these different channels and their relative importance. Think about how you would construct a taxonomy of these channels and come up with a list. Record that list here. (Ideally you will do this before completing the relevant reading. To give away that game, the Fullerton article provides one such taxonomy, which we will discuss.)

## Response
### (a) Channels through which pollution impacts welfare
(+) benefits of the service/good that is causing the pollution (like the 
    electricity from the coal plant down the street)
(+) cheaper prices for those goods
(-) costs of defensive investments (buying an air purifier, reduced 
    outdoor activities)
(-) health impacts (decreased lung capacity, high blood pressure, asthma, 
    reduced health from reduced outdoor activities, etc)
(-) costs of avoiding pollution / search costs (searching for alternative bike routes to work to
    avoid car exhaust near the busy roads, vacationing somewhere without smoke)
(-) increased uncertainty about your health (more doctor visits, more stress)

### (b) Channels through which mitigating policies impacts welfare
(-) direct cost of the policy (taxes, tolls, etc)
(-) indirect cost of policy (increased prices of regulated or downstream goods)
(+) decreased cost of individual mitigations
(+) health impacts
(+) long-run decreased costs of cleaner alternatives due to shifts in investments
    away from the polluting technology
(+) increase the number of alternatives in the market 
(+) decreased health uncertainty

