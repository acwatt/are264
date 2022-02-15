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

