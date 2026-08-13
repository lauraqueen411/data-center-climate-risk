# Physical Climate Risk to Data Centers and their Critical Dependencies Across Oregon

## Introduction

Data centers are the physical infrastructure behind digital services, like cloud computing and the training and running of artificial intelligence (AI) models, that human systems increasingly rely on. They range from small server rooms serving a single business to hyperscale campuses covering millions of square feet. At their core, they house the same basic components regardless of scale: servers and storage hardware, networking equipment and cabling, and the power and cooling systems needed to keep it all running. Oregon is home to one of the earliest hyperscale data centers in the U.S. - Google's facility in The Dalles, which opened in 2006 and currently occupies over 1.3 million square feet.

AI adoption is accelerating both data center construction and the energy and water demand that comes with it. National electricity demand from data centers is projected to grow from 4.4% of total U.S. consumption in 2023 to as much as 12% by 2028, alongside a comparable rise in water use: direct water consumption for cooling is projected to roughly double to quadruple, from about 17 billion gallons in 2023 to 38-73 billion gallons by 2028 [1]. Data centers rely on energy and water to keep running: constant server load needs continuous, redundant power, and the heat that load generates needs substantial water and energy to remove, whether through water-cooled chillers, cooling towers, or humidification systems [2,3].

Climate change poses material risk to data center operations and the infrastructure systems that support them. The IPCC defines climate risk as the potential for adverse consequences to human or ecological systems from climate hazards interacting with exposure and vulnerability [4]. It breaks the concept into three parts: hazard (a physical event or trend, like a heatwave or flood), exposure (people, assets, or systems present in places that could be affected), and vulnerability (a system's predisposition to be harmed, often from a lack of capacity to cope or adapt). Financial and corporate frameworks built on this science, such as the Task Force on Climate-related Financial Disclosures [5], split physical risk into two types: acute risk from short-term extreme events like cyclones or flash floods, and chronic risk from long-term shifts like sustained temperature increase, drought, or sea-level rise.

Quantifying physical climate risk has become a fast-growing industry, led largely by climate tech firms (XDI, First Street) and financial institutions, but much of this work relies on proprietary data, especially at the asset and vulnerability level. Academic literature, by contrast, has focused mostly on the reverse relationship: how data centers affect the climate, through electricity consumption, water withdrawal, and emissions (e.g., [6,7,8]). This study focuses on how climate change affects data center operations and the infrastructure they depend on, rather than how data centers contribute to climate change.

Two recent studies have assessed data center climate risk at the scale of the contiguous U.S. Kollar & Grady [2] evaluated 2,445 data center locations for heat wave exposure, water stress, and the emissions intensity of the local electricity grid, using chi-squared tests and Moran's I to identify spatial clustering. They found heat risk was spatially correlated with both emissions and water stress, with Texas showing consistently high risk across all three. Their own framing describes data centers as embedded in a network of water, electricity, and climate feedbacks, yet the risk metrics they compute stay anchored to each facility's own coordinates; the authors acknowledge that their approach misses operator-level network diversification, which could offset locational risk.

Esparza et al. [3] took a related but distinct approach, pairing FEMA National Risk Index data across six hazard types with EAGLE-I historical power outage records for 2,660 data centers, using LISA and bivariate LISA analysis to identify hazard hotspots and rank counties for future siting. Including historical outage frequency and duration is a step toward capturing grid-mediated vulnerability, but the outage data are aggregated at the county level and reflect past history rather than a traced pathway from a specific hazard to a specific upstream asset - a substation or generation plant - that a given data center depends on. The authors flag this gap themselves, calling for future work on "how other infrastructure systems interact with a DC" at a more granular, asset-specific level.

Together, these two studies show that data center hazard exposure varies meaningfully and clusters spatially across the country, and both point to infrastructure dependency as a critical, unaddressed piece of vulnerability without actually measuring it. Both are also grounded in historical hazard data rather than future climate projections, so neither captures how this exposure is expected to change as the climate continues to warm.

Two constraints shape this study's approach. First, as-built engineering data for individual facilities - the exact cooling system model, mechanical specifications, or transformer loading history for a given data center or substation - is proprietary and not publicly available. In its place, this study applies published engineering standards and vulnerability functions: ASHRAE's thermal operating thresholds, IEEE and ISO derating curves for transformers, transmission lines, and generators, and FEMA's Hazus depth-damage functions. These are genuine engineering-based vulnerability models, just generic ones substituted for facility- and asset-specific specifications that aren't publicly obtainable.

Second, facility-level engineering vulnerability alone doesn't capture how climate hazards actually threaten data center operations, since a facility can fail because a hazard hits an upstream asset it depends on rather than the building itself. This motivates extending the analysis to the power infrastructure - substations, transmission lines, and generators - each facility depends on, following the infrastructure interdependency framework of Rinaldi, Peerenboom & Kelly [9]. These constraints and complexities reflect the interdisciplinary nature of the problem, drawing on climate and hazard science, infrastructure engineering, and critical infrastructure systems research to characterize risk that no single field fully captures alone.

This study assesses two hazards - heat and flood - at both the facility level and the level of the power infrastructure (substations, transmission lines, and backup generators) each facility depends on. For heat, this means acute risk (equipment failure above a temperature threshold) and chronic risk (rising cooling costs), plus the temperature-driven capacity loss of dependent power infrastructure. For flood, this means facility-level exposure and expected structural damage, an average annual damage estimate integrated across return periods, and the exposure of dependent power infrastructure to the same flood hazard. Water and connectivity dependencies, and wildfire hazard, are natural extensions of this framework but fall outside the current scope [for now]. The goal is to move beyond facility-level hazard exposure toward a more systemic view of climate vulnerability - one that reflects how risk actually reaches data center operations through the infrastructure they depend on - and to provide a methodological foundation for both future research and practical planning for Oregon's data center sector.

## Methods

### Data Sources

| Category | Dataset / Source | Key Specs | Used For |
|---|---|---|---|
| Data Center Exposure | IM3/PNNL Data Center Inventory | List of current and projected future U.S. data center locations, sizes, and status, compiled by Pacific Northwest National Laboratory | Master list of data center sites analyzed throughout the study (1.1, 1.2, 1.3, 2.2) |
| | HIFLD | Federal database of power substation and transmission line locations | Identifies the power grid infrastructure each data center depends on (1.3, 2.5) |
| | EIA-860 | Federal database of power plant and backup generator locations, sizes, and fuel types | Identifies the backup generators each data center depends on (1.3, 2.5) |
| Heat Hazard | LOCA2 | Statistically downscaled CMIP6 future temperature projections; daily high and low temperature (tasmax, tasmin), historical and future periods | Provides the temperature data used throughout the heat analysis (1.1, 1.2, 1.3) |
| | MACA | Statistically downscaled CMIP6 humidity projections; daily high and low relative humidity (rhsmin, rhsmax), historical and future periods | Provides the humidity data needed to estimate cooling system performance (1.2) |
| Heat Vulnerability | ASHRAE TC9.9 | Industry engineering guideline for the temperature range data center cooling systems are designed to handle | Sets the temperature threshold used to estimate when a facility's cooling would fail (1.1) |
| | Lei & Masanet PUE model | Open-source engineering model estimating data center cooling efficiency (Power Usage Effectiveness, or PUE) under different weather conditions and cooling technologies | Estimates how much less efficient data center cooling becomes as the climate warms (1.2) |
| Flood Hazard | JRC flood depth grids | Global historic flood depths at ~90m resolution for seven return periods ranging from 1-in-10 to 1-in-500 years | Provides the baseline flood depth used to assess current flood exposure (2.1) |
| | ISIMIP flood depth | Lower-resolution (~25km) flood depth maps covering both past and future conditions | Estimates how flood depths may change in the future, applied to the more detailed baseline maps (2.1) |
| | USGS 3DEP | High-resolution (1-meter or better) digital elevation model | Provides precise ground elevation at each facility, used to calculate flood depth relative to the facility itself (2.2, 2.3) |
| Flood Vulnerability | FEMA Hazus depth-damage functions | Standard federal formulas relating flood depth to expected building damage | Converts flood depth into an estimate of expected building damage, and into an average annual damage estimate (2.3, 2.4) |

### 1. Heat Risk

#### 1.1 Acute Heat Risk: Facility Operational Failure

We use daily high and low temperatures from LOCA2 climate projections at each Oregon data center, for a historical baseline period and future period (SSP5-8.5, mid-century horizon, ~30-year window). We track daily lows alongside highs because cooling systems need cool nights to recover during multi-day heat events. For each facility and year, we identify the hottest day and count the number of days exceeding a fixed temperature threshold, comparing historical and future periods.

The threshold comes from ASHRAE's published temperature ranges for data center cooling equipment [10]: a recommended range of 18-27°C across all classes, and a wider allowable range that varies by equipment class. We use the Class A2 allowable range (10-35°C), the range most modern equipment is built to handle [10], and adjust it downward at higher-elevation sites where cooling performance drops. We assume Class A2 equipment across all facilities. Crossing the threshold does not mean physical damage - it means cooling can no longer keep pace, putting the facility at risk of failure.

Output: exceedance-day counts and peak severity, historical vs. future, showing how failure risk changes over time.

#### 1.2 Chronic Heat Risk: Cooling Load

Data centers use extra electricity to run their cooling systems, and this cost rises as outside temperatures increase. We measure this using Power Usage Effectiveness (PUE): the ratio of total facility electricity use to the electricity used by computing equipment alone. A higher PUE means more electricity is spent on cooling and overhead rather than computing.

We estimate PUE with an open-source engineering model (Lei & Masanet [8,11]) that calculates cooling efficiency from local weather and cooling system type. This model is well established in the field: Lei and Masanet are themselves co-authors on the 2024 LBNL Data Center Energy Usage Report [1], which incorporates the same climate- and technology-specific PUE approach at the national scale. The model needs both temperature and humidity, because evaporative cooling - a common energy-saving method - works less effectively in humid air. We pair LOCA2 temperature data with humidity data from MACA, a comparable regional climate dataset, to calculate this combined measure for each facility.

The model supports several cooling system designs; we use the version for outside-air and evaporative cooling, the design most Oregon facilities are known to use [1]. We check the model's historical output against real-world efficiency values reported for Oregon's two main electric utility regions [12] to confirm it produces realistic results before projecting forward.

Output: for each facility, the projected change in PUE from historical to future conditions - a direct estimate of how much less efficient cooling becomes as the region warms. This feeds into Section 1.3: a higher PUE means more electricity draw for the same computing load, adding pressure to the same power infrastructure already strained by heat during extreme events.

#### 1.3 Critical Infrastructure Heat Dependency: Substation and Generator Derating

[this section I am less confident in / not sure the network-approach will be as valuable as flood + wildfire, where risk by location could meaningfully change.]

Power infrastructure loses capacity in extreme heat - the same conditions that push facility cooling demand to its highest. We apply LOCA2 temperature projections to the substations, transmission lines, and backup generators each data center depends on, to estimate how much capacity they lose during future heat events.

For substations and transformers, we use published derating rates from IEEE C57.91 [13]: capacity drops by about 1.5% per 1°C above a 30°C reference temperature for self-cooled units, or 1% per °C for forced-air-cooled units. For transmission lines, we use a simplified version of the standard current-rating method (IEEE 738 [14]), which estimates lost current-carrying capacity as ambient temperature rises above a line's design rating; this simplified version is less precise than a full engineering calculation and is documented as a limitation. For backup generators, capacity loss is estimated from standard reference conditions (ISO 8528-1 [15]) and typical manufacturer derating rates for temperature and altitude; exact rates vary by generator make and model, so this is flagged as an assumption to refine with manufacturer-specific data where available.

Output: for each facility, whether projected extreme heat days coincide with reduced capacity at its dependent substations, transmission lines, or generators - a compound risk beyond facility-level cooling failure alone.

### 2. Flood Risk

#### 2.1 Flood Hazard Data and Future Scaling

We characterize flood hazard using two datasets. JRC flood depth maps provide flood depth at ~90m resolution for seven return periods (1-in-10 to 1-in-500 years) during the historical period, and serve as our baseline. ISIMIP flood depth data, available at much coarser resolution (~25km) for both historical and future periods, is used only to calculate how flood depths are expected to change - not as a standalone hazard layer, since its resolution is too coarse for site-level analysis.

For each data center location and return period, we calculate a change factor: the ratio between ISIMIP's future and historical flood depth at the nearest grid cell. We apply this change factor to the detailed JRC baseline map to produce a future flood depth estimate that keeps JRC's spatial detail while reflecting ISIMIP's projected change.

#### 2.2 Facility-Level Flood Exposure

We extract the scaled flood depth at each data center's location, for each return period and time horizon. We use a high-resolution digital elevation model (USGS 3DEP) to determine each facility's precise ground elevation, so flood depth is calculated relative to the facility site itself rather than a coarse or generic reference point.

Output: flood depth by return period and time horizon for each facility, and a flag for whether depth exceeds a minimum damage threshold.

#### 2.3 Vulnerability: Depth-Damage Translation

Flood depth alone doesn't tell us how much damage to expect. We convert it to an expected loss using FEMA's Hazus depth-damage functions [16], which relate flood depth to percent structural damage for different building types. Data centers aren't a specific Hazus building category, so we use the closest available category (e.g., industrial or essential-facility) as a stand-in, and test this choice against at least one alternative category.

Hazus damage functions are measured relative to a building's first floor, not ground level. We combine our elevation data with a standard first-floor height offset to make this conversion, since as-built floor elevations aren't available for most facilities; this offset is flagged as an assumption to refine if facility-specific values become available. We apply the selected damage function to the resulting depth at each facility, for each return period and time horizon, to estimate percent structural damage. We use only the structural damage curve, not the separate contents-damage curve, consistent with excluding equipment like servers from the facility's replacement-cost basis; adding contents damage is a possible future extension.

Output: percent structural damage per facility, by return period and time horizon.

#### 2.4 Average Annual Damage

A single return period only describes damage from one flood severity. To get a single number describing overall flood risk, we combine all return periods into average annual damage (AAD) - a standard flood-risk metric estimating the damage expected in any given year, averaged across all possible flood sizes and their likelihoods [16,17].

Each return period corresponds to an annual probability of occurring (a 100-year flood has a 1% chance of happening in a given year). We plot the damage estimate from Section 2.3 against this probability for each return period, then calculate the area under that curve - integrating damage over probability - to get AAD. Because we only have damage estimates at a handful of return periods, we approximate this area using the trapezoidal rule between each pair of points. We assume no damage at the most frequent end of the curve (an event with roughly a 1-year return period) and hold damage constant beyond the largest return period we have data for; both are documented as simplifications.

Output: average annual damage per facility, calculated separately for historical and future flood conditions - a single comparable measure of how flood risk changes over time, the flood equivalent of the PUE trend calculated for heat risk in Section 1.2.

#### 2.5 Critical Infrastructure Flood Dependency

We apply the same future flood depth estimate to the substations, generators, and transmission lines each data center depends on. This catches a failure mode that facility-level analysis alone would miss: a data center building sitting outside the flood zone but losing power because a substation floods instead.

Output: which data centers have at least one critical dependency - a substation, generator, or transmission line - exposed to flood depths above the damage threshold, regardless of whether the facility itself is directly exposed.

## References

[1] Shehabi, A., et al. (2024). 2024 United States Data Center Energy Usage Report. Lawrence Berkeley National Laboratory.
[2] Kollar, A., & Grady, C. (2025). The relationship between data centers and the climate is a systems challenge: a spatial analysis of United States data centers. Environmental Research Communications, 7(11), 111005. https://doi.org/10.1088/2515-7620/ae193a
[3] Esparza, M. T., Li, B., Ma, J., & Mostafavi, A. (2025). AI meets natural hazard risk: A nationwide vulnerability assessment of data centers to natural hazards and power outages. International Journal of Disaster Risk Reduction, 126, 105616.
[4] IPCC. (2022). Climate Change 2022: Impacts, Adaptation and Vulnerability. Contribution of Working Group II to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change. Cambridge University Press.
[5] Task Force on Climate-related Financial Disclosures (TCFD). (2017). Final Report: Recommendations of the Task Force on Climate-related Financial Disclosures. Financial Stability Board.
[6] Siddik, M. A. B., Shehabi, A., & Marston, L. (2021). The environmental footprint of data centers in the United States. Environmental Research Letters, 16(6), 064017. https://doi.org/10.1088/1748-9326/abfba1
[7] Shehabi, A., et al. (2011). [Citation not independently verified - see note below.]
[8] Lei, N., & Masanet, E. (2022). Climate- and technology-specific PUE and WUE estimations for U.S. data centers using a hybrid statistical and thermodynamics-based approach. Resources, Conservation and Recycling, 182, 106323. https://doi.org/10.1016/j.resconrec.2022.106323
[9] Rinaldi, S. M., Peerenboom, J. P., & Kelly, T. K. (2001). Identifying, understanding, and analyzing critical infrastructure interdependencies. IEEE Control Systems Magazine, 21(6), 11-25. https://doi.org/10.1109/37.969131
[10] ASHRAE Technical Committee 9.9. Thermal Guidelines for Data Processing Environments, Table 2 (2015 edition, SI version).
[11] Lei, N., & Masanet, E. (2020). Statistical analysis for predicting location-specific data center PUE and its improvement potential. Energy, 197, 117556. https://doi.org/10.1016/j.energy.2020.117556
[12] Guidi, et al. (2026). Assessing the Carbon Emissions and Energy Consumption of U.S. Hyperscale Data Centers. arXiv:2606.05420.
[13] IEEE Std C57.91-2011. IEEE Guide for Loading Mineral-Oil-Immersed Transformers and Step-Voltage Regulators. Institute of Electrical and Electronics Engineers.
[14] IEEE Std 738-2023. IEEE Standard for Calculating the Current-Temperature Relationship of Bare Overhead Conductors. Institute of Electrical and Electronics Engineers.
[15] ISO 8528-1:2018. Reciprocating Internal Combustion Engine Driven Alternating Current Generating Sets - Part 1: Application, Ratings and Performance. International Organization for Standardization.
[16] Federal Emergency Management Agency (FEMA). Hazus Flood Model Technical Manual, Hazus 6.1.
[17] U.S. Army Corps of Engineers. Risk-Based Analysis for Flood Damage Reduction Studies, EM 1110-2-1619.
