# Boiling Point for Depth (BPD) Curve

Construction of BPD curves with salinity and gas corrections. These curves define the temperature of saturation (i.e., when liquid and vapour water are present) for a given pressure.

This project was initiated during the World Geothermal Congress May 2021 Hackathon organised by [Agile](https://agilescientific.com/). The project was scoped by Irene Wallis, [Cubic Earth](https://www.cubicearth.nz/). Development is led by [Thorsten Hörbrand](https://github.com/thoerbr) with contributions from [Jan Niederau](https://github.com/Japhiolite), and [Irene Wallis](https://github.com/ICWallis). This BPD tool is in its infancy and contributions to its continued development are welcome. 

***
## BPD Curves in Geothermal Resource and Well Analysis

Development of the BPD curve concept, as it is used in geothermal resource analysis, was informed by early drilling at Yellowstone National Park in Wyoming and at Steamboat Springs in Nevada. In these systems, the reservoir temperature was found to be close to boiling for the corresponding hydrostatic pressure (White, 1968b; White et al., 1968). Since then, the BPD curve has become a valuable tool for both reservoir and well analysis in two-phase liquid or vapour dominated geothermal resources. 

_Classification for geothermal systems by enthalpy (from Kaya et al., 2011)_

![system-types](https://github.com/Geothermal-Fluids/bpd_curve/blob/main/system-types.png)

During exploration, the BPD curve may used to approximate the upper limit of temperature with depth by assuming a hydrostatic pressure gradient, which is a reasonable assumption for geothermal reservoirs. Prior to drilling and if springs or fumaroles are available for sampling, the reservoir temperature is estimated using geothermometers (i.e., fluid-mineral equilibria that are modelled or empirically fit to derive a numerical expression that relates temperature to constituent concentrations). The combination of these temperature estimates and the BPD curve will inform the pre-drilling conceptual model and support the evaluation of drilling risk from steam/two-phase conditions. 

As seen in the figre below, the shallowest point of the BPD curve is tied to the water table which can be 10s or 100s of meters below the ground surface. Typically the elevation of discharging chloride hot springs are used to ascertain the reservoir liquid level prior to drilling.

_Sketch BPD curve adapted from Nicholson (1993)_

![BPD-Curve-Concept](https://github.com/Geothermal-Fluids/bpd_curve/blob/main/bpd_concept.png)

Assuming adiabatic flow, a hypothetical reservoir where 250°C waters rise buoyantly through a permeable rock would be liquid until meeting saturation conditions at around 400 m below the water table. In the natural (pre-development) state, where the rock is in thermal equilibrium with the fluid, we can assume that fluid flow is adiabatic (i.e., that there is no energy loss). In reality, because flow through the reservoir involves friction, there is some energy loss. Once the rising fluids meet saturation conditions, boiling ensues. As fluid continues to rise, the temperature will change and this change follows the sweep of the BPD curve until the water table is met. In some high-enthalpy reservoirs, a steam zone with associated isobaric pressures may be present above the liquid level and below the reservoir cap.  

The BPD curve informs analysis of temperature logs acquired in geothermal wells. During well analysis, the measured pressure is used to determine saturation conditions. Just like the reservoir, the liquid level inside a well may be a long way below the wellhead. Depending on the well condition while the pressure - temperature log was acquired, the liquid level inside the well may or may not approximate the liquid level inside the reservoir. For further information on how BPD are used during temperature log analysis, refer to [this](https://github.com/ICWallis/T21-Tutorial-WellTestAnalysis) tutorial on geothermal well test analysis. 

Increasing salinity (Na-Ca-K-Cl) or gas (dominantly CO2 and H2S) concentrations will have opposing effects on depth to saturation, and the methods in this repository enable correction for these effects. Increased gas concentration will allow boiling to start at greater depths whereas an increased salinity suppresses the boiling point, requiring the fluid to ascend to shallower depths before boiling. Most geothermal reservoirs have low concentrations of salts and small changes in salinity have little impact on the BPD curve. In contrast, a shift in gas concentration will have a great impact on BPD. Subsequently, boiling zones can be much deeper in high-gas geothermal reservoirs than in those that are gas-poor (Nicholson, 1994). 

Steam bubbles in a column of fluid will also result in a steeper BPD gradient, and therefore deeper boiling, because they lower the weight (pressure) of the fluid above (Haas, 1971). Canet et al. (2011) modelled the impact of excluding vapour bubbles in a BPD curve for the shallow subsurface (<550 m). They found that for temperatures between 200 and 225°C, as would be expected in the shallow subsurface of high-temperature reservoirs, exclusion of vapour bubbles when deriving the BPD condition leads to a c. 40-50% underestimate of the depth to boiling. Although worthy of consideration, the effect of steam bubbles is beyond the scope of this repository.

For the majority of high-temperature geothermal reservoirs, uncertainty generated by unknown concentrations of gas and salts have little effect on reservoir-scale analysis, such as conceptual model development. Furthermore, reservoir temperature estimates derived using geothermometers will typically have a greater uncertainty than the effect that minor variations in salt or gas concentrations would have on the BPD curve.  

There are, however, a small number of high-temperature geothermal reservoirs that contain sufficient salt and gas concentrations that it is useful to include the effect in pre-drilling BPD curve estimates, such as the Salton Sea, California where salt concentrations are greater than sea water.

In contrast, detailed analysis of geothermal wells completed in two-phase liquid or vapour dominated geothermal resources benefits greatly from accounting for gas and salt concentrations. This is especially true where inflections in the temperature profile are used to support interpretation of well processes and the connection with the reservoir (feed zones).  

TBC - provide well example(s)

### References TBC

Canet, C., Franco, S. I., Prol-Ledesma, R. M., González-Partida, E., and Villanueva-Estrada, R. E., 2011, A model of boiling for fluid inclusion studies: Application to the Bolaños Ag–Au–Pb–Zn epithermal deposit, Western Mexico: Journal of Geochemical Exploration, v. 110, no. 2, p. 118-125.

Haas, J. L., Jr., 1971, The effect of salinity on the maximum thermal gradient of a hydrothermal system at hydrostatic pressure: Economic Geology and the Bulletin of the Society of Economic Geologists, v. 66, no. 6, p. 940-946.

Kaya, E., Zarrouk, S. J., and O'Sullivan, M. J., 2011, Reinjection in geothermal fields: A review of worldwide experience: Renewable and Sustainable Energy Reviews, v. 15, no. 1, p. 47-68.

Nicholson, K. N., 1993, Geothermal Fluids: Chemistry and Exploration Techniques, Berlin Heidelberg, Springer-Verlag.

White, D. E., Muffler, L. J. P., Truesdell, A. H., and Fournier, R. O., 1968, Preliminary results of research drilling in Yellowstone thermal areas: Transactions - American Geophysical Union, v. 49, no. 1, p. 358.

White, 1968, Hydrology, Activity, and Heat Flow of the Steamboat Springs Thermal System, Washoe County Nevada.

