# Simulation_SEmismatch

# Background
Striped bass *Morone saxatilis* have undergone many dynamic changes as a result of changing environmental conditions, particularly in the Chesapeake Bay. For example, water temperature has been shown to affect striped bass spawning behavior, recruitment, movement, and natural mortality. We developed a simulation framework to evaluate how potential climate-change induced variability in the Chesapeake Bay striped bass stock may affect the performance of a spatially-explicit stock assessment model. We simulated five scenarios with alternative assumptions of potential climate effects on striped bass. 

This repository has the final versions of five operating models (developed in R studio) and 3 estimation model code (.tpl) and dat file (.dat) for  Atlantic Striped Bass that were developed by Samara Nehemiah during her Ph.D. research at the Chesapeake Biological Laboratory. All estimation models were developed in ADMB (https://www.admb-project.org/). 

# Files included

**Operating models**
This study considered five different operating models that included alternative assumptions of potential climate change effects on striped bass. 

*OM-B* was the base scenario. In the base scenario, population dynamics for natural mortality and movement were constant over time and recruitment followed a stationary stochastic process. OM-B followed the same dynamics described in the simulation (https://github.com/snehemiah/SE-model-performance) and used parameter estimates from the original spatially-explicit model developed for striped bass (https://github.com/snehemiah/StripedBass-spatiallyexplicit-models). <br />
*OM-M*   assumed that natural mortality changed over time based on the results of the Brownie mark-recapture model that estimated natural mortality from 1990-2020 (Schonfeld, 2023).<br />
*OM-P* assumed migration of striped bass changed over time. OM-P assumed migration into and out of the Chesapeake Bay occurred more gradually throughout the year than at 6-month intervals to better match the observed movement of Chesapeake Bay striped bass (Secor et al., 2020).<br />
*OM-R* represented changes to recruitment and productivity over time. Average recruitment was simulated to decline over time for the Chesapeake Bay stock.<br />
*OM-MPR* assumed the same changes in mortality as OM-M, occupancy as OM-p, and recruitment as OM-R.<br />

Full descriptions of the operating models can be found in Nehemiah (2024). 

**Estimation models**
Three estimation models were evaluated for each data generating scenario. The EM folder contains all .dat and .tpl files for the estimation models in ADMB. 

*FAA* was a similar structure to the current stock assessment model used to inform management of striped bass (Northeast Fisheries Science Center, 2019). This model was a spatially-implicit statical catch-at-age model that assumes fleets-as-areas. <br />
*SE* was a multi-stock, spatially-explicit population model that assumed two regions, two stocks, and had two 6-month time-steps.<br />
*SE-S* was a multi-stock, spatially-explicit model that included data with additional stock composition for the last ten years.<br />

The "scaa-stripedbass-faa" files correspond to the FAA models; the "scaa-stripedbass-se" correspond to the SE model; the "scaa-stripedbass-ses" correspond to the SE-S model with added stock data.

Full descriptions of the operating models can be found in Nehemiah (2024). 


# Partner Institutions
This work was funded by the NOAA Chesapeake Bay Office and the NMFS/Sea Grant Population and Ecosystems Dynamic Fellowship. Collaborators include researchers from the Chesapeake Biological Laboratory, Virginia Institute of Marine Science, and the NMFS Southeast Fisheries Science Cener. 

# Manuscript Citation
Nehemiah, S., R. Latour, A. Schonfeld, A. Schueller, and M.J. Wilberg. Testing the performance of a spatially-explicit population model for striped bass under different climate change scenarios. In preparation.

# Contact
If you have any questions, please contact Samara (she/her). Current email: snehemiah@asmfc.org. 

# References
Northeast Fisheries Science Center, 2019. 66th Northeast Regional Stock Assessment Workshop (66th SAW) Assessment Report. Northeast Fisheries Science Center Reference Document 457–1170.
Nehemiah, S. 2024.DEVELOPMENT AND EVALUATION OF SPATIALLY-EXPLICIT POPULATION MODELS FOR ESTIMATING THE ABUNDANCE OF CHESAPEAKE BAY FISHES. [http://hdl.handle.net/1903/33663.](https://doi.org/10.13016/ct4p-ym7b)
Schonfeld, A.J., 2023. Climate Impacts on Spatiotemporal Habitat Usage of Mid-Atlantic Fishes. William & Mary Ph.D. Dissertation.
Secor, D.H., O’Brien, M.H.P., Gahagan, B.I., Carter Watterson, J., Fox, D.A., 2020. Differential migration in Chesapeake Bay striped bass. PLoS One 15, 1–19. https://doi.org/10.1371/journal.pone.0233103

<img src="https://www.umces.edu/sites/default/files/UMCES-CBL-logo.jpg" jsaction="" class="sFlh5c pT0Scc iPVvYb" style="max-width: 600px; height: 221px; margin: 0px; width: 557px;" alt="UMCES CBL logo.jpg | University of Maryland Center for Environmental Science" jsname="kn3ccd" aria-hidden="false">
