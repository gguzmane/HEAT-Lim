# ASU/USYD HEAT-Lim (Human/Environmental Adaptation and Threshold Limits Model) 

**Physiological-based estimations of human survivability and liveability under extreme heat in a changing climate.**

This repository hosts the Python implementation of the ASU/USYD Human/Environmental Adaptation and Threshold Limits Model (HEAT-Lim), formerly known as PyHHB. HEAT-Lim supports assessments of heat survivability and heat liveability by simulating human-environment heat exchange using a biophysical model. It is designed for use with climate model output, following the approach developed in Vanos et al. (2023). 

If you use this code, we kindly ask that you cite both the study and the code repository as follows:

**First mention in your work:**  ASU/USYD Human/Environmental Adaptation and Threshold Limits Model (HEAT-Lim), and cite the original 2023 study (https://doi.org/10.1038/s41467-023-43121-5) and the zenodo repository of the code (https://zenodo.org/records/10020137).
**Subsequent mentions**:  Simply refer to the model as HEAT-Lim! thanks! :)

**UPDATE NOTICE:**
Originally, HEAT-Lim’s liveability analyses we applied the method limited to air temperatures above 25°C, to avoid cases where the maximum metabolic heat production (`Mmax`) becomes unrealistically high under normothermic or cold conditions. We have since updated the model to cap `Mmax` at 10 METs, which reflects a reasonable upper limit of sustainable activity intensity for a healthy adult under thermal comfort or cold stress.

# What does this repository contains?
This code supports replication of the:

- HEAT-Lim Survivability matrices (Figures 2, S2, S3 in the paper)
- HEAT-Lim Liveability matrices (Figures 4, S6, S7)

...based on the environmental and personal characteristics evaluated in the study.

In Vanos et al. (2023), HEAT-Lim builds on the wet-bulb temperature threshold of 35°C as a limit for human heat-stroke survivability (Sherwood & Huber, 2010). The model expands upon this by providing a more flexible and physiologically grounded approach that:

- Applies to diverse climate regimes
- Can be customized for population groups with potential co-morbidities
- Integrates well-established and fundamental principles from thermal physiology and human biophysics
- Aligns with 3-hour and 6-hour exposure windows from climate model outputs and previous studies

For detailed information about the differences between HEAT-Lim and the wet bulb temperature limit of 35°C, read the session *"Considerations of New Model Estimating Physiological Survivability Limits and Liveability"* in the paper's Supplemental material.

## How is defined the survivability limit?

The limit of survivability to heat stroke death is determined by detecting the conditions in which a person would reach a core temperature of 43°C in 3- or 6-hour exposure windows to allow for comparison with the Tw of 35°C assumption (heat stroke death after 6 hours) used in previous survivability assessments based on Sherwood and Huber, (2010). 

## How is the maximum internal heat production (Mmax) defined as a liveability metric?

Mmax is the maximum safe internal heat production, or level of physical activity, that a person can generate without a sustained positive rate of heat storage in the prevailing environment, thus allowing safe, sustained work and play for an extended period.

## What does this repository folders contain?

- *Ancillary:* This folder contains necessary supporting images for the Mmax scale illustrating different activity levels based on metabolic equivalent units (Ainsworth et al., 2011). It also includes the matrices with pre-calculated Tw using the Davies-Jones method (Davies-Jones, 2008) to allow direct graphic comparison with the isothermal Tw=35°C in the tutorial for survivability matrix estimation.
- *Codes:*  This folder contains a Python script called HEATLim.py to be invoked in further routines as the HEAT-Lim module (import HEATLim as HEATLim), and two Jupyter notebooks with tutorials to replicate the survivability (Figures 2, S2, S3) and liveability (Figure 4, S6, S7) matrices for the environmental and personal features included in the study. The equations to apply the human-environment exchange model in HHB.py follow the notation and rationale from Vanos et al. (2023) and Cramer and Jay (2019).
- *Outputs:* This folder contains the output files generated after running the tutorials.
- *Personal profiles:*  This folder contains text files that include the information required for setting up anthropometrics, activities, and clothing for a target population.


**NOTE: *If you want to run this on your computer, we recommend running the survivability tutorial first once the survivability outputs are needed for liveability analysis*.** 

**Bibliography:**

- Ainsworth, B. E., W. L. Haskell, S. D. Herrmann, N. Meckes, D. R. Bassett, C. Tudor-Locke, J. L. Greer, J. Vezina, M. C. Whitt-Glover, and A. S. Leon. 2011. 2011 compendium of physical activities: A second update of codes and MET values. Medicine and Science in Sports and Exercise, 43(8):1575–1581, https://doi.org/10.1249/MSS.0b013e31821ece12.
Cramer, M. N., and O. Jay. 2019. Cores of reproducibility in physiology partitional calorimetry. Journal of Applied Physiology, 126(2):267–277, https://doi.org/10.1152/japplphysiol.00191.2018.
- Cramer, M. N., & Jay, O. (2019). Cores of reproducibility in physiology partitional calorimetry. Journal of Applied Physiology, 126(2), 267–277. https://doi.org/10.1152/japplphysiol.00191.2018
- Davies-Jones, R. 2008. An Efficient and Accurate Method for Computing the Wet-Bulb Temperature along Pseudoadiabats. Monthly Weather Review, 136(7):2764–2785, https://doi.org/https://doi.org/10.1175/2007MWR2224.1.
-Sherwood, S. C., and M. Huber. 2010. An adaptability limit to climate change due to heat stress. Proceedings of the National Academy of Sciences of the United States of America, 107(21):9552–9555, https://doi.org/10.1073/pnas.0913352107.
- Vanos, J., Guzman-Echavarria, G., Baldwin, J. W., Bongers, C., Ebi, K. L., & Jay, O. (2023). A physiological approach for assessing human survivability and liveability to heat in a changing climate. Nat Commun 14, 7653 (2023). https://doi.org/10.1038/s41467-023-43121-5
