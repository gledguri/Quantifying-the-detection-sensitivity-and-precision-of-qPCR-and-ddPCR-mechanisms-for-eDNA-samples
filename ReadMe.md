# Quantifying the Detection Sensitivity and Precision of qPCR and ddPCR Mechanisms for eDNA Samples

## Overview
This repository supports the research article titled **"Quantifying the detection sensitivity and precision of qPCR and ddPCR mechanisms for eDNA samples."** The study investigates and compares the sensitivity and quantification precision of quantitative PCR (qPCR) and droplet digital PCR (ddPCR) in environmental DNA (eDNA) surveys, focusing on three teleost fish species assays.

## Research aim

The aim of this study is to empirically assess the sensitivity and precision of these methods through Bayesian analysis applied to known DNA concentrations (standards) and environmental samples. This helps determine the conditions under which ddPCR outperforms qPCR and provides insights into optimizing quantitative eDNA protocols.

### Key Findings
- **Higher Sensitivity of ddPCR**: ddPCR demonstrated greater sensitivity than qPCR at low DNA concentrations ($10^{-2}$ to $10^0$ copies/µL), typical of eDNA samples. This improved performance can be attributed to ddPCR’s end-point detection, which mitigates PCR inhibition effects and amplification efficiency variability.
- **Precision in Quantification**: ddPCR provided more precise estimates across all concentrations, typically by half to one order of magnitude. This precision is essential for studies where accurate quantification at low concentrations is critical.
- **Lower Limit of Quantification for ddPCR**: We propose an alternative method for estimating the lower limit of quantification (Clow-threshold, Clt) for ddPCR, defined by the presence of a single positive droplet in the generated droplet pool. This approach suggests that increasing the number of droplets analyzed can further reduce the Clt, enhancing detection and quantification limits.
- **Correlation Between Methods**: Despite the differences in sensitivity and precision, both methods showed a positive correlation in DNA quantification, aligning with previous research. However, ddPCR consistently outperformed qPCR in precision at low concentrations.
- **Assay-Specific Observations**: Differences were noted between the assays used, with the cod assay showing less variation between qPCR and ddPCR compared to the herring assay.
- **Role of Standard Samples**: Including standard samples during ddPCR optimization was found to be beneficial for understanding assay-specific behavior and improving quantification accuracy.

### Discussion and Recommendations
The study highlights the strengths and limitations of qPCR and ddPCR in the context of eDNA research. ddPCR’s higher sensitivity and precision at low concentrations make it advantageous for environmental monitoring where DNA is often present at trace levels. The findings suggest that while qPCR remains a useful tool, ddPCR should be the preferred method for low-concentration applications.

The inclusion of standard samples for calibration, even in ddPCR, provides additional accuracy and aids in troubleshooting assay performance. The study also recommends using technical replicates to improve precision and reduce detection uncertainty, especially in low-concentration scenarios.

### Conclusion
This research provides critical insights into the comparative performance of qPCR and ddPCR, with empirical evidence supporting the use of ddPCR for more sensitive and precise eDNA quantification at low concentrations. The study recommends ddPCR for eDNA concentrations below 1 copy/µL and underscores the value of including standard samples during optimization. The two-step qPCR model introduced here can enhance quantification precision and guide future eDNA survey methodologies.

## Citation
If you use any part of this repository in your work, please cite:
**Guri, G., Ray, J.L., Shelton, A.O., Kelly, R.P., Præbel, K., et al. (2024). Quantifying the detection sensitivity and precision of qPCR and ddPCR mechanisms for eDNA samples.**

## Acknowledgments
We acknowledge Agneta Hansen for significant contributions to field and lab work. This research was funded by the project FISHDIV (NRC301691) funded by the Norwegian Research Council (NRC).

## Contact
For questions or further information, please contact:
**Gledis Guri**  
Email: gled.guri@gmail.com  
Phone: +46 700658104  

---

**Keywords**: eDNA, qPCR, ddPCR, sensitivity, quantification precision, fish monitoring.
