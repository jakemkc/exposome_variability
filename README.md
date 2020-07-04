# Toward Capturing the Exposome: Exposure Biomarker Variability and Co-Exposure Patterns in the Shared Environment


## July 4, 2020 update:
I have created simple R wrapper functions for plotting the exposome globe. Now you can adjust commonly used parameters without going into the coding detail.
All materials are in the directory "R_function".

You can provide one sample or two samples (e.g., twins, couples) as the input.
<img src="https://github.com/jakemkc/exposome_variability/blob/master/R_function/chord_2_samples.png">

### Instruction
For MacOS:
1. Download the fake correlation matrix file (fakeData.Rdata), which contains two samples as used in the paper, and save it into the "Downloads" folder
2. Download the function script (chord_function.R), and save it into the "Downloads" folder
3. Download the example script, and run in Rstudio (example_1_sample.R or example_2_samples.R)
4. You can see the globe figure in Rstudio, and a pdf copy will also be saved to your "Downloads" folder
5. You can edit the pdf in any vector graphic editor (e.g., Inkscape) to fine tune different aspects of the figure

For other OSes:
- The steps are basically the same as above, but you have to modify the paths for 1) and 2) in the example script. Also, you need to run "dev.copy2pdf" instead of "quartz.save" to export the figure. I have included the saving command in the example script too.

### Go further
You can adjust the settings in the wrapper functions to see how they will affect the look of the globe. Alternatively, you can check the circlize package ebook provided by the author (https://jokergoo.github.io/circlize_book/book/)

### Citation
Please cite:
Chung, M. K., Kannan, K., Louis, G. M. & Patel, C. J. Toward capturing the exposome: exposure biomarker variability and coexposure patterns in the shared environment. Environmental Science & Technology 52, 8801–8810 (2018). http://pubs.acs.org/doi/10.1021/acs.est.8b01467

1.Gu, Z., Gu, L., Eils, R., Schlesner, M. & Brors, B. circlize Implements and enhances circular visualization in R. Bioinformatics 30, 2811–2812 (2014). https://doi.org/10.1093/bioinformatics/btu393.


------ 
## Authors
- Ming Kei (Jake) Chung
  - github: [\@jakemkc](http://github.com/jakemkc)
  - twitter: [\@jakekei](http://twitter.com/jakekei)
  - email: jake_chung[at]hms[dot]harvard[dot]edu
- Kurunthachalam Kannan
  - email: kurunthachalam[dot]kannan[at]health[dot]ny[dot]gov
- Germaine M. Buck Louis
  - email: louisg[at]mail[dot]nih[dot]gov
- Chirag J. Patel
  - github: [\@chiragjp](http://github.com/chiragjp)
  - web: [www.chiragjpgroup.org](http://www.chiragjpgroup.org)
  
## Abstract
### *BACKGROUND*
Along with time, variation in the exposome is dependent on the location and sex of study participants. One specific factor that may influence exposure co-variations is a shared household environment.

### *OBJECTIVES*
To examine the influence of shared household and partner’s sex in relation to the variation in 128 endocrine disrupting chemical (EDC) exposures among couples.

### *METHODS* 
In a cohort comprising 501 couples trying for pregnancy, we measured 128 (13 chemical classes) persistent and non persistent EDCs and estimated 1) sex-specific differences; 2) variance explained by shared household; and 3) Spearman's rank correlation coefficients (rs) for females, males, and couples.

## Getting Started
- This github repository analytical materials for characterzing the correlations in 501 subjects recruited in the LIFE Study.
- Code for data clean-up, processing, analysis, and visulization can be found on [GitHub](https://github.com/jakemkc/exposome_variability)
- Main findings of our study can be found [here](results/results.md)
- Preprint paper can be downloaded from [BioRxiv](https://doi.org/10.1101/175513)
- Please contact us for further information, thank you.

