## Initialize the workspace
rm(list=ls())
cat("\014")


library(circlize)

## --------------------------------
## Notes on using your own data to generate exposome globe
#' 1) correlation matrix
#' ¥ou can generate correlation matrix using cor() function, e.g., cor(mtcars).
#' If you have 2 samples and the data is stored in the dataframe "sample1" and "sample2", you can run cor(sample1, sample2)
#' "sample1" and "sample2" need to have the same number of variables
#' 
#' 2) name list
#' variables are grouped in the name list
#' You can define your own group name in the list, and it will be shown in the figure
#' However, the variable names in each group must be the same as those used in corrleion matrix
#' Names and variables are shown in the figure according to the order of the names and elements in the list
#' 
#' 3) color of the track
#' The number of colors in the "trackColorList" is the same as the number of the groups
#' You can generate a set of colors easily at https://medialab.github.io/iwanthue/
#' simply providing the number of colors wanted in the color palette and click "make a palette"
## --------------------------------



## LOAD function
source("~/Downloads/chord_function.R")

## LOAD data
load("~/Downloads/fakeData.Rdata")
#' The fake correlation data is using 2 samples as the published paper (male and female of a couple)
#' corrMatrixF is the correlation matrix of 128 chemicals within females
#' corrMatrixM is the correlation matrix of 128 chemicals within males
#' corrMatrixMF is the correlation matrix of 128 chemicals between females and males (note that diagnoal is not 1)


## Create list
#' All 128 chemicals are classified into groups and are specified in the list
#' Please note that the chemicals have a suffix of "_f" or "_m" to differentiate whether it is from females or males
#' We are using the correlation matrix for females only as the example.

clf = list(
    #PCBs
    PCBs = c(
        "PCB_28_f", "PCB_44_f", "PCB_49_f", "PCB_52_f", "PCB_66_f", "PCB_74_f", "PCB_87_f", "PCB_99_f", "PCB_101_f", "PCB_105_f", "PCB_110_f", "PCB_118_f", "PCB_114_f", "PCB_128_f", "PCB_138_f", "PCB_146_f", "PCB_149_f", "PCB_151_f", "PCB_153_f", "PCB_156_f", "PCB_157_f", "PCB_167_f", "PCB_170_f", "PCB_172_f", "PCB_177_f", "PCB_178_f", "PCB_180_f", "PCB_183_f", "PCB_187_f", "PCB_189_f", "PCB_194_f", "PCB_195_f", "PCB_196_f", "PCB_201_f", "PCB_206_f", "PCB_209_f" 
    ),
    #OCPs
    OCPs = c(
        "HCB_f", "b_HCB_f", "g_HCB_f", "op_DDT_f", "pp_DDE_f", "pp_DDT_f", "oxychlordane_f", "tr_nonachlor_f", "mirex_f"
    ),
    #PBC
    Polybrominated_cpds = c(
        "BB_153_f", "BDE_17_f", "BDE_28_f", "BDE_47_f", "BDE_66_f", "BDE_85_f", "BDE_99_f", "BDE_100_f", "BDE_153_f", "BDE_154_f", "BDE_183_f"
    ),
    #PFASs
    PFASs = c(
        "Et_PFOSA_AcOH_f", "Me_PFOSA_AcOH_f", "PFDeA_f", "PFNA_f", "PFOSA_f", "PFOS_f", "PFOA_f"
    ),
    #blood metals
    Blood_metals = c(
        "blood_Cd_f", "blood_Pb_f", "blood_Hg_f"
    ),
    #cotinine
    Cotinine = c(
        "cotinine_f"),
    #phytoestrogens
    Phytoestrogens = c(
        "genistein_f", "daidzein_f", "O_DMA_f", "equol_f", "enterodiol_f", "enterolactone_f"
    ),
    #phthalates
    Phthalates = c(
        "mMP_f", "mEP_f", "mCPP_f", "mBP_f", "miBP_f", "mECPP_f", "mCMHP_f", "mEHHP_f", "mEOHP_f", "mCHP_f", "mBzP_f", "mEHP_f", "mOP_f", "mNP_f"
    ),
    #phenols
    Phenols = c(
        "BPA_f", "2_OH_4MeO_BP_f", "4_OH_BP_f", "24_OH_BP_f", "22_OH_4MeO_BP_f", "2244_OH_BP_f"
    ),
    #anti microbial
    Anti_microbial_cpds = c(
        "MP_f", "EP_f", "PP_f", "BP_f", "BzP_f", "HP_f", "4_HB_f", "34_DHB_f", "OH_Me_P_f", "OH_Et_P_f", "TCS_f", "TCC_f"
    ),
    #paracetamol
    Paracetamols = c(
        "paracetamol_f", "4_aminophenol_f"
    ),
    #urine metals
    Urine_metals = c(
        "manganese_f", "chromium_f", "beryllium_f", "cobalt_f", "molybdenum_f", "cadmium_f", "tin_f", "caesium_f", "barium_f", "nickel_f", "copper_f", "zinc_f", "tungsten_f", "platinum_f", "thallium_f", "lead_f", "uranium_f"
    ),
    #urine metalloids
    Urine_metalloids = c(
        "selenium_f", "arsenic_f", "antimony_f", "tellurium_f"
    )
)


# color of the track
trackColorList <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                    "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1")
#' The color of the chemical group shown in the figure
#' Since we have 13 chemical classes in the example, you have to pick 13 colors



chord_1_sample(
    ## Inputs
    tracksample1 = clf, 
    correlationsample1 = corrMatrixF, 
    trackColor = trackColorList,
    # Track options
    trackFontSize = 0.65, # font size of the chemical class
    textDistance = 5, # distance between the chemical class to the track of chemicals
    # Correlation matrix options  
    corrLineColor = c("darkgreen", "white", "red"),  # The color of the correlation lines
    corrLineColorRampLv = c(-0.75, 0, 0.75), # What are the correlation values to build the color gradient based on the above
    corrNALv = 0.25, # To avoid plotting a hair ball, weaker correlations are sub. to NA.
    # Line options
    LineWidthMultipler = 2, # correlation line width is proportion to the correlation value. This is a multipler to enhance the contrast
    WithinCorrHeightMultiplier_inclass = -0.35, # the "within chemical class" correlation line height for "correlationsample1" and "correlationsample2" inputs, i.e., those lines drawn in the outter side of the tracks. Note that the heights are proportional to the correlation values. This is a multipler to enhance the contrast
    WithinCorrHeightMultiplier_betweenclass = 0.45, # similar to above but for the correlation lines connecting across chemical classes
    WithinCorrHeightOffset_betweenclass = 0.01) # A fixed value added to the proportional height to the above setting. Adjusting this setting could make it easier to observe the correlation patterns between chemical classes.



## For MacOS user:
# Save to pdf
quartz.save("~/Downloads/chord_1_sample.pdf", type = "pdf", width = 7, height = 7, device = dev.cur(), dpi = 400)

# Save to png
# quartz.save("~/Downloads/chord_1_sample.png", type = "png", width = 7, height = 7, device = dev.cur(), dpi = 400)


## For Windows user:
# dev.copy2pdf(width = 7, height = 7, file = "~/Downloads/chord_1_sample.pdf")


    
