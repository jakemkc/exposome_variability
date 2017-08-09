## Nov 28 2016
## Goal: "intra category" correlation

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L


# Load data (Total lipid and creatinine adjusted spearman r)
load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata")

# ******** -----
# A. Extract spearman coef ------------------------------------------------

# within female spearman matrix
indF <- which(upper.tri(corrubinf, diag = F) , arr.ind = T)  # arr.ind: return array indices

spearF <- data.frame(col = dimnames(corrubinf)[[2]][indF[,2]],  # picking the selected (with repeat) from the list, NOT extract
                     row = dimnames(corrubinf)[[1]][indF[,1]],
                     val = corrubinf[indF])  # Extracting the df with the indF matrix's coordinated pair into a vector

# Create pair ID
spearF$col2 <- gsub("_f$", "", spearF$col, ignore.case = T)

spearF$row2 <- gsub("_f$", "", spearF$row, ignore.case = T)

spearF$col_row_pair_2 <- sprintf("%s-%s", spearF$col2, spearF$row2)



# within male spearman matrix
indM <- which(upper.tri(corrubinm, diag = F) , arr.ind = TRUE )

spearM <- data.frame(col = dimnames(corrubinm)[[2]][indM[,2]],
                     row = dimnames(corrubinm)[[1]][indM[,1]],
                     val = corrubinm[indM])
# OK with manual inspection

# Create pair ID
spearM$col2 <- gsub("_m$", "", spearM$col, ignore.case = T)

spearM$row2 <- gsub("_m$", "", spearM$row, ignore.case = T)

spearM$col_row_pair_2 <- sprintf("%s-%s", spearM$col2, spearM$row2)



# within couple spearman matrix 
    ## full matrix
indMF <- which(corrubinmf==corrubinmf, arr.ind = TRUE)

spearMF <- data.frame(col = dimnames(corrubinmf)[[2]][indMF[,2]],
                     row = dimnames(corrubinmf)[[1]][indMF[,1]],
                     val = corrubinmf[indMF])
# OK with manual inspection

# Create pair ID 
spearMF$col2 <- gsub("_f$", "", spearMF$col, ignore.case = T)

spearMF$row2 <- gsub("_m$", "", spearMF$row, ignore.case = T)

spearMF$col_row_pair_2 <- sprintf("%s-%s", spearMF$col2, spearMF$row2)

# ******** -----
# B. Setup within chem category r df ------------------------------------------------

# category. 13 categories, 1 only with 1 chemicals (cotinine), 1 with only 2 (paracemtacl)
cl_f = list(
    #PCBs
    PCBs = c(
        "PCB_28_f", 
        "PCB_44_f", 
        "PCB_49_f", 
        "PCB_52_f", 
        "PCB_66_f", 
        "PCB_74_f", 
        "PCB_87_f", 
        "PCB_99_f", 
        "PCB_101_f", 
        "PCB_105_f", 
        "PCB_110_f", 
        "PCB_118_f", 
        "PCB_114_f", 
        "PCB_128_f", 
        "PCB_138_f", 
        "PCB_146_f", 
        "PCB_149_f", 
        "PCB_151_f", 
        "PCB_153_f", 
        "PCB_156_f", 
        "PCB_157_f", 
        "PCB_167_f", 
        "PCB_170_f", 
        "PCB_172_f", 
        "PCB_177_f", 
        "PCB_178_f", 
        "PCB_180_f", 
        "PCB_183_f", 
        "PCB_187_f", 
        "PCB_189_f", 
        "PCB_194_f", 
        "PCB_195_f", 
        "PCB_196_f", 
        "PCB_201_f", 
        "PCB_206_f", 
        "PCB_209_f" 
    ),
    #OCPs
    OCPs = c(
        "HCB_f",
        "b_HCB_f",
        "g_HCB_f",
        "op_DDT_f",
        "pp_DDE_f",
        "pp_DDT_f",
        "oxychlordane_f",
        "tr_nonachlor_f",
        "mirex_f"),
    #PBC
    Polybrominated_cpds = c(
        "BB_153_f", 
        "BDE_17_f", 
        "BDE_28_f", 
        "BDE_47_f", 
        "BDE_66_f", 
        "BDE_85_f", 
        "BDE_99_f", 
        "BDE_100_f", 
        "BDE_153_f", 
        "BDE_154_f", 
        "BDE_183_f"),
    #PFASs
    PFASs = c(
        "Et_PFOSA_AcOH_f", 
        "Me_PFOSA_AcOH_f", 
        "PFDeA_f", 
        "PFNA_f", 
        "PFOSA_f", 
        "PFOS_f", 
        "PFOA_f"),
    #blood metals
    Blood_metals = c(
        "blood_Cd_f", 
        "blood_Pb_f", 
        "blood_Hg_f"),
    #cotinine
    Cotinine = c(
        "cotinine_f"),
    #phytoestrogens
    Phytoestrogens = c(
        "genistein_f", 
        "daidzein_f", 
        "O_DMA_f", 
        "equol_f", 
        "enterodiol_f", 
        "enterolactone_f"),
    #phthalates
    Phthalates = c(
        "mMP_f", 
        "mEP_f", 
        "mCPP_f", 
        "mBP_f", 
        "miBP_f", 
        "mECPP_f", 
        "mCMHP_f", 
        "mEHHP_f", 
        "mEOHP_f", 
        "mCHP_f", 
        "mBzP_f", 
        "mEHP_f", 
        "mOP_f", 
        "mNP_f"),
    #phenols
    Phenols = c(
        "BPA_f", 
        "2_OH_4MeO_BP_f", 
        "4_OH_BP_f", 
        "24_OH_BP_f", 
        "22_OH_4MeO_BP_f", 
        "2244_OH_BP_f"),
    #anti microbial
    Anti_microbial_cpds = c(
        "MP_f", 
        "EP_f", 
        "PP_f", 
        "BP_f", 
        "BzP_f", 
        "HP_f", 
        "4_HB_f", 
        "34_DHB_f", 
        "OH_Me_P_f", 
        "OH_Et_P_f", 
        "TCS_f", 
        "TCC_f"),
    #paracetamol
    Paracetamols = c(
        "paracetamol_f", 
        "4_aminophenol_f"),
    #urine metals
    Urine_metals = c(
        "manganese_f", 
        "chromium_f", 
        "beryllium_f", 
        "cobalt_f", 
        "molybdenum_f", 
        "cadmium_f", 
        "tin_f", 
        "caesium_f", 
        "barium_f", 
        "nickel_f", 
        "copper_f", 
        "zinc_f", 
        "tungsten_f", 
        "platinum_f", 
        "thallium_f", 
        "lead_f", 
        "uranium_f"),
    #urine metalloids
    Urine_metalloids = c(
        "selenium_f", 
        "arsenic_f", 
        "antimony_f", 
        "tellurium_f")
)

cl_m <- cl_f

# replace the _f in the cl_m to _m
for(i in 1:length(cl_m)){
    cl_m[[i]] <- gsub("_f$", "_m", cl_m[[i]], ignore.case = T)
    # print(cl_m[[i]])
}

######################################################### #
# \\ B2. Setup boxplot----
######################################################### #

# Create categories

# female
spearF$colcateg <- "n____"
spearF$rowcateg <- "n____"


for(i in 1:length(cl_f)) {
    mtchC <- which(spearF$col %in% cl_f[[i]])
    spearF$colcateg[mtchC] <- names(cl_f)[i]
    mtchR <- which(spearF$row %in% cl_f[[i]])
    spearF$rowcateg[mtchR] <- names(cl_f)[i]
}

spearF$categlab <- "Heterogenous"
replistF <- spearF$colcateg == spearF$rowcateg 
spearF$categlab[replistF] <- spearF$colcateg[replistF]

spearF_boxp <- subset(spearF, c(!(categlab == "Heterogenous"))) # wo "heterogenous"
spearF_boxp_tmp <- spearF
spearF_boxp_tmp$categlab <- "All"
spearF_boxp <- rbind(spearF_boxp, spearF_boxp_tmp)
spearF_boxp$absR <- abs(spearF_boxp$val)



# male
spearM$colcateg <- "n____"
spearM$rowcateg <- "n____"


for(i in 1:length(cl_m)) {
    mtchC <- which(spearM$col %in% cl_m[[i]])
    spearM$colcateg[mtchC] <- names(cl_m)[i]
    mtchR <- which(spearM$row %in% cl_m[[i]])
    spearM$rowcateg[mtchR] <- names(cl_m)[i]
}

spearM$categlab <- "Heterogenous"
replistM <- spearM$colcateg == spearM$rowcateg 
spearM$categlab[replistM] <- spearM$colcateg[replistM]

# Misc setup
spearM_boxp <- subset(spearM, c(!(categlab == "Heterogenous")))
spearM_boxp_tmp <- spearM
spearM_boxp_tmp$categlab <- "All"
spearM_boxp <- rbind(spearM_boxp, spearM_boxp_tmp)
spearM_boxp$absR <- abs(spearM_boxp$val)


# couple
# male
spearMF$colcateg <- "n____"
spearMF$rowcateg <- "n____"



for(i in 1:length(cl_m)) {
    mtchC <- which(spearMF$col %in% cl_f[[i]])  #match column
    spearMF$colcateg[mtchC] <- names(cl_f)[i]
    mtchR <- which(spearMF$row %in% cl_m[[i]])  # match row
    spearMF$rowcateg[mtchR] <- names(cl_m)[i]
}

spearMF$categlab <- "Heterogenous"
spearMF$colcateg <- gsub("f$", "", spearMF$colcateg, ignore.case = T)
spearMF$rowcateg <- gsub("m$", "", spearMF$rowcateg, ignore.case = T)

replistMF <- spearMF$colcateg == spearMF$rowcateg 
spearMF$categlab[replistMF] <- spearMF$colcateg[replistMF]

# Misc setup
spearMF_boxp <- subset(spearMF, c(!(categlab == "Heterogenous")))
spearMF_boxp_tmp <- spearMF
spearMF_boxp_tmp$categlab <- "All"
spearMF_boxp <- rbind(spearMF_boxp, spearMF_boxp_tmp)
spearMF_boxp$absR <- abs(spearMF_boxp$val)


######################################################### #
# \\ B2. Boxplots----
######################################################### #

# ggplot, female
library(ggplot2)

# quartz(height=6, width=8)

p <- ggplot(data = spearF_boxp, aes(x = reorder(categlab, -absR, FUN = median), y = absR)) # -absR as revered order; absR as normal order
p <- p + geom_boxplot()
p <- p + ggtitle("Intra category corr (Spearman, abs r), Within Female") +
        xlab("Chemical classes") + 
        ylab("Absolute correlation")
p <- p + ylim(0, 1)
p <- p + geom_hline(yintercept = 0.08763547)
pf <- p + theme(axis.text.x  = element_text(angle=45, vjust=0.7, hjust= 0.7, size=10)) # move to the left a bit (using vjust since the labels are rotated)
pf

ggsave("results/boxplot_r_f_within_category_ggplot_v2.png", scale=1, dpi=400)

# IQRs
do.call("rbind", tapply(spearF_boxp$absR, spearF_boxp$categlab, quantile))

#                           0%        25%        50%       75%      100%
# All                 1.147983e-05 0.02373718 0.05358920 0.1136130 0.9664594
# Anti_microbial_cpds 2.488508e-03 0.06079075 0.10586855 0.1875692 0.8080839
# Blood_metals        8.428001e-02 0.15199999 0.21971998 0.2571962 0.2946725
# OCPs                7.894854e-03 0.14672326 0.33702398 0.4432452 0.8008421
# Paracetamols        1.137124e-01 0.11371239 0.11371239 0.1137124 0.1137124
# PCBs                3.769511e-03 0.18594190 0.37648044 0.5967019 0.9664594
# PFASs               2.308017e-02 0.07602142 0.17748399 0.5611388 0.8703481
# Phenols             1.060026e-02 0.09780848 0.19032106 0.2907887 0.8880725
# Phthalates          2.034649e-03 0.04777790 0.14760981 0.3376916 0.8426250
# Phytoestrogens       7.433051e-03 0.17173156 0.19540116 0.3962171 0.8635675
# Polybrominated_cpds 3.441057e-03 0.29335001 0.45413341 0.6835102 0.9074570
# Urine_metalloids    1.724230e-02 0.03540358 0.09154023 0.1722652 0.2109672
# Urine_metals        2.460949e-03 0.10204321 0.19094770 0.2902442 0.6267938

# ggplot, male
library(ggplot2)

# get the order from female plot
    # test1 <- reorder(spearF_boxp$categlab, -spearF_boxp$absR, FUN = median) # -spearF_boxp$absR: reverse median from high to low
    # dput(levels(test1))
female_lv <- c("Polybrominated_cpds", "PCBs", "OCPs", "Blood_metals", "Phytoestrogens", 
               "Urine_metals", "Phenols", "PFASs", "Phthalates", "Paracetamols", 
               "Anti_microbial_cpds", "Urine_metalloids", "All")

p <- ggplot(data = spearM_boxp, aes(x = factor(categlab, levels = female_lv), y = absR)) # -absR as revered order; absR as normal order
p <- p + geom_boxplot()
p <- p + ggtitle("Intra category corr (Spearman, abs r), Within male") +
    xlab("Chemical classes") + 
    ylab("Absolute correlation")
p <- p + ylim(0, 1)
p <- p + geom_hline(yintercept = 0.08762994)
pm <- p + theme(axis.text.x  = element_text(angle=45, vjust=0.7, hjust= 0.7, size=10)) # move to the left a bit (using vjust since the labels are rotated)
pm

ggsave("results/boxplot_r_m_within_category_ggplot_v2.png", scale=1, dpi=400)

# IQRs
do.call("rbind", tapply(spearM_boxp$absR, spearM_boxp$categlab, quantile))

#                           0%        25%        50%       75%      100%
# All                 9.313644e-06 0.02441996 0.05349058 0.1059788 0.9541171
# Anti_microbial_cpds 2.495408e-03 0.05441808 0.08931950 0.1632549 0.6620441
# Blood_metals        2.581070e-02 0.10143594 0.17706118 0.2488405 0.3206198
# OCPs                5.079724e-03 0.14697284 0.29084436 0.3738199 0.8233692
# Paracetamols        2.292490e-02 0.02292490 0.02292490 0.0229249 0.0229249
# PCBs                2.781649e-03 0.18905841 0.36916280 0.5855265 0.9541171
# PFASs               1.912452e-02 0.08308655 0.18551436 0.4259897 0.8377637
# Phenols             6.706766e-02 0.09241150 0.17177420 0.3259517 0.8336033
# Phthalates          1.051173e-03 0.03416967 0.12606872 0.2671517 0.8343612
# Phytoestrogens       1.184393e-02 0.14972916 0.17702124 0.3709395 0.8482781
# Polybrominated_cpds 3.629525e-02 0.29040887 0.41269203 0.6501801 0.9087436
# Urine_metalloids    2.225244e-02 0.04771969 0.07806818 0.1369303 0.2193285
# Urine_metals        1.342424e-03 0.08798792 0.16149378 0.2663486 0.5830269


# ggplot, couple
library(ggplot2)

female_lv_2 <- c("Polybrominated_cpds", "PCBs", "OCPs", "Blood_metals", "Phytoestrogens", 
               "Urine_metals", "Phenols", "PFASs", "Phthalates", "Paracetamols", 
               "Anti_microbial_cpds", "Urine_metalloids", "Cotinine", "All")

p <- ggplot(data = spearMF_boxp, aes(x = factor(categlab, levels = female_lv_2), y = absR)) # -absR as revered order; absR as normal order
p <- p + geom_boxplot()
p <- p + ggtitle("Intra category corr (Spearman, abs r), Within couple") +
    xlab("Chemical classes") + 
    ylab("Absolute correlation")
p <- p + ylim(0, 1)
p <- p + geom_hline(yintercept = 0.08761084)
pc <- p + theme(axis.text.x  = element_text(angle=45, vjust=0.7, hjust= 0.7, size=10)) # move to the left a bit (using vjust since the labels are rotated)
pc

ggsave("results/boxplot_r_mf_within_category_ggplot_v2.png", scale=1, dpi=400)

# IQRs
do.call("rbind", tapply(spearMF_boxp$absR, spearMF_boxp$categlab, quantile))


#                            0%        25%        50%        75%      100%
# All                 1.030608e-06 0.01908198 0.04042582 0.07319634 0.7286553
# Anti_microbial_cpds 9.486938e-04 0.02277612 0.05466717 0.10009212 0.4137254
# Blood_metals        2.251812e-02 0.14283550 0.19242320 0.43242201 0.5947251
# Cotinine            6.153042e-01 0.61530415 0.61530415 0.61530415 0.6153042
# OCPs                4.497030e-04 0.04323014 0.09178286 0.13881228 0.3824282
# Paracetamols        1.984245e-02 0.02885775 0.13216637 0.24178844 0.2697441
# PCBs                8.518740e-05 0.08802201 0.17061857 0.23845482 0.4980726
# PFASs               2.170265e-03 0.06723019 0.13297604 0.25115162 0.7286553
# Phenols             4.899969e-04 0.05402743 0.13499856 0.22468195 0.4842363
# Phthalates          1.707660e-04 0.02008516 0.04565449 0.09346579 0.3285960
# Phytoestrogens       1.611232e-03 0.01836758 0.06743420 0.19727749 0.2841075
# Polybrominated_cpds 1.002362e-04 0.12364484 0.21041181 0.34549444 0.5055719
# Urine_metalloids    7.873779e-03 0.02519231 0.04008211 0.08680822 0.2936487
# Urine_metals        1.030608e-06 0.01618411 0.03267931 0.05719062 0.4784089


