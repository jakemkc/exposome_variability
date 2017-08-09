## Dec 4 2016
## Goal: Meff (effective no of variables)

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L

# Load
load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata") 


# ******** -----
# B. On imputed DF  -------------------------------------------

######################################################### #
# \\ B1. Preparing df -------------------------------------------------------

## ++ B1.1 for females ----

## These 133 name are not R fault-proof (e.g.4_HB_f). Order of (papername) content is matched with (name(papername))
papername_f <- c("pcb101amt", "pcb105amt", "pcb110amt", "pcb114amt", "pcb118amt", "pcb128amt", "pcb138amt", "pcb146amt", "pcb149amt", "pcb151amt", "pcb153amt", "pcb156amt", 
                 "pcb157amt", "pcb167amt", "pcb170amt", "pcb172amt", "pcb177amt", "pcb178amt", "pcb180amt", "pcb183amt", "pcb187amt", "pcb189amt", "pcb194amt", "pcb195amt", 
                 "pcb196amt", "pcb201amt", "pcb206amt", "pcb209amt", "pcb028amt", "pcb044amt", "pcb049amt", "pcb052amt", "pcb066amt", "pcb074amt", "pcb087amt", "pcb099amt", 
                 "pfcepahamt", "pfcmpahamt", "pfcpfdeamt", "pfcpfnaamt", "pfcpfsaamt", "pfcpfosamt", "pfcpfoaamt", "metbcdamt", "metbpbamt", "metthgamt", "popbb1amt", "popbhcamt", 
                 "popghcamt", "pophcbamt", "popmiramt", "popodtamt", "popoxyamt", "poppdeamt", "poppdtamt", "poptnaamt", "pbdebr1amt", "pbdebr2amt", "pbdebr3amt", "pbdebr4amt", 
                 "pbdebr5amt", "pbdebr6amt", "pbdebr7amt", "pbdebr8amt", "pbdebr9amt", "pbdebr66amt", 
                 "DAZAMOUNT", "DMAAMOUNT", "EQUAMOUNT", "ETDAMOUNT", "ETLAMOUNT", "GNSAMOUNT", "CREAMOUNT", "fcholamt", "cholamt", "trigamt", "phosamt", "cotamt", "Selenium", 
                 "Arsenic", "Manganese", "Chromium", "Beryllium", "Cobalt", "Molybdenum", "Cadmium_Corrected", "Tin", "Antimony", "Tellurium", "Caesium", "Barium", "Nickel", 
                 "Copper", "Zinc", "Tungsten", "Platinum", "Thallium", "Lead", "Uranium", "mMethylPhthalate", "mEthylPhthalate", "mCarboxyPropylPhthalate", "mButylPhthalate", 
                 "mIsobutylPhthalate", "mCarboxyEthylPentylPhthalate", "mCarboxyMethylHexylPhthalate", "mEthylHydroxyHexylPhthalate", "mEthylOxoHexylPhthalate", "mCycloHexylPhthalate", 
                 "mBenzylPhthalate", "mEthylHexylPhthalate", "mOctylPhthalate", "mIsononylPhthalate", "BPA", "HydroxyMethoxyBenzoPhenone", "HydroxyBenzoPhenone", "DiHydroxyBenzoPhenone", 
                 "DiHydroxyMethoxyBenzoPhenone", "TetraHydroxyBenzoPhenone", "MeP", "EtP", "PrP", "BuP", "BzP", "HeP", "X_4_HB", "X_3_4_DHB", "OH_MeP", "OH_EtP", "TCS", "TCC", "PAP", "APAP")


names(papername_f) <- c("PCB_101_f", "PCB_105_f", "PCB_110_f", "PCB_118_f", "PCB_114_f", "PCB_128_f", "PCB_138_f", "PCB_146_f", "PCB_149_f", "PCB_151_f", "PCB_153_f", "PCB_156_f", "PCB_157_f", "PCB_167_f", "PCB_170_f", "PCB_172_f", "PCB_177_f", "PCB_178_f", "PCB_180_f", "PCB_183_f", "PCB_187_f", "PCB_189_f", "PCB_194_f", "PCB_195_f", "PCB_196_f", "PCB_201_f", "PCB_206_f", "PCB_209_f", "PCB_28_f", "PCB_44_f", "PCB_49_f", "PCB_52_f", "PCB_66_f", "PCB_74_f", "PCB_87_f", "PCB_99_f", "Et_PFOSA_AcOH_f", "Me_PFOSA_AcOH_f", "PFDeA_f", "PFNA_f", "PFOSA_f", "PFOS_f", "PFOA_f", "blood_Cd_f", "blood_Pb_f", "blood_Hg_f", "BB_153_f", "b_HCB_f", "g_HCB_f", "HCB_f", "mirex_f", "op_DDT_f", "oxychlordane_f", "pp_DDE_f", "pp_DDT_f", "tr_nonachlor_f", "BDE_17_f", "BDE_28_f", "BDE_47_f", "BDE_85_f", "BDE_99_f", "BDE_100_f", "BDE_153_f", "BDE_154_f", "BDE_183_f", "BDE_66_f", "daidzein_f", "O_DMA_f", "equol_f", "enterodiol_f", "enterolactone_f", "genistein_f", "CREAMOUNT_f", "fcholamt_f", "cholamt_f", "trigamt_f", "phosamt_f", "cotinine_f", "selenium_f", "arsenic_f", "manganese_f", "chromium_f", "beryllium_f", "cobalt_f", "molybdenum_f", "cadmium_f", "tin_f", "antimony_f", "tellurium_f", "caesium_f", "barium_f", "nickel_f", "copper_f", "zinc_f", "tungsten_f", "platinum_f", "thallium_f", "lead_f", "uranium_f", "mMP_f", "mEP_f", "mCPP_f", "mBP_f", "miBP_f", "mECPP_f", "mCMHP_f", "mEHHP_f", "mEOHP_f", "mCHP_f", "mBzP_f", "mEHP_f", "mOP_f", "mNP_f", "BPA_f", "2_OH_4MeO_BP_f", "4_OH_BP_f", "24_OH_BP_f", "22_OH_4MeO_BP_f", "2244_OH_BP_f", "MP_f", "EP_f", "PP_f", "BP_f", "BzP_f", "HP_f", "4_HB_f", "34_DHB_f", "OH_Me_P_f", "OH_Et_P_f", "TCS_f", "TCC_f", "paracetamol_f", "4_aminophenol_f")

# remove _f
no_f <- gsub("_f$", "", colnames(resid_impute_f)[-1], ignore.case = T)

# match and get the index order by female_no_f
index_f <- vector(mode = "numeric")
for(i in 1:length(no_f)) {
    mtchF <- which(papername_f %in% no_f[i])
    print(mtchF) # error check
    index_f <- append(index_f, mtchF)
}

# rename residual df with my short-form
colnames(resid_impute_f)[-1] <- names(papername_f[index_f])

## List in right order (13 classes, 128 chemicals wo lipids and creatinine)
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


# Special: drop cotinine
cl_f$Cotinine <- NULL


## ++ B1.1 for females ----

## These 133 name are not R fault-proof (e.g.4_HB_f). Order of (papername) content is matched with (name(papername))
papername_m <- c("pcb101amt", "pcb105amt", "pcb110amt", "pcb114amt", "pcb118amt", "pcb128amt", "pcb138amt", "pcb146amt", "pcb149amt", "pcb151amt", "pcb153amt", "pcb156amt", 
                 "pcb157amt", "pcb167amt", "pcb170amt", "pcb172amt", "pcb177amt", "pcb178amt", "pcb180amt", "pcb183amt", "pcb187amt", "pcb189amt", "pcb194amt", "pcb195amt", 
                 "pcb196amt", "pcb201amt", "pcb206amt", "pcb209amt", "pcb028amt", "pcb044amt", "pcb049amt", "pcb052amt", "pcb066amt", "pcb074amt", "pcb087amt", "pcb099amt", 
                 "pfcepahamt", "pfcmpahamt", "pfcpfdeamt", "pfcpfnaamt", "pfcpfsaamt", "pfcpfosamt", "pfcpfoaamt", "metbcdamt", "metbpbamt", "metthgamt", "popbb1amt", "popbhcamt", 
                 "popghcamt", "pophcbamt", "popmiramt", "popodtamt", "popoxyamt", "poppdeamt", "poppdtamt", "poptnaamt", "pbdebr1amt", "pbdebr2amt", "pbdebr3amt", "pbdebr4amt", 
                 "pbdebr5amt", "pbdebr6amt", "pbdebr7amt", "pbdebr8amt", "pbdebr9amt", "pbdebr66amt", 
                 "DAZAMOUNT", "DMAAMOUNT", "EQUAMOUNT", "ETDAMOUNT", "ETLAMOUNT", "GNSAMOUNT", "CREAMOUNT", "fcholamt", "cholamt", "trigamt", "phosamt", "cotamt", "Selenium", 
                 "Arsenic", "Manganese", "Chromium", "Beryllium", "Cobalt", "Molybdenum", "Cadmium_Corrected", "Tin", "Antimony", "Tellurium", "Caesium", "Barium", "Nickel", 
                 "Copper", "Zinc", "Tungsten", "Platinum", "Thallium", "Lead", "Uranium", "mMethylPhthalate", "mEthylPhthalate", "mCarboxyPropylPhthalate", "mButylPhthalate", 
                 "mIsobutylPhthalate", "mCarboxyEthylPentylPhthalate", "mCarboxyMethylHexylPhthalate", "mEthylHydroxyHexylPhthalate", "mEthylOxoHexylPhthalate", "mCycloHexylPhthalate", 
                 "mBenzylPhthalate", "mEthylHexylPhthalate", "mOctylPhthalate", "mIsononylPhthalate", "BPA", "HydroxyMethoxyBenzoPhenone", "HydroxyBenzoPhenone", "DiHydroxyBenzoPhenone", 
                 "DiHydroxyMethoxyBenzoPhenone", "TetraHydroxyBenzoPhenone", "MeP", "EtP", "PrP", "BuP", "BzP", "HeP", "X_4_HB", "X_3_4_DHB", "OH_MeP", "OH_EtP", "TCS", "TCC", "PAP", "APAP")


names(papername_m) <- c("PCB_101_m", "PCB_105_m", "PCB_110_m", "PCB_118_m", "PCB_114_m", "PCB_128_m", "PCB_138_m", "PCB_146_m", "PCB_149_m", "PCB_151_m", "PCB_153_m", "PCB_156_m", "PCB_157_m", "PCB_167_m", "PCB_170_m", "PCB_172_m", "PCB_177_m", "PCB_178_m", "PCB_180_m", "PCB_183_m", "PCB_187_m", "PCB_189_m", "PCB_194_m", "PCB_195_m", "PCB_196_m", "PCB_201_m", "PCB_206_m", "PCB_209_m", "PCB_28_m", "PCB_44_m", "PCB_49_m", "PCB_52_m", "PCB_66_m", "PCB_74_m", "PCB_87_m", "PCB_99_m", "Et_PFOSA_AcOH_m", "Me_PFOSA_AcOH_m", "PFDeA_m", "PFNA_m", "PFOSA_m", "PFOS_m", "PFOA_m", "blood_Cd_m", "blood_Pb_m", "blood_Hg_m", "BB_153_m", "b_HCB_m", "g_HCB_m", "HCB_m", "mirex_m", "op_DDT_m", "oxychlordane_m", "pp_DDE_m", "pp_DDT_m", "tr_nonachlor_m", "BDE_17_m", "BDE_28_m", "BDE_47_m", "BDE_85_m", "BDE_99_m", "BDE_100_m", "BDE_153_m", "BDE_154_m", "BDE_183_m", "BDE_66_m", "daidzein_m", "O_DMA_m", "equol_m", "enterodiol_m", "enterolactone_m", "genistein_m", "CREAMOUNT_m", "fcholamt_m", "cholamt_m", "trigamt_m", "phosamt_m", "cotinine_m", "selenium_m", "arsenic_m", "manganese_m", "chromium_m", "beryllium_m", "cobalt_m", "molybdenum_m", "cadmium_m", "tin_m", "antimony_m", "tellurium_m", "caesium_m", "barium_m", "nickel_m", "copper_m", "zinc_m", "tungsten_m", "platinum_m", "thallium_m", "lead_m", "uranium_m", "mMP_m", "mEP_m", "mCPP_m", "mBP_m", "miBP_m", "mECPP_m", "mCMHP_m", "mEHHP_m", "mEOHP_m", "mCHP_m", "mBzP_m", "mEHP_m", "mOP_m", "mNP_m", "BPA_m", "2_OH_4MeO_BP_m", "4_OH_BP_m", "24_OH_BP_m", "22_OH_4MeO_BP_m", "2244_OH_BP_m", "MP_m", "EP_m", "PP_m", "BP_m", "BzP_m", "HP_m", "4_HB_m", "34_DHB_m", "OH_Me_P_m", "OH_Et_P_m", "TCS_m", "TCC_m", "paracetamol_m", "4_aminophenol_m")

# remove _m
no_m <- gsub("_m$", "", colnames(resid_impute_m)[-1], ignore.case = T)

# match and get the index order by female_no_f
index_m <- vector(mode = "numeric")
for(i in 1:length(no_m)) {
    mtchM <- which(papername_m %in% no_m[i])
    print(mtchM) # error check
    index_m <- append(index_m, mtchM)
}

# rename residual df with my short-form
colnames(resid_impute_m)[-1] <- names(papername_m[index_m])

## List in right order (13 classes, 128 chemicals wo lipids and creatinine)
cl_m <- cl_f

# replace the _f in the cl_m to _m
for(i in 1:length(cl_m)){
    cl_m[[i]] <- gsub("_f$", "_m", cl_m[[i]], ignore.case = T)
    # print(cl_m[[i]])
}


# Special: drop cotinine
cl_m$Cotinine <- NULL



######################################################### #
# \\ B2. Meff -------------------------------------------------------
######################################################### #

######################################################### #
## ++ B2.1 for females ----

library(psych)

Meff_f <- data.frame() 

for(i in 1:length(cl_f)){  
meff_list_f <- vector()
    for(j in 1:10) {
        tmpf <- subset(resid_impute_f, X_Imputation_ == j)[2:ncol(resid_impute_f)]
        cor_r <- as.matrix(corr.test(tmpf[, cl_f[[i]]], method="spearman")$r)  # corr matric of selected class
        eigen_f <- eigen(cor_r)$values # extract eigenvalue
        meff_list_f[j] <- 1+(ncol(cor_r)-1)*(1-var(eigen_f)/ncol(cor_r))  # get Meff per imputed frame (x10)
    }
frm_row_f <- as.data.frame(mean(meff_list_f)) #data.frame only exist in each i lv loop
frm_row_f$orig_n <- length(cl_f[[i]])
frm_row_f$class <- names(cl_f[i])
Meff_f <- rbind(Meff_f, frm_row_f)
}

Meff_f <- rbind(Meff_f, sapply(Meff_f[, 1:2], sum))
Meff_f[13,3] <- "sum"

# > Meff_f
# mean(meff_list_f) orig_n               class
# 1          28.381161     36                PCBs
# 2           7.943579      9                OCPs
# 3           8.222635     11 Polybrominated_cpds
# 4           6.134291      7               PFASs
# 5           2.905104      3        Blood_metals
# 6           5.341413      6       Phytoestrogens
# 7          12.719671     14          Phthalates
# 8           5.478211      6             Phenols
# 9          11.381480     12 Anti_microbial_cpds
# 10          1.986892      2        Paracetamols
# 11         16.128138     17        Urine_metals
# 12          3.947973      4    Urine_metalloids
# 13        110.570548    127                 sum


######################################################### #
## ++ B2.2 for males ----

library(psych)

Meff_m <- data.frame() 

for(i in 1:length(cl_m)){  
    meff_list_m <- vector()
    for(j in 1:10) {
        tmpm <- subset(resid_impute_m, X_Imputation_ == j)[2:ncol(resid_impute_m)]
        cor_r <- as.matrix(corr.test(tmpm[, cl_m[[i]]], method="spearman")$r)  # corr matric of selected class
        eigen_m <- eigen(cor_r)$values # extract eigenvalue
        meff_list_m[j] <- 1+(ncol(cor_r)-1)*(1-var(eigen_m)/ncol(cor_r))  # get Meff per imputed frame (x10)
    }
    frm_row_m <- as.data.frame(mean(meff_list_m)) #data.frame only exist in each i lv loop
    frm_row_m$orig_n <- length(cl_m[[i]])
    frm_row_m$class <- names(cl_m[i])
    Meff_m <- rbind(Meff_m, frm_row_m)
}

Meff_m <- rbind(Meff_m, sapply(Meff_m[, 1:2], sum))
Meff_m[13,3] <- "sum"


# > Meff_m
# mean(meff_list_m) orig_n               class
# 1          28.691510     36                PCBs
# 2           8.170960      9                OCPs
# 3           8.382288     11 Polybrominated_cpds
# 4           6.253780      7               PFASs
# 5           2.910057      3        Blood_metals
# 6           5.379011      6       Phytoestrogens
# 7          12.868670     14          Phthalates
# 8           5.509800      6             Phenols
# 9          11.603682     12 Anti_microbial_cpds
# 10          1.999182      2        Paracetamols
# 11         16.220216     17        Urine_metals
# 12          3.954558      4    Urine_metalloids
# 13        111.943716    127                 sum


