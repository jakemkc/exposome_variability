## Mar 28 2017
    ## Goal: % var explained by shared and non shared (ANOVA)
    ## Goal: paired t test between gender (taking out the shared/non shared effect)

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L



# ******** -----
# A. Load ready data ---------------------------------------------------------------------------

# Load data from corr_chem_m_f_mf_w...
## LIFE df ready with Y binary as numeric
load("results/LIFE_df_X_Y_imputed_factor_rdy.Rdata") 



# ******** -----
# B. Prep df for ANOVA & ttest ---------------------------------------------------------------------------

# >>> 1
# LIFE name check
colnames(life)

# remove male and female chems in LIFE
life2 <- life

intersect(colnames(life2), colnames(impute_m[3:length(colnames(impute_m))])) 
    # 133 chem, ignoring ID and imputation columes
intersect(colnames(life2), colnames(impute_f[3:length(colnames(impute_f))])) 
    # 133 chem, ignoring ID and imputation columes

life_183 <- life2[, setdiff(colnames(life2), colnames(impute_f[3:length(colnames(impute_f))]))]
life_50 <- life_183[, setdiff(colnames(life_183), colnames(impute_m[3:length(colnames(impute_m))]))]



# >>> 2
# Separate life_50 colname into female and male
dput(colnames(life_50))


life_50_f_name <- c("ID", "Age_f", "RaceEthnicity_f", "DFEDU", "DFHSINCO", "parity_f", "LFSMKNOW", "catbmi_f", "conditionalParity_f", "conceptiondelay", "INFERTILE", "LFEXERCS", 
                    "LFSMKCGD", "LFAL12DR", "LFETOHTY", "MFHIGHBP", "MFHGHCHL", "MFDIABTS", "amy_mean", "cort_mean", "SEVOLUME", "SPCOUNT", "TOTCNT", 
                    "PERMOT", "STCRIT", "WHONORM", "catsevolume", "catspcount", "cattotcnt", 
                    "catstcrit", "catwhonorm", "SCSADFI", "SCSAHDS", "catscsadfi", "catscsahds")

life_50_m_name <- c("ID", "Age_m", "RaceEthnicity_m", "DMEDU", "DMHSINCO", "parity_m", 
                    "LMSMKNOW", "catbmi_m","conditionalParity_m", "conceptiondelay", "INFERTILE", "LMEXERCS", "LMSKCIGD", "LMETOHBE", "LMETOHTY", "MMHIGHBP", "MMHIGHCH", 
                    "MMDIABET", "amy_mean", "cort_mean", "SEVOLUME", "SPCOUNT", "TOTCNT", 
                    "PERMOT", "STCRIT", "WHONORM", "catsevolume", "catspcount", "cattotcnt", "catstcrit", "catwhonorm", "SCSADFI", "SCSAHDS", "catscsadfi", "catscsahds")
    # Manually confirmed the life_50_name lists are correponding to each other



# >>> 3
# Separate life_50 df into female and male dfs, create x11 repeated
life_50_f <- life_50[, life_50_f_name]
life_50_m <- life_50[, life_50_m_name]

# Use female frame colnames
colnames(life_50_m) <- colnames(life_50_f)

# sort by ID
life_50_f <- life_50_f[order(life_50_f$ID),]
life_50_m <- life_50_m[order(life_50_m$ID),]

# x11 repeats
life_50_f_11x <- do.call("rbind", replicate(11, life_50_f, simplify = FALSE))
life_50_m_11x <- do.call("rbind", replicate(11, life_50_m, simplify = FALSE))



# >>> 4
# Setup NA and new vars
life_50_m_11x$amy_mean <- NA
life_50_m_11x$cort_mean <- NA

dput(which(colnames(life_50_f_11x) %in% c("SEVOLUME", "SPCOUNT", "TOTCNT", "PERMOT", "STCRIT", "WHONORM", "catsevolume", "catspcount", "cattotcnt", "catstcrit", "catwhonorm", "SCSADFI", "SCSAHDS", "catscsadfi", "catscsahds")))
    # 15 semen vars
life_50_f_11x[, 21:35] <- NA

# rename family var
colnames(life_50_f_11x)[1] <- "family"
colnames(life_50_m_11x)[1] <- "family"

# add gender var
life_50_f_11x$gender <- "F"
life_50_m_11x$gender <- "M"



# >>> 5
# create impute = 0 df from LIFE df
imp0_f <- life[, which(colnames(life) %in% colnames(impute_f))] 
imp0_f$X_Imputation_ <- 0
imp0_f <- imp0_f[, c(135, 1:134)]

imp0_m <- life[, which(colnames(life) %in% colnames(impute_m))] 
imp0_m$X_Imputation_ <- 0
imp0_m <- imp0_m[, c(135, 1:134)]
colnames(imp0_m) <- colnames(imp0_f)

# sort by ID
imp0_f <- imp0_f[order(imp0_f$ID),]
imp0_m <- imp0_m[order(imp0_m$ID),]


# >>> 6
# Confirm impute_f and impute_m chem names are corresponding & male uses female colnames
no_f <- gsub("_f$", "", colnames(impute_f), ignore.case = T)
no_m <- gsub("_m$", "", colnames(impute_m), ignore.case = T)
all(no_f == no_m)

impute_f_2 <- impute_f
impute_m_2 <- impute_m

colnames(impute_m_2) <- colnames(impute_f_2)

# sort by imputation then ID
impute_f_2 <- impute_f_2[order(impute_f_2$X_Imputation_, impute_f_2$ID),]
impute_m_2 <- impute_m_2[order(impute_m_2$X_Imputation_, impute_m_2$ID),]


# >>> 7
# rbind, cbind
impute_f_2 <- rbind(imp0_f, impute_f_2)
impute_m_2 <- rbind(imp0_m, impute_m_2)

tmp1 <- cbind(impute_f_2, life_50_f_11x)
tmp2 <- cbind(impute_m_2, life_50_m_11x)

lifemerg <- rbind(tmp1, tmp2)

# sort by imputation then ID
lifemerg <- lifemerg[order(lifemerg$X_Imputation_, lifemerg$ID),]
row.names(lifemerg) <- c(1:nrow(lifemerg))


# >>> 8
# prep factor and numeric vars
    # str(lifemerg[,133:171])

lifemerg$family <- as.factor(lifemerg$family)
    # contrasts(lifemerg$family)
    # table(lifemerg$fammily)

lifemerg$gender <- as.factor(lifemerg$gender)
    # lm(Age_f ~ gender, data = lifemerg) # average age = table 2 in paper, level as string or number doesn't matter but order of level

# save(file = "testmixed.Rdata", lifemerg)



# ******** -----
# C. Get lipid gender residuals ---------------------------------------------------------------------------

# >>> 1

## Female
lm_resid_fl <- function(indvar="lipids_f", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log(%s+1))) ~ scale(%s) + gender + scale(Age_f)", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}

lm_resid_fc <- function(indvar="CREAMOUNT_F", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log(%s+1))) ~ scale(%s) + gender + scale(Age_f)", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}

# log10 maybe better for more powerful transform in env data

# >>> 2
# chemical looping lists

# Create a full chem category and chemical list for both females and males

cl = list(
    #PCBs
    PCBsF = c("pcb028amt_F", "pcb044amt_F", "pcb049amt_F", 
              "pcb052amt_F", "pcb066amt_F", "pcb074amt_F", "pcb087amt_F", "pcb099amt_F", 
              "pcb101amt_F", "pcb105amt_F", "pcb110amt_F", "pcb114amt_F", "pcb118amt_F", 
              "pcb128amt_F", "pcb138amt_F", "pcb146amt_F", "pcb149amt_F", "pcb151amt_F", 
              "pcb153amt_F", "pcb156amt_F", "pcb157amt_F", "pcb167amt_F", "pcb170amt_F", 
              "pcb172amt_F", "pcb177amt_F", "pcb178amt_F", "pcb180amt_F", "pcb183amt_F", 
              "pcb187amt_F", "pcb189amt_F", "pcb194amt_F", "pcb195amt_F", "pcb196amt_F", 
              "pcb201amt_F", "pcb206amt_F", "pcb209amt_F"),
    #PFCs
    PFCsF = c("pfcepahamt_F", 
              "pfcmpahamt_F", "pfcpfdeamt_F", "pfcpfnaamt_F", "pfcpfsaamt_F", 
              "pfcpfosamt_F", "pfcpfoaamt_F"),
    #cotinine
    cotF = c("cotamt_F"),
    #blood_metals
    BmetF = c("metbcdamt_F", "metbpbamt_F", 
              "metthgamt_F"),
    #POPs
    POPsF = c("pophcbamt_F", "popmiramt_F", "popbb1amt_F", "popbhcamt_F", 
              "popghcamt_F", "popodtamt_F", "popoxyamt_F", "poppdeamt_F", "poppdtamt_F", 
              "poptnaamt_F"),
    #PBDEs
    PBDEsF = c("pbdebr1amt_F", "pbdebr2amt_F", "pbdebr3amt_F", 
               "pbdebr4amt_F", "pbdebr5amt_F", "pbdebr6amt_F", "pbdebr7amt_F", 
               "pbdebr8amt_F", "pbdebr9amt_F", "pbdebr66amt_F"),
    #serum_lipids
    serlipF = c("fcholamt_F", "cholamt_F", "trigamt_F", "phosamt_F"),
    #phytoestrogens
    phytoestF = c("DAZAMOUNT_F", "DMAAMOUNT_F", 
                  "EQUAMOUNT_F", "ETDAMOUNT_F", "ETLAMOUNT_F", "GNSAMOUNT_F"),
    #creatinine
    creatF = c("CREAMOUNT_F"),
    #phthalates
    phthaF = c("mMethylPhthalate_F", "mEthylPhthalate_F", "mCarboxyPropylPhthalate_F", 
               "mButylPhthalate_F", "mIsobutylPhthalate_F", "mCarboxyEthylPentylPhthalate_F", 
               "mCarboxyMethylHexylPhthalate_F", "mEthylHydroxyHexylPhthalate_F", 
               "mEthylOxoHexylPhthalate_F", "mCycloHexylPhthalate_F", "mBenzylPhthalate_F", 
               "mEthylHexylPhthalate_F", "mOctylPhthalate_F", "mIsononylPhthalate_F"),
    #BPA_benzophenone
    BPA_bzF = c("BPA_F", "HydroxyMethoxyBenzoPhenone_F", "HydroxyBenzoPhenone_F", 
                "DiHydroxyBenzoPhenone_F", "DiHydroxyMethoxyBenzoPhenone_F", 
                "TetraHydroxyBenzoPhenone_F"),
    #metalloids
    metalloF = c("Selenium_F", "Arsenic_F", "Antimony_F", "Tellurium_F"),
    #metals
    metalsF = c("Manganese_F", "Chromium_F", "Beryllium_F", "Cobalt_F", "Molybdenum_F", "Cadmium_Corrected_F", 
                "Tin_F", "Caesium_F", "Barium_F", "Nickel_F", "Copper_F", "Zinc_F", "Tungsten_F", "Platinum_F", 
                "Thallium_F", "Lead_F", "Uranium_F"),
    #paraben
    parabF = c("MeP_f", "EtP_f", "PrP_f", 
               "BuP_f", "BzP_f", "HeP_f", "X_4_HB_f", "X_3_4_DHB_f", "OH_MeP_f", 
               "OH_EtP_f", "TCS_f", "TCC_f", "PAP_f", "APAP_f")
)



# >>> 3
# Operators for lipids and creatinine
    # males using females labels in long format
list_f_oper <- list(
    F_creat = c("phytoestF", "phthaF", "BPA_bzF", "metalloF", "metalsF", "parabF"),
    F_lipid = c("PCBsF", "POPsF", "PBDEsF"),
    F_null = c("PFCsF", "cotF", "BmetF", "serlipF", "creatF"))



# >>> 4
# Getting residuals

# Take out impute=0 parts, loop x10 frames
lifemerg_tmp <- subset(lifemerg, X_Imputation_ > 0)

# loop
resid_retFrame <- data.frame(row.names = 1:nrow(lifemerg_tmp)) # for cbinding residual columns

for(i in 1:length(unlist(cl))){
    lifemerg_tmp$lipids_f <- (1.494*lifemerg_tmp$cholamt_F) + lifemerg_tmp$trigamt_F + lifemerg_tmp$phosamt_F  # new var (total lipids Philips formula)
    tmpfvar <- unlist(cl)[i]
    tmpfoper <- names(cl[sapply(cl, "%in%", x = tmpfvar)])
    finalfoper <- names(list_f_oper[sapply(list_f_oper, "%in%", x = tmpfoper)])
    
    if (finalfoper == "F_null"){
        matchindex <- which(colnames(lifemerg_tmp) %in% tmpfvar)
        frm <- scale(log(lifemerg_tmp[, matchindex]+1))
    } else if (finalfoper == "F_creat"){
        mod <- lm_resid_fc("CREAMOUNT_F", lifemerg_tmp, tmpfvar)
        frm <- as.data.frame(summary(mod)$residuals)
    } else {
        mod <- lm_resid_fl("lipids_f", lifemerg_tmp, tmpfvar)
        frm <- as.data.frame(summary(mod)$residuals)
    }
    colnames(frm) <- tmpfvar
    resid_retFrame <- cbind(resid_retFrame, frm)
}



# >>> 5
# rebuilding lifemerg with residual

dput(which(colnames(lifemerg_tmp) %in% colnames(resid_retFrame)))
lifemerg_tmp[, 3:135] <- resid_retFrame
lifemerg_tmp[, 172] <- NULL  # lipid var

index0 <- which(lifemerg$X_Imputation_ == 0, arr.ind = TRUE)
lifemerg_anova <- rbind(lifemerg[index0, ], lifemerg_tmp)

# sort by imputation then ID
lifemerg_anova <- lifemerg_anova[order(lifemerg_anova$X_Imputation_, lifemerg_anova$ID),]
row.names(lifemerg_anova) <- NULL



# ******** -----
# C2. MICE & LM adjusted R2 ---------------------------------------------------------------------------

library(mice)
library(miceadds)

# imp column as factor for MICE
lifemerg_anova$X_Imputation_ <- as.factor(lifemerg_anova$X_Imputation_)  # previous has problem if not a factor var
# str(lifemerg_anova)  # OK

# read in master df for MICE
lifemice_anova <- as.mids(data = lifemerg_anova, .imp=1, .id = NULL)  # if id not as NULL, will have import error



# customized function of pool.r.square
pool_r_squared_jake <- function (object, adjusted = FALSE) # don't provide adjusted = TRUE
{
    call <- match.call()
    if (!is.mira(object)) 
        stop("The object must have class 'mira'")
    if ((m <- length(object$analyses)) < 2) 
        stop("At least two imputations are needed for pooling.\n")
    if (class((object$analyses[[1]]))[1] != "lm") 
        stop("r^2 can only be calculated for results of the 'lm' modelling function")
    analyses <- object$analyses
    m <- length(analyses)
    r2 <- matrix(NA, nrow = m, ncol = 3, dimnames = list(1:m, 
                                                         c("R^2", "Fisher trans F^2", "se()")))
    for (i in 1:m) {
        fit <- analyses[[i]]
        if (adjusted == FALSE) 
            r2[i, 1] <- sqrt(summary(fit)$r.squared)
        else r2[i, 1] <- sqrt(summary(fit)$adj.r.squared)
        r2[i, 2] <- 0.5 * log((r2[i, 1] + 1)/(1 - r2[i, 1]))
        r2[i, 3] <- 1/(length(summary(fit)$residuals) - 3)
    }
    r2 <- r2[complete.cases(r2), ]
    fit <- pool.scalar(r2[, 2], r2[, 3])
    table <- array(((exp(2 * fit$qbar) - 1)/(1 + exp(2 * fit$qbar)))^2, 
                   dim = c(1, 4))
    if (adjusted == FALSE) 
        dimnames(table) <- list("R^2", c("est", "lo 95", "hi 95", 
                                         "fmi"))
    else dimnames(table) <- list("adj R^2", c("est", "lo 95", 
                                              "hi 95", "fmi"))
    table[, 2] <- ((exp(2 * (fit$qbar - 1.96 * sqrt(fit$t))) - 
                        1)/(1 + exp(2 * (fit$qbar - 1.96 * sqrt(fit$t)))))^2
    table[, 3] <- ((exp(2 * (fit$qbar + 1.96 * sqrt(fit$t))) - 
                        1)/(1 + exp(2 * (fit$qbar + 1.96 * sqrt(fit$t)))))^2
    table[, 4] <- fit$f
    return(table)
}


# model
lm_shared_adjr <- function(dat, depvar="pcb028amt_F") {
    setform <- sprintf("%s ~ family", depvar)
    mod2 <- with(data = dat, lm(formula = as.formula(setform)))
}

# Loop
retFrame_adjr <- data.frame()

# list (take out error one)
cl_adj <- cl
cl_adj$metalsF <- cl_adj$metalsF[1:16] # remove uranium. only has 1 valid case and an error to stop the loop
cl_adj$parabF <- cl_adj$parabF[-3] # PrP_f just all NAN


for(i in 1:length(unlist(cl_adj))) {
    DV <- unlist(cl_adj)[i]
    frm <- lm_shared_adjr(lifemice_anova, DV)
    r <- pool_r_squared_jake(frm)
    adjr <- pool_r_squared_jake(frm, adjust = TRUE)
    est <- as.data.frame(r[,1])
    names(est) <- "r2"
    est$adjr <- adjr[,1]
    est$DV <- DV
    print(DV)
    retFrame_adjr <- rbind(retFrame_adjr, est)
}


# Write
write.csv(retFrame_adjr, file = "results/ANOVA_share_non_shared_coef_age_adj_adjr_removed.csv")

# Read
# retFrame_adjr = read.csv("results/ANOVA_share_non_shared_coef_age_adj_adjr_removed.csv", stringsAsFactors = FALSE)  # read csv file 
# retFrame_adjr = read.csv("results/ANOVA_share_non_shared_coef_age_adj_adjr_Nan_to_1eminus10", stringsAsFactors = FALSE)  # read csv file 


# ******** -----
# C2. LM adj r & Boxplot ---------------------------------------------------------------------------

# drop covariates
drop_anova <- which(retFrame_adjr$DV %in% c("CREAMOUNT_F", "fcholamt_F", "cholamt_F", "trigamt_F", "phosamt_F"))
retFrame_adjr <- retFrame_adjr[-drop_anova, ]


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
no_f <- gsub("_f$", "", retFrame_adjr$DV, ignore.case = T)

# match and get the index order by female_no_f
index_f <- vector(mode = "numeric")
for(i in 1:length(no_f)) {
    mtchF <- which(papername_f %in% no_f[i])
    print(mtchF) # error check
    index_f <- append(index_f, mtchF)
}

# create new DV2
retFrame_adjr$DV2 <- names(papername_f[index_f])

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

# Match to categories
retFrame_adjr$categ <- "n____" # dummary for partial rename in the loop

for(i in 1:length(cl_f)) {
    mtchC <- which(retFrame_adjr$DV2 %in% cl_f[[i]])
    retFrame_adjr$categ[mtchC] <- names(cl_f)[i]
}


# ggplot
library(ggplot2)


# make same ranking level as female plot
female_lv_2 <- c("Polybrominated_cpds", "PCBs", "OCPs", "Blood_metals", "Phytoestrogens", 
                 "Urine_metals", "Phenols", "PFASs", "Phthalates", "Paracetamols", 
                 "Anti_microbial_cpds", "Urine_metalloids", "Cotinine")

# quartz(height=6, width=8)
p <- ggplot(data = retFrame_adjr, aes(x = factor(categ, levels = female_lv_2), y = adjr*100)) # -absR as revered order; absR as normal order
p <- p + geom_boxplot()
p <- p + ggtitle("% variance of biomarkers explained by the shared environment, adjusted R^2") +
    xlab("Chemical classes") + 
    ylab("% variance explained")
p <- p + ylim(0, 100)
# p <- p + geom_hline(yintercept = 0.08763547)
pf <- p + theme(axis.text.x  = element_text(angle=45, vjust=0.7, hjust= 0.7, size=10)) # move to the left a bit (using vjust since the labels are rotated)
pf


ggsave("results/boxplot_r2_varexplained_within_category_ggplot_adj_r2.png", scale=1, dpi=400)

# IQRs
do.call("rbind", tapply(retFrame_adjr$adjr, retFrame_adjr$categ, quantile))

#                     0%            25%        50%        75%      100%
# Anti_microbial_cpds 0.0138839046 0.04804445 0.14542115 0.20519344 0.3068417
# Blood_metals        0.1781315477 0.29343137 0.40873119 0.45638482 0.5040384
# Cotinine            0.2069823656 0.20698237 0.20698237 0.20698237 0.2069824
# OCPs                0.0652476456 0.09435673 0.11462487 0.18614548 0.6920669
# Paracetamols        0.0520965982 0.07845530 0.10481401 0.13117272 0.1575314
# PCBs                0.0250984199 0.05748563 0.08072945 0.11381976 0.7636008
# PFASs               0.0822843903 0.15153277 0.42524645 0.46799526 0.6702116
# Phenols             0.0383111863 0.05013000 0.06546339 0.12103391 0.1527078
# Phthalates          0.0008392527 0.02556641 0.04398021 0.08881305 0.2951607
# Phytoestrogens       0.0026057983 0.05829084 0.14917943 0.24099166 0.2989792
# Polybrominated_cpds 0.0378591859 0.05545279 0.11908613 0.17588867 0.5128096
# Urine_metalloids    0.0365316315 0.05046975 0.08037397 0.10649169 0.1090703
# Urine_metals        0.0090155687 0.05776970 0.10004113 0.19048747 0.3285290

# average of the median: 0.157282596 (no outlier effect here


# ******** -----
# D. Get lipid residuals ---------------------------------------------------------------------------

# >>> 1
# Set residual formula for males and females
# males using females labels in long format

## Female
lm_resid_fl_ttest <- function(indvar="lipids_f", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log(%s+1))) ~ scale(%s) + scale(Age_f)", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}

lm_resid_fc_ttest <- function(indvar="CREAMOUNT_F", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log(%s+1))) ~ scale(%s) + scale(Age_f)", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}



# >>> 2
# chemical looping lists

# Create a full chem category and chemical list for both females and males

cl = list(
    #PCBs
    PCBsF = c("pcb028amt_F", "pcb044amt_F", "pcb049amt_F", 
              "pcb052amt_F", "pcb066amt_F", "pcb074amt_F", "pcb087amt_F", "pcb099amt_F", 
              "pcb101amt_F", "pcb105amt_F", "pcb110amt_F", "pcb114amt_F", "pcb118amt_F", 
              "pcb128amt_F", "pcb138amt_F", "pcb146amt_F", "pcb149amt_F", "pcb151amt_F", 
              "pcb153amt_F", "pcb156amt_F", "pcb157amt_F", "pcb167amt_F", "pcb170amt_F", 
              "pcb172amt_F", "pcb177amt_F", "pcb178amt_F", "pcb180amt_F", "pcb183amt_F", 
              "pcb187amt_F", "pcb189amt_F", "pcb194amt_F", "pcb195amt_F", "pcb196amt_F", 
              "pcb201amt_F", "pcb206amt_F", "pcb209amt_F"),
    #PFCs
    PFCsF = c("pfcepahamt_F", 
              "pfcmpahamt_F", "pfcpfdeamt_F", "pfcpfnaamt_F", "pfcpfsaamt_F", 
              "pfcpfosamt_F", "pfcpfoaamt_F"),
    #cotinine
    cotF = c("cotamt_F"),
    #blood_metals
    BmetF = c("metbcdamt_F", "metbpbamt_F", 
              "metthgamt_F"),
    #POPs
    POPsF = c("pophcbamt_F", "popmiramt_F", "popbb1amt_F", "popbhcamt_F", 
              "popghcamt_F", "popodtamt_F", "popoxyamt_F", "poppdeamt_F", "poppdtamt_F", 
              "poptnaamt_F"),
    #PBDEs
    PBDEsF = c("pbdebr1amt_F", "pbdebr2amt_F", "pbdebr3amt_F", 
               "pbdebr4amt_F", "pbdebr5amt_F", "pbdebr6amt_F", "pbdebr7amt_F", 
               "pbdebr8amt_F", "pbdebr9amt_F", "pbdebr66amt_F"),
    #serum_lipids
    serlipF = c("fcholamt_F", "cholamt_F", "trigamt_F", "phosamt_F"),
    #phytoestrogens
    phytoestF = c("DAZAMOUNT_F", "DMAAMOUNT_F", 
                  "EQUAMOUNT_F", "ETDAMOUNT_F", "ETLAMOUNT_F", "GNSAMOUNT_F"),
    #creatinine
    creatF = c("CREAMOUNT_F"),
    #phthalates
    phthaF = c("mMethylPhthalate_F", "mEthylPhthalate_F", "mCarboxyPropylPhthalate_F", 
               "mButylPhthalate_F", "mIsobutylPhthalate_F", "mCarboxyEthylPentylPhthalate_F", 
               "mCarboxyMethylHexylPhthalate_F", "mEthylHydroxyHexylPhthalate_F", 
               "mEthylOxoHexylPhthalate_F", "mCycloHexylPhthalate_F", "mBenzylPhthalate_F", 
               "mEthylHexylPhthalate_F", "mOctylPhthalate_F", "mIsononylPhthalate_F"),
    #BPA_benzophenone
    BPA_bzF = c("BPA_F", "HydroxyMethoxyBenzoPhenone_F", "HydroxyBenzoPhenone_F", 
                "DiHydroxyBenzoPhenone_F", "DiHydroxyMethoxyBenzoPhenone_F", 
                "TetraHydroxyBenzoPhenone_F"),
    #metalloids
    metalloF = c("Selenium_F", "Arsenic_F", "Antimony_F", "Tellurium_F"),
    #metals
    metalsF = c("Manganese_F", "Chromium_F", "Beryllium_F", "Cobalt_F", "Molybdenum_F", "Cadmium_Corrected_F", 
                "Tin_F", "Caesium_F", "Barium_F", "Nickel_F", "Copper_F", "Zinc_F", "Tungsten_F", "Platinum_F", 
                "Thallium_F", "Lead_F", "Uranium_F"),
    #paraben
    parabF = c("MeP_f", "EtP_f", "PrP_f", 
               "BuP_f", "BzP_f", "HeP_f", "X_4_HB_f", "X_3_4_DHB_f", "OH_MeP_f", 
               "OH_EtP_f", "TCS_f", "TCC_f", "PAP_f", "APAP_f")
)



# >>> 3
# Operators for lipids and creatinine
# males using females labels in long format
list_f_oper <- list(
    F_creat = c("phytoestF", "phthaF", "BPA_bzF", "metalloF", "metalsF", "parabF"),
    F_lipid = c("PCBsF", "POPsF", "PBDEsF"),
    F_null = c("PFCsF", "cotF", "BmetF", "serlipF", "creatF"))



# >>> 4
# Getting residuals

# Take out impute=0 parts, loop x10 frames
lifemerg_tmp_ttest <- subset(lifemerg, X_Imputation_ > 0)

# loop
resid_retFrame_ttest <- data.frame(row.names = 1:nrow(lifemerg_tmp_ttest)) # for cbinding residual columns

for(i in 1:length(unlist(cl))){
    lifemerg_tmp_ttest$lipids_f <- (1.494*lifemerg_tmp_ttest$cholamt_F) + lifemerg_tmp_ttest$trigamt_F + lifemerg_tmp_ttest$phosamt_F  # new var (total lipids Philips formula)
    tmpfvar <- unlist(cl)[i]
    tmpfoper <- names(cl[sapply(cl, "%in%", x = tmpfvar)])
    finalfoper <- names(list_f_oper[sapply(list_f_oper, "%in%", x = tmpfoper)])
    
    if (finalfoper == "F_null"){
        matchindex <- which(colnames(lifemerg_tmp_ttest) %in% tmpfvar)
        frm <- scale(log(lifemerg_tmp_ttest[, matchindex]+1))
    } else if (finalfoper == "F_creat"){
        mod <- lm_resid_fc_ttest("CREAMOUNT_F", lifemerg_tmp_ttest, tmpfvar)
        frm <- as.data.frame(summary(mod)$residuals)
    } else {
        mod <- lm_resid_fl_ttest("lipids_f", lifemerg_tmp_ttest, tmpfvar)
        frm <- as.data.frame(summary(mod)$residuals)
    }
    colnames(frm) <- tmpfvar
    resid_retFrame_ttest <- cbind(resid_retFrame_ttest, frm)
}



# >>> 5
# rebuilding lifemerg with residual

dput(which(colnames(lifemerg_tmp_ttest) %in% colnames(resid_retFrame_ttest)))
lifemerg_tmp_ttest[, 3:135] <- resid_retFrame
lifemerg_tmp_ttest[, 172] <- NULL  # lipid var

index0 <- which(lifemerg$X_Imputation_ == 0, arr.ind = TRUE)
lifemerg_ttest <- rbind(lifemerg[index0, ], lifemerg_tmp_ttest)

# sort by imputation then ID
lifemerg_ttest <- lifemerg_ttest[order(lifemerg_ttest$X_Imputation_, lifemerg_ttest$ID),]
row.names(lifemerg_ttest) <- NULL


# >>> 6
# prep amelia df
lifemerg_tmp_ttest <- lifemerg_tmp_ttest[order(lifemerg_tmp_ttest$X_Imputation_, lifemerg_tmp_ttest$ID),]

lamelia_f <- subset(lifemerg_tmp_ttest, gender == "F")

lamelia_m <- subset(lifemerg_tmp_ttest, gender == "M")

# all(colnames(lamelia_f) == colnames(lamelia_m))  # pairwise name check
# [1] TRUE

la_m_names <- colnames(lamelia_m)

la_m_names_tmp <- sapply(la_m_names, function(x) paste0(x,"_m"))  
colnames(lamelia_m) <- la_m_names_tmp

lamelia_mf <- cbind(lamelia_f, lamelia_m)

lamelia_tmp <-list()

for(i in 1:10) {
    lamelia_tmp[[i]] <- subset(lamelia_mf, X_Imputation_ == i)
}

lamelia_ttest <- list("imputations" = lamelia_tmp)

class(lamelia_ttest) <- "amelia"

# ******** -----
# D. Amelia & Paired t ---------------------------------------------------------------------------
library(MKmisc)

# test
# trick: it's calling the imputation list directly. No need to build amelia object indeed
# mi.t.test(lamelia_ttest$imputations, x = "pcb028amt_F", y = "pcb028amt_F_m", paired = TRUE)

retFrame_ttest <- data.frame()

# for loop
for (i in 1:length(unlist(cl))){
    FV <- unlist(cl)[i]
    MV <- sprintf("%s_m", FV)
    frm <- mi.t.test(lamelia_ttest$imputations, x = FV, y = MV, paired = TRUE)
    est <- data.frame(frm$p.value)
    names(est) <- "pvalue"
    est$t <- frm$statistic
    est$df <- frm$parameter
    est$meandiff <- frm$estimate[1]
    est$var <- FV
    retFrame_ttest <- rbind(retFrame_ttest, est)
}

# Write
write.csv(retFrame_ttest, file = "results/ttest_gender_shared_out_age_adj.csv")

# Read
# retFrame_ttest = read.csv("results/ttest_gender_shared_out_age_adj.csv", stringsAsFactors=FALSE)  # read csv file 

# drop covariates
drop_ttest <- which(retFrame_ttest$var %in% c("CREAMOUNT_F", "fcholamt_F", "cholamt_F", "trigamt_F", "phosamt_F"))
retFrame_ttest <- retFrame_ttest[-drop_ttest, ]

# fdr  # Use BY (Benjamini & Yekutieli) to control the FDR for positive dependence if needed
retFrame_ttest$fdr <- p.adjust(retFrame_ttest$pvalue, method = "fdr")

# Get the chemical names by p value
table(retFrame_ttest$pvalue <= 0.05) # 9/128 ~7%
retFrame_ttest$var[retFrame_ttest$pvalue <= 0.05]  # what they are

# [1] "pfcpfdeamt_F" "pfcpfnaamt_F" "pfcpfsaamt_F" "pfcpfosamt_F" "pfcpfoaamt_F" "cotamt_F"     "metbcdamt_F" 
# [8] "metbpbamt_F"  "metthgamt_F" 


## Create proper names
# match and get the index order by female_no_f
index_f_ttest <- vector(mode = "numeric")
for(i in 1:length(no_f)) {
    mtchF <- which(papername_f %in% no_f[i])
    print(mtchF) # error check
    index_f_ttest <- append(index_f_ttest, mtchF)
}

# create new var2
retFrame_ttest$var2 <- names(papername_f[index_f_ttest])

# Get the proper chemical names by fdr
table(retFrame_ttest$fdr <= 0.05) # 8/128 ~ 6%
retFrame_ttest$var[retFrame_ttest$fdr <= 0.05]

# [1] "pfcpfdeamt_F" "pfcpfnaamt_F" "pfcpfsaamt_F" "pfcpfosamt_F" "pfcpfoaamt_F" "cotamt_F"     "metbpbamt_F" 
# [8] "metthgamt_F" 

# plot
library(ggplot2)

retFrame_ttest$group <- "t_test"

# quartz(height=6, width=8)
# p <- ggplot(retFrame_ttest, aes(fdr)) + stat_ecdf()
p <- ggplot(retFrame_ttest, aes(fdr, colour = group)) + stat_ecdf()
p <- p + ggtitle("Cumulative Density plot of FDR adjusted p values") +
    xlab("p") + 
    ylab("cumulative density")
p

ggsave("results/Cum_freq_ggplot_fdr_t_test_gender_shared.png", scale=1, dpi=400)


