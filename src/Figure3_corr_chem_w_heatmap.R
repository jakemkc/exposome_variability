## Dec 14 2016
## Goal: correlation of chemicals within male, female and couple


rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings

# install.packages("haven")
library(haven)

life <- read_sas("data/LIFE_exposome.sas7bdat", catalog_file = "data/formats.sas7bcat")
# ?read_sas, read both data file and catalog file. The latter is for recording value labels
impute_m <- read_sas("data/imputedchemicals_male.sas7bdat")
impute_f <- read_sas("data/imputedchemicals_female.sas7bdat")

tlife <- as.data.frame(t(head(life))) # for View()


# ******** -----
# A. Cleanup ---------------------------------------------------------------------------

######################################################### #
## \\ A1. Name validity ----

# Question: Are the variable names valid?
# Answer: No, some are invalid and replace
######################################################### #


# Equal after variable check?
colnames(life) == make.names(colnames(life))
# Which ones?
colnames(life[!colnames(life) == make.names(colnames(life))])
# [1] "_4_HB_f"    "_3_4_DHB_f" "_4_HB_m"    "_3_4_DHB_m"

colnames(life) <- make.names(colnames(life))
# as "X_4_HB_m"

# update
tlife <- as.data.frame(t(head(life))) # for View()

# Equal after variable check?
colnames(impute_f) == make.names(colnames(impute_f))
# Which ones?
colnames(impute_f[!colnames(impute_f) == make.names(colnames(impute_f))])
# [1] "_Imputation_" "_4_HB_f"      "_3_4_DHB_f"  

colnames(impute_f) <- make.names(colnames(impute_f))

# Equal after variable check?
colnames(impute_m) == make.names(colnames(impute_m))
# Which ones?
colnames(impute_m[!colnames(impute_m) == make.names(colnames(impute_m))])
# [1] "_Imputation_" "_4_HB_m"      "_3_4_DHB_m"  
# update
colnames(impute_m) <- make.names(colnames(impute_m))


######################################################### #
## \\ A2. "" (blk) & " " (spc) as NA ----
######################################################### #

life <- as.data.frame(apply(life, 2, function(x) gsub("^$|^ $", NA, x)), stringsAsFactors = FALSE)

# attribute all gone and var become chr
str(life); class(life); mode(life)

# to numeric
life <- sapply(life, as.numeric)
str(life); class(life); mode(life)

# from matrix numeric to data.frame numeric 
life <- as.data.frame(life, stringsAsFactors = FALSE)

summary(life)



# ******** -----
# B. Setup vars ---------------------------------------------------------------------------


######################################################### #
## \\ B1. X adj. covar imputation ----

## Y:       Use complete.case
## X chem:  MI impute
## X cov:   Imputing here (see desc and cleanup.R)
######################################################### #


######################################################### #
## Female
######################################################### #

xcateg <- c("catbmi_f", "RaceEthnicity_f", "DFEDU", "DFHSINCO", 
            "conditionalParity_f", "parity_f", "LFSMKNOW", "LFEXERCS", 
            "LFAL12DR", "MFHIGHBP", "MFHGHCHL", "MFDIABTS")

# apply table function to each of the extracted columns
xcategout <- apply(life[, xcateg], 2, table, useNA = "always")

# Prep condensed table: Combine list to df
indxy <- sapply(xcategout, length) # length of each sublist

xcaetegtable <- as.data.frame(do.call(rbind, lapply(xcategout, "length<-", max(indxy))))


# Impute with most freq. value

for(i in 1:length(xcateg)) {
    life[, xcateg[i]][is.na(life[, xcateg[i]])] <- tail(names(sort(table(life[, xcateg[i]]))), 1) 
    # life$catbmi_f[is.na(life$catbmi_f)] <- tail(names(sort(table(life$catbmi_f))), 1) 
}


######################################################### #
## Male
######################################################### #

xcateg <- c("catbmi_m", "RaceEthnicity_m", "DMEDU", "DMHSINCO", 
            "conditionalParity_m", "parity_m", "LMSMKNOW", "LMEXERCS", 
            "LMETOHBE", "MMDIABET", "MMHIGHBP", "MMHIGHCH")

# apply table function to each of the extracted columns
xcategout <- apply(life[, xcateg], 2, table, useNA = "always")

# Prep condensed table: Combine list to df
indxy <- sapply(xcategout, length) # length of each sublist

xcaetegtable <- as.data.frame(do.call(rbind, lapply(xcategout, "length<-", max(indxy))))


# Impute with most freq. value

for(i in 1:length(xcateg)) {
    life[, xcateg[i]][is.na(life[, xcateg[i]])] <- tail(names(sort(table(life[, xcateg[i]]))), 1) 
    # life$catbmi_f[is.na(life$catbmi_f)] <- tail(names(sort(table(life$catbmi_f))), 1) 
}


######################################################### #
## \\ B2. Var types and str prep  ----
######################################################### #


######################################################### #
## \\ B2.1. Y as numeric in life (male, female)
######################################################### #

# binary var, no matter Y or X, doesn't matter to be factor or numeric variables
# >2 lv cateogical as X, then it matters

ylist9 <- c("conceptiondelay", "INFERTILE", "catsevolume", "catspcount", "cattotcnt", "catstcrit", "catwhonorm", "catscsadfi", "catscsahds")

for(i in 1:length(ylist9)) {
    ind <- match(ylist9[i], colnames(life))
    life[,ind] <- as.numeric(life[,ind])
}

######################################################### #
## \\ B2.2. X categorical covariates as factor vars (male, female)
######################################################### #

## female categorical x12
ffactlist <- c("catbmi_f", "RaceEthnicity_f", "DFEDU", "DFHSINCO", "conditionalParity_f", "parity_f", "LFSMKNOW", "LFEXERCS", "LFAL12DR", "MFHIGHBP", "MFHGHCHL", "MFDIABTS")

for(i in 1:length(ffactlist)) {
    ind <- match(ffactlist[i], colnames(life))
    life[,ind] <- as.factor(life[,ind])
}


## male categorical x12
mfactlist <- c("catbmi_m", "RaceEthnicity_m", "DMEDU", "DMHSINCO", "conditionalParity_m", "parity_m", "LMSMKNOW", "LMEXERCS", "LMETOHBE", "MMDIABET", "MMHIGHBP", "MMHIGHCH")

for(i in 1:length(mfactlist)) {
    ind <- match(mfactlist[i], colnames(life))
    life[,ind] <- as.factor(life[,ind])
}

# str(life) # check 



# ******** -----
## B-1. Saving Rdata ----

outdirectory <- "results"
outfilename <- "LIFE_df_X_Y_imputed_factor_rdy.Rdata"
# outfilename <- sprintf("%s_reg_7.Rdata", depVariable)
save(file=file.path(outdirectory,outfilename), 
     life, impute_f, impute_m)


# Load data
# load("results/LIFE_df_X_Y_imputed_factor_rdy.Rdata") 



# ******** -----
# B. On imputed DF  -------------------------------------------


## 2. Sorted by imputation & ID
library(dplyr)
impute_f <- arrange(impute_f, X_Imputation_, ID) # sort by X_imputation then ID. use desc(mpg) for descending
impute_m <- arrange(impute_m, X_Imputation_, ID) # sort by X_imputation then ID



# 3. Model regresses against covariate(s)

## Female
lm_resid_fl <- function(indvar="lipids_f", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log10(%s+1))) ~ I(scale(%s))", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}

# Some cmd takes with "", some don't for formula
# > setform;    # [1] "pcb028amt_F ~ family"
# > as.formula(setform); # pcb028amt_F ~ family

lm_resid_fc <- function(indvar="CREAMOUNT_F", dat, depvar="pcb028amt_F") {
    setform <- sprintf("I(scale(log10(%s+1))) ~ I(scale(%s))", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}


## Male
lm_resid_ml <- function(indvar="lipids_m", dat, depvar="pcb028amt_M") {
    setform <- sprintf("I(scale(log10(%s+1))) ~ I(scale(%s))", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}
    
lm_resid_mc <- function(indvar="CREAMOUNT_M", dat, depvar="pcb028amt_M") {
    setform <- sprintf("I(scale(log10(%s+1))) ~ I(scale(%s))", depvar, indvar)
    mod <- lm(as.formula(setform), data = dat)
}


# 3.1. Create list for matching in the loop to choose model

## Create a full chem category and chemical list
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
               "OH_EtP_f", "TCS_f", "TCC_f", "PAP_f", "APAP_f"),
    
    ## for male
    #PCBs
    PCBsM = c("pcb028amt_M", "pcb044amt_M", "pcb049amt_M", 
              "pcb052amt_M", "pcb066amt_M", "pcb074amt_M", "pcb087amt_M", "pcb099amt_M", 
              "pcb101amt_M", "pcb105amt_M", "pcb110amt_M", "pcb114amt_M", "pcb118amt_M", 
              "pcb128amt_M", "pcb138amt_M", "pcb146amt_M", "pcb149amt_M", "pcb151amt_M", 
              "pcb153amt_M", "pcb156amt_M", "pcb157amt_M", "pcb167amt_M", "pcb170amt_M", 
              "pcb172amt_M", "pcb177amt_M", "pcb178amt_M", "pcb180amt_M", "pcb183amt_M", 
              "pcb187amt_M", "pcb189amt_M", "pcb194amt_M", "pcb195amt_M", "pcb196amt_M", 
              "pcb201amt_M", "pcb206amt_M", "pcb209amt_M"),
    #PFCs
    PFCsM = c("pfcepahamt_M", 
              "pfcmpahamt_M", "pfcpfdeamt_M", "pfcpfnaamt_M", "pfcpfsaamt_M", 
              "pfcpfosamt_M", "pfcpfoaamt_M"),
    #cotinine
    cotM = c("cotamt_M"),
    #blood_metals
    BmetM = c("metbcdamt_M", "metbpbamt_M", 
              "metthgamt_M"),
    #POPs
    POPsM = c("pophcbamt_M", "popmiramt_M", "popbb1amt_M", "popbhcamt_M", 
              "popghcamt_M", "popodtamt_M", "popoxyamt_M", "poppdeamt_M", "poppdtamt_M", 
              "poptnaamt_M"),
    #PBDEs
    PBDEsM = c("pbdebr1amt_M", "pbdebr2amt_M", "pbdebr3amt_M", 
               "pbdebr4amt_M", "pbdebr5amt_M", "pbdebr6amt_M", "pbdebr7amt_M", 
               "pbdebr8amt_M", "pbdebr9amt_M", "pbdebr66amt_M"),
    #serum_lipids
    serlipM = c("fcholamt_M", 
                "cholamt_M", "trigamt_M", "phosamt_M"),
    #phytoestrogens
    phytoestM = c("DAZAMOUNT_M", "DMAAMOUNT_M", 
                  "EQUAMOUNT_M", "ETDAMOUNT_M", "ETLAMOUNT_M", "GNSAMOUNT_M"),
    #creatinine
    creatM = c("CREAMOUNT_M"),
    #phthalates
    phthaM = c("mMethylPhthalate_M", "mEthylPhthalate_M", "mCarboxyPropylPhthalate_M", 
               "mButylPhthalate_M", "mIsobutylPhthalate_M", "mCarboxyEthylPentylPhthalate_M", 
               "mCarboxyMethylHexylPhthalate_M", "mEthylHydroxyHexylPhthalate_M", 
               "mEthylOxoHexylPhthalate_M", "mCycloHexylPhthalate_M", "mBenzylPhthalate_M", 
               "mEthylHexylPhthalate_M", "mOctylPhthalate_M", "mIsononylPhthalate_M"),
    #BPA_benzophenone
    BPA_bzM = c("BPA_M", "HydroxyMethoxyBenzoPhenone_M", "HydroxyBenzoPhenone_M", 
                "DiHydroxyBenzoPhenone_M", "DiHydroxyMethoxyBenzoPhenone_M", 
                "TetraHydroxyBenzoPhenone_M"),
    #metalloids
    metalloM = c("Selenium_M", "Arsenic_M", "Antimony_M", "Tellurium_M"),
    #metals
    metalsM = c("Manganese_M", 
                "Chromium_M", "Beryllium_M", "Cobalt_M", "Molybdenum_M", "Cadmium_Corrected_M", 
                "Tin_M", "Caesium_M", "Barium_M", "Nickel_M", "Copper_M", "Zinc_M", "Tungsten_M", "Platinum_M", 
                "Thallium_M", "Lead_M", "Uranium_M"),
    #paraben
    parabM = c("MeP_m", "EtP_m", "PrP_m", "BuP_m", "BzP_m", "HeP_m", "X_4_HB_m", "X_3_4_DHB_m", "OH_MeP_m", 
               "OH_EtP_m", "TCS_m", "TCC_m", "PAP_m", "APAP_m")
)


# Operators 
list_f_oper <- list(
F_creat = c("phytoestF", "phthaF", "BPA_bzF", "metalloF", "metalsF", "parabF"),
F_lipid = c("PCBsF", "POPsF", "PBDEsF"),
F_null = c("PFCsF", "cotF", "BmetF", "serlipF", "creatF"))

list_m_oper <- list(
M_creat = c("phytoestM", "phthaM", "BPA_bzM", "metalloM", "metalsM", "parabM"),
M_lipid = c("PCBsM", "POPsM", "PBDEsM"),
M_null = c("PFCsM", "cotM", "BmetM", "serlipM", "creatM"))

# match: names(cl[sapply(cl, '%in%', x="pbdebr66amt_F")])


# 4. on nested loop to create the impute_F like DF with residual

## Female
resid_retFrame <- data.frame()  # the frame to combine for 1st lv loop
retFrame <- data.frame(row.names = 1:nrow(life)) # cbind can't take empty frame as rbind, for 2nd lv loop

for(i in 1:10){
    tmpf <- subset(impute_f, X_Imputation_ == i) # get an impute df out
    tmpf$lipids_f <- (1.494*tmpf$cholamt_F)+tmpf$trigamt_F+tmpf$phosamt_F  # new var (total lipids Philips formula)
    tmpf <- tmpf[, c(136, 1:135)]  # Rearrange the new var
        
        for(j in 1:(length(colnames(tmpf))-3)){         #Looping getting residual for each col of an selected imputed df
            tmpfvar <- colnames(tmpf)[j+3]
            tmpfoper <- names(cl[sapply(cl, "%in%", x = tmpfvar)])
            finalfoper <- names(list_f_oper[sapply(list_f_oper, "%in%", x = tmpfoper)])
            
            if (finalfoper == "F_null"){
                frm <- scale(log10(tmpf[j+3]+1))
            } else if (finalfoper == "F_creat"){
                mod <- lm_resid_fc("CREAMOUNT_F", tmpf, tmpfvar)
                frm <- as.data.frame(summary(mod)$residuals)
            } else {
                mod <- lm_resid_fl("lipids_f", tmpf, tmpfvar)
                frm <- as.data.frame(summary(mod)$residuals)
            }
            
            colnames(frm) <- tmpfvar 
            retFrame <- cbind(retFrame, frm)  
        }
    
    retFrame$X_Imputation_ <- tmpf$X_Imputation_     # Add an imputation column
    resid_retFrame <- rbind(resid_retFrame, retFrame)   # rind the residual df, imputation = 1:10
    retFrame <- data.frame(row.names = 1:nrow(life))  # reset the retFrame for another round of imputed df loop
}

resid_impute_f <- resid_retFrame
resid_impute_f <- resid_impute_f[, c(134, 1:133)]
    
    
## male
resid_retFrame <- data.frame()  # the frame to combine for 1st lv loop
retFrame <- data.frame(row.names = 1:nrow(life)) # cbind can't take empty frame as rbind, for 2nd lv loop


for(i in 1:10){
    tmpm <- subset(impute_m, X_Imputation_ == i) # get an impute df out
    tmpm$lipids_m <- (1.494*tmpm$cholamt_M)+tmpm$trigamt_M+tmpm$phosamt_M  # new var (total lipids Philips formula)
    tmpm <- tmpm[, c(136, 1:135)]  # Rearrange the new var
    
        for(j in 1:(length(colnames(tmpm))-3)){         #Looping getting residual for each col of an selected imputed df
            tmpmvar <- colnames(tmpm)[j+3]
            tmpmoper <- names(cl[sapply(cl, "%in%", x = tmpmvar)])
            finalmoper <- names(list_m_oper[sapply(list_m_oper, "%in%", x = tmpmoper)])
            
            if (finalmoper == "M_null"){
                frm <- scale(log10(tmpm[j+3]+1))
            } else if (finalmoper == "M_creat"){
                mod <- lm_resid_mc("CREAMOUNT_M", tmpm, tmpmvar)
                frm <- as.data.frame(summary(mod)$residuals)
            } else {
                mod <- lm_resid_ml("lipids_m", tmpm, tmpmvar)
                frm <- as.data.frame(summary(mod)$residuals)
            }
            
            colnames(frm) <- tmpmvar
            retFrame <- cbind(retFrame, frm)  
        }
    
    retFrame$X_Imputation_ <- tmpm$X_Imputation_     # Add an imputation column
    resid_retFrame <- rbind(resid_retFrame, retFrame)   # rind the residual df, imputation = 1:10
    retFrame <- data.frame(row.names = 1:nrow(life))  # reset the retFrame for another round of imputed df loop
}

resid_impute_m <- resid_retFrame
resid_impute_m <- resid_impute_m[, c(134, 1:133)]





######################################################### #
# \\ B2. Corr matrices (133 chems) -------------------------------------------------------
######################################################### #



######################################################### #
## ++ B2.1 for within impute_f ----

library(psych)
corrlistf <- list()

for(i in 1:10) {
    tmpf <- subset(resid_impute_f, X_Imputation_ == i)[2:ncol(resid_impute_f)]
    cor_r <- as.matrix(corr.test(tmpf, method="spearman")$r)   # "cor" can take 2 matrices too
    corrlistf[[i]] <- cor_r
}

# Pooling estimate of r by Rubin method
bindlistf <- do.call(cbind, corrlistf) # corrlist elements must be in matrix format
bindarrayf <- array(bindlistf, dim=c(dim(corrlistf[[1]]), length(corrlistf)))  
corrubinf <- apply(bindarrayf, c(1, 2), mean, na.rm = TRUE)  # c(1, 2) indicates rows and columns; apply mean rowwise, columnwise

colnames(corrubinf) <- colnames(resid_impute_f[,2:ncol(resid_impute_f)])
rownames(corrubinf) <- colnames(resid_impute_f[,2:ncol(resid_impute_f)])


######################################################### #
##   ++ B2.2 for within impute_m ----

library(psych)
corrlistm <- list()

for(i in 1:10) {
    tmpm <- subset(resid_impute_m, X_Imputation_ == i)[2:ncol(resid_impute_m)]
    cor_r <- as.matrix(corr.test(tmpm, method="spearman")$r) 
    corrlistm[[i]] <- cor_r
}

# Pooling estimate of r by Rubin method
bindlistm <- do.call(cbind, corrlistm) # corrlist elements must be in matrix format
bindarraym <- array(bindlistm, dim=c(dim(corrlistm[[1]]), length(corrlistm)))  
corrubinm <- apply(bindarraym, c(1, 2), mean, na.rm = TRUE)  # c(1, 2) indicates rows and columns

colnames(corrubinm) <- colnames(resid_impute_m[,2:ncol(resid_impute_m)])
rownames(corrubinm) <- colnames(resid_impute_m[,2:ncol(resid_impute_m)])



######################################################### #
## ++ B2.3 for within imputed couple ----
library(psych)
corrlistmf <- list()

for(i in 1:10) {
    tmpm <- subset(resid_impute_m, X_Imputation_ == i)[2:ncol(resid_impute_m)]
    tmpf <- subset(resid_impute_f, X_Imputation_ == i)[2:ncol(resid_impute_f)]
    cor_r <- as.matrix(corr.test(tmpm, tmpf, method = "spearman")$r) 
    corrlistmf[[i]] <- cor_r
}

# Pooling estimate of r by Rubin method
bindlistmf <- do.call(cbind, corrlistmf) # corrlist elements must be in matrix format
bindarraymf <- array(bindlistmf, dim=c(dim(corrlistmf[[1]]), length(corrlistmf)))  
corrubinmf <- apply(bindarraymf, c(1, 2), mean, na.rm = TRUE)  # c(1, 2) indicates rows and columns


colnames(corrubinmf) <- colnames(resid_impute_f[,2:ncol(resid_impute_f)]) #check the cor_r to know f as colname
rownames(corrubinmf) <- colnames(resid_impute_m[,2:ncol(resid_impute_m)])







# ******** -----
## C. Set variables for figures ----

######################################################### #
## \\ C0-1. females ----
######################################################### #


## 1. Replace the variable names in corrubin to mine short_form reporting in the paper

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
no_f <- gsub("_f$", "", colnames(corrubinf), ignore.case = T)

# match and get the index order by female_no_f
index_f <- vector(mode = "numeric")
for(i in 1:length(no_f)) {
    mtchF <- which(papername_f %in% no_f[i])
    print(mtchF) # error check
    index_f <- append(index_f, mtchF)
}


# rename corrubin with my short-form
colnames(corrubinf) <- names(papername_f[index_f])
rownames(corrubinf) <- names(papername_f[index_f])


## Drop creatinine and lipids in row and column of the corrbuinf
dput(which(rownames(corrubinf) %in% c("CREAMOUNT_f", "fcholamt_f", "cholamt_f", "trigamt_f", "phosamt_f")))
corrubinf <- corrubinf[-c(68L, 69L, 70L, 71L, 78L), -c(68L, 69L, 70L, 71L, 78L)]

# B. resort the variable so they will be more consistent with the class order in the paper

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

# match and get the index
index_f_2 <- vector(mode = "numeric")
for(i in 1:length(unlist(cl_f))) {
    mtchF <- which(colnames(corrubinf) %in% unlist(cl_f)[i])
    print(mtchF) # error check
    index_f_2 <- append(index_f_2, mtchF)
}

# rearrange the order of the 128 corrubin matrix
corrubinf <- corrubinf[index_f_2, index_f_2] 
# diag(corrubinf)  #OK





######################################################### #
## \\ C0-2. males ----
######################################################### #


## 1. Replace the variable names in corrubin to mine short_form reporting in the paper

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
no_m <- gsub("_m$", "", colnames(corrubinm), ignore.case = T)

# match and get the index order by female_no_f
index_m <- vector(mode = "numeric")
for(i in 1:length(no_m)) {
    mtchM <- which(papername_m %in% no_m[i])
    print(mtchM) # error check
    index_m <- append(index_m, mtchM)
}

# rename corrubin with my short-form
colnames(corrubinm) <- names(papername_m[index_m])
rownames(corrubinm) <- names(papername_m[index_m])


## Drop creatinine and lipids in row and column of the corrbuin
dput(which(rownames(corrubinm) %in% c("CREAMOUNT_m", "fcholamt_m", "cholamt_m", "trigamt_m", "phosamt_m")))
corrubinm <- corrubinm[-c(68L, 69L, 70L, 71L, 78L), -c(68L, 69L, 70L, 71L, 78L)]

# B. resort the variable so they will be more consistent with the class order in the paper

## List in right order (13 classes, 128 chemicals wo lipids and creatinine)
cl_m <- cl_f

# replace the _f in the cl_m to _m
for(i in 1:length(cl_m)){
    cl_m[[i]] <- gsub("_f$", "_m", cl_m[[i]], ignore.case = T)
    # print(cl_m[[i]])
}


# match and get the index
index_m_2 <- vector(mode = "numeric")
for(i in 1:length(unlist(cl_m))) {
    mtchM <- which(colnames(corrubinm) %in% unlist(cl_m)[i])
    print(mtchM) # error check
    index_m_2 <- append(index_m_2, mtchM)
}

# rearrange the order of the 128 corrubin matrix
corrubinm <- corrubinm[index_m_2, index_m_2] 
# diag(corrubinm)  #OK





######################################################### #
## \\ C0-3. couples ----
######################################################### #


## 1. Replace the variable names in corrubin to mine short_form reporting in the paper

# rename corrubin with my short-form
colnames(corrubinmf) <- names(papername_f[index_f])
rownames(corrubinmf) <- names(papername_m[index_m])


## Drop creatinine and lipids in row and column of the corrbuin
dput(which(rownames(corrubinmf) %in% c("CREAMOUNT_m", "fcholamt_m", "cholamt_m", "trigamt_m", "phosamt_m")))
dput(which(colnames(corrubinmf) %in% c("CREAMOUNT_f", "fcholamt_f", "cholamt_f", "trigamt_f", "phosamt_f")))
corrubinmf <- corrubinmf[-c(68L, 69L, 70L, 71L, 78L), -c(68L, 69L, 70L, 71L, 78L)]

# B. resort the variable so they will be more consistent with the class order in the paper

# rearrange the order of the 128 corrubin matrix
corrubinmf <- corrubinmf[index_m_2, index_f_2] 



# ******** -----
## C-1. Saving Rdata ----

outdirectory <- "results"
outfilename <- "corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata"
# outfilename <- sprintf("%s_reg_7.Rdata", depVariable)
save(file=file.path(outdirectory,outfilename), 
     life, impute_f, impute_m,
     resid_impute_f, resid_impute_m,
     corrubinf_133, corrubinm_133, corrubinmf_133,
     corrubinf, corrubinm, corrubinmf)


# Load data
# load("results/corr_chems_heatmap__resid_lipid_creat_adj.Rdata") 

# v2
# load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata") 



# ******** -----
# D. Heatmaps of 133 chems -------------------------------------------------------------



######################################################### #
## \\ D1. within female ----
######################################################### #


library(gplots)

# Color function  # use cp group's common scheme for faster interpretation
heatmapColors <- function(numColors=16) {
    c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
    c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
    c3 <- c(c1,rev(c2));  # rev: reverse element
    return(c3)
}



######################################################### #
# color strip of category
######################################################### #

######################################################### #
# Set color

colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                "#673770", "#D3D93E", "#508578", "#D7C1B1", "#689030") # 13 color


# Create chem class name and color character vector
names(colorCodes) <- names(cl_f)


## create a chemical and chemical category directory vector
cd_f = structure(rep(names(cl_f), times = sapply(cl_f, length)),
                 names = unlist(cl_f, use.names = F))

# put conversion vector cd_f into corrubin order
heatorder_f <- unlist(lapply(colnames(corrubinf), function(x) which(names(cd_f) %in% x)))
heatorder_f_2 <- cd_f[heatorder_f]

# put conversion vector "colorcodes" into heatorder_f_2
heatmap_color_f <- unlist(lapply(heatorder_f_2, function(x) which(names(colorCodes) %in% x)))
heatmap_color_f_2 <- colorCodes[heatmap_color_f]

#
# quartz(height=6, width=8)
quartz(height=6, width=6) # square


## No cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinf,
          main = "Females. Spearman r. \n Lipids and Creatinine adj",
          trace = "none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          dendrogram = "none",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = F,  # row dendrogram should be reordered?
          Colv = F,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_f_2,  # color strip
          ColSideColors = heatmap_color_f_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
            # breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = T)   

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_f_chems_lipids_creat_adj_v2_1.png", type = "png", device = dev.cur(), dpi = 400)
# dev.copy(png, filename="results/heatmap_r_f_chems_lipids_creat_adj.png", width=1024, height=768)
# dev.off()

# Plot legend separately
plot.new()          # create blank plot, given a device is opened and will be cleared.  Not equal to "dev.new()"
par(lend = 1)       # square line or circular etc ends for the color strip indicator in lengend
# par("usr")        # get plot coordinte information of a graphic device
legend(0.0, 0.5,      # e.g. "topright". location of the legend on the heatmap plot
       legend = c("PCBs", "OCPs", "Polybrominated_cpds", "PFASs", "blood_metals", 
                  "Cotinine", "Phytoestrogens", "Phthalates", "Phenols", "Anti_microbial_cpds", 
                  "Paracetamols", "Urine_metals", "Urine_metalloids"), # category labels
       col = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
               "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", 
               "#D7C1B1", "#689030"),  # color key
       lty= 1,             # line style of the color strip
       lwd = 5,            # line width of the color strip
       box.lty = 1,        # line type of the box outline
       box.lwd = 0.5,      # line wide of the box outline
       cex = 0.45           # size of text. The box size will be adj auto.
)

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_f_chems_lipids_creat_adj_v2_1_legend.png", type = "png", device = dev.cur(), dpi = 400)



# quantile, history w desntiy line
quantile(corrubinf[upper.tri(corrubinf)])

hist(corrubinf[upper.tri(corrubinf)],
     main="females, spearman r",
     xlab="r",
     border="blue",
     col="green",
     ylim = c(0, 6),
     # xlim=c(100,700),
     # las=1,  # labels are perpendicular to axis, 0 = parallel, 1 = perpendicular
     # breaks=5,
     prob = TRUE  #probability density expressed through the y-axis instead of the regular frequency.
     )

lines(density(corrubinf[upper.tri(corrubinf)])) #Get a density curve over hist. prob = TRUE MUST be on in hist


##  With cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinf,
          main = "Females. Spearman r. 1-abs(r) \n Lipids and Creatinine adj",
          trace="none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          distfun = function(x) as.dist(1-abs(corrubinf)),  #using 1-abs(r) as the dist matrix
          hclustfun = hclust, # default complete linkage in hclust
          dendrogram = "both",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = TRUE,  # row dendrogram should be reordered?
          Colv = TRUE,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_f_2,  # color strip
          ColSideColors = heatmap_color_f_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          # breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = T)  

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj__cluster.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_f_chems_lipids_creat_adj_v2_1_minus_abs_r.png", type = "png", device = dev.cur(), dpi = 400)




######################################################### #
## \\ D2. within male ----
######################################################### #

######################################################### #
# Set color

colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                "#673770", "#D3D93E", "#508578", "#D7C1B1", "#689030") # 13 color


# Create chem class name and color character vector
names(colorCodes) <- names(cl_m)

## create a chemical and chemical category directory vector
cd_m = structure(rep(names(cl_m), times = sapply(cl_m, length)),
                 names = unlist(cl_m, use.names = F))

# put conversion vector cd_m into corrubin order
heatorder_m <- unlist(lapply(colnames(corrubinm), function(x) which(names(cd_m) %in% x)))
heatorder_m_2 <- cd_m[heatorder_m]

# put conversion vector "colorcodes" into heatorder_m_2
heatmap_color_m <- unlist(lapply(heatorder_m_2, function(x) which(names(colorCodes) %in% x)))
heatmap_color_m_2 <- colorCodes[heatmap_color_m]


# quartz(height=6, width=6)

## No cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinm,
          main = "Males. Spearman r. \n Lipids and Creatinine adj",
          trace = "none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          dendrogram = "none",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = F,  # row dendrogram should be reordered?
          Colv = F,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_m_2,  # color strip
          ColSideColors = heatmap_color_m_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          # breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = T)   

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_m_chems_lipids_creat_adj_v2_1.png", type = "png", device = dev.cur(), dpi = 400)
# dev.copy(png, filename="results/heatmap_r_f_chems_lipids_creat_adj.png", width=1024, height=768)
# dev.off()

# Plot legend separately
plot.new()          # create blank plot, given a device is opened and will be cleared.  Not equal to "dev.new()"
par(lend = 1)       # square line or circular etc ends for the color strip indicator in lengend
# par("usr")        # get plot coordinte information of a graphic device
legend(0.0, 0.5,      # e.g. "topright". location of the legend on the heatmap plot
       legend = c("PCBs", "OCPs", "Polybrominated_cpds", "PFASs", "blood_metals", 
                  "Cotinine", "Phytoestrogens", "Phthalates", "Phenols", "Anti_microbial_cpds", 
                  "Paracetamols", "Urine_metals", "Urine_metalloids"), # category labels
       col = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
               "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", 
               "#D7C1B1", "#689030"),  # color key
       lty= 1,             # line style of the color strip
       lwd = 5,            # line width of the color strip
       box.lty = 1,        # line type of the box outline
       box.lwd = 0.5,      # line wide of the box outline
       cex = 0.45           # size of text. The box size will be adj auto.
)

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_m_chems_lipids_creat_adj_v2_1_legend.png", type = "png", device = dev.cur(), dpi = 400)



# quantile, history w desntiy line
quantile(corrubinm[upper.tri(corrubinm)])

hist(corrubinm[upper.tri(corrubinm)],
     main="males, spearman r",
     xlab="r",
     border="blue",
     col="green",
     ylim = c(0, 6),
     # xlim=c(100,700),
     # las=1,  # labels are perpendicular to axis, 0 = parallel, 1 = perpendicular
     # breaks=5,
     prob = TRUE  #probability density expressed through the y-axis instead of the regular frequency.
)

lines(density(corrubinm[upper.tri(corrubinm)])) #Get a density curve over hist. prob = TRUE MUST be on in hist


##  With cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinm,
          main = "Males. Spearman r. 1-abs(r) \n Lipids and Creatinine adj",
          trace="none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          distfun = function(x) as.dist(1-abs(corrubinm)),  #using 1-abs(r) as the dist matrix
          hclustfun = hclust, # default complete linkage in hclust
          dendrogram = "both",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = TRUE,  # row dendrogram should be reordered?
          Colv = TRUE,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_m_2,  # color strip
          ColSideColors = heatmap_color_m_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          # breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = T)  

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj__cluster.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_m_chems_lipids_creat_adj_v2_1_minus_abs_r.png", type = "png", device = dev.cur(), dpi = 400)



######################################################### #
## \\ D3. within couples ----
######################################################### #


# quartz(height=6, width=6)

## No cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinmf,
          main = "Couples Spearman r. \n Lipids and Creatinine adj",
          trace = "none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          dendrogram = "none",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = F,  # row dendrogram should be reordered?
          Colv = F,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_m_2,  # color strip
          ColSideColors = heatmap_color_f_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = F)  # should be T but not sure why got a warning and 1 v unfilled color key

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_mf_chems_lipids_creat_adj_v2_1.png", type = "png", device = dev.cur(), dpi = 400)
# dev.copy(png, filename="results/heatmap_r_f_chems_lipids_creat_adj.png", width=1024, height=768)
# dev.off()

# Plot legend separately
plot.new()          # create blank plot, given a device is opened and will be cleared.  Not equal to "dev.new()"
par(lend = 1)       # square line or circular etc ends for the color strip indicator in lengend
# par("usr")        # get plot coordinte information of a graphic device
legend(0.0, 0.5,      # e.g. "topright". location of the legend on the heatmap plot
       legend = c("PCBs", "OCPs", "Polybrominated_cpds", "PFASs", "blood_metals", 
                  "Cotinine", "Phytoestrogens", "Phthalates", "Phenols", "Anti_microbial_cpds", 
                  "Paracetamols", "Urine_metals", "Urine_metalloids"), # category labels
       col = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
               "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", 
               "#D7C1B1", "#689030"),  # color key
       lty= 1,             # line style of the color strip
       lwd = 5,            # line width of the color strip
       box.lty = 1,        # line type of the box outline
       box.lwd = 0.5,      # line wide of the box outline
       cex = 0.45           # size of text. The box size will be adj auto.
)

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_mf_chems_lipids_creat_adj_v2_1_legend.png", type = "png", device = dev.cur(), dpi = 400)



# quantile, history w desntiy line
quantile(corrubinmf[upper.tri(corrubinmf)])

hist(corrubinmf[upper.tri(corrubinmf)],
     main="couples, spearman r",
     xlab="r",
     border="blue",
     col="green",
     ylim = c(0, 8),
     # xlim=c(100,700),
     # las=1,  # labels are perpendicular to axis, 0 = parallel, 1 = perpendicular
     # breaks=5,
     prob = TRUE  #probability density expressed through the y-axis instead of the regular frequency.
)

lines(density(corrubinmf[upper.tri(corrubinmf)])) #Get a density curve over hist. prob = TRUE MUST be on in hist


##  With cluster to show couples' corr of each chemical in the diagonal 
heatmap.2(corrubinmf,
          main = "Couples Spearman r. 1-abs(r) \n Lipids and Creatinine adj",
          trace="none", 
          margins = c(4, 4), 
          col = heatmapColors(4),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.23,  # row label size
          cexCol = 0.23,  # column label size
          distfun = function(x) as.dist(1-abs(corrubinmf)),  #using 1-abs(r) as the dist matrix
          hclustfun = hclust, # default complete linkage in hclust
          dendrogram = "both",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = TRUE,  # row dendrogram should be reordered?
          Colv = TRUE,  # column dendrogram should be reordered?
          RowSideColors = heatmap_color_m_2,  # color strip
          ColSideColors = heatmap_color_f_2,  # color strip
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), # force break like this from -1 to 1 with stated interval width
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = F)  # should be T but not sure why got a warning and 1 v unfilled color key

# quartz.save("results/heatmap_r_f_chems_lipids_creat_adj__cluster.pdf", type = "pdf", device = dev.cur(), dpi = 200)
quartz.save("results/heatmap_r_mf_chems_lipids_creat_adj_v2_1_minus_abs_r.png", type = "png", device = dev.cur(), dpi = 400)


