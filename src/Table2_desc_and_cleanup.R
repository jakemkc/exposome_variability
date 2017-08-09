## Oct 17 2016
## Goal: explore the dataset a bit

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L

# install.packages("haven")
library(haven)

life <- read_sas("data/LIFE_exposome.sas7bdat", catalog_file = "data/formats.sas7bcat")
# ?read_sas, read both data file and catalog file. The latter is for recording value labels
impute_m <- read_sas("data/imputedchemicals_male.sas7bdat")
impute_f <- read_sas("data/imputedchemicals_female.sas7bdat")

tlife <- as.data.frame(t(head(life))) # for View()

# ******** -----
# A. Initial exploration ------------------------------------------------------------

class(life) # tibble
str(life) # lots of attributes
mode(life)
dim(life)
    # n = 501. we have 501 pairs of couples, so it's in partial short format
summary(life) # interpret with the cookbook xls file


# To numeric (since all are numbers)
life <- sapply(life, as.numeric)
str(life); class(life); mode(life)

# From matrix numeric to data.frame numeric 
life <- as.data.frame(life, stringsAsFactors = FALSE)

# Skim data
lifehead <- t(head(life))

# Make a better table
sumstat1 <- do.call('rbind', sapply(life, function(x) if(is.numeric(x)) summary(x)))

sum <- sapply(life, summary)
lengthindex <- sapply(sum, length)
sumstat2 <- as.data.frame(do.call(rbind, lapply(sum, "length<-", max(lengthindex))), stringsAsFactors = FALSE)
sumstat2    # v7 is NA column



# ******** -----
# B. Cleanup ----------------------------------------------------- 

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
## \\ B2. "" (blk) & " " (spc) as NA ----
######################################################### #

life <- as.data.frame(apply(life, 2, function(x) gsub("^$|^ $", NA, x)), stringsAsFactors = FALSE)

# attribute all gone and var become chr
str(life); class(life); mode(life)

# to numeric
life <- sapply(life, as.numeric)
str(life); class(life); mode(life)

# from matrix numeric to data.frame numeric 
life <- as.data.frame(life, stringsAsFactors = FALSE)


# ******** -----
# C. Special value char  -----------------------------------------


######################################################### #
# \\ C1. Negative count in chems ----
######################################################### #


############################ #
# \C1.1 Gen Chem index 
############################ #

## Select the index range for chemicals for plots
which(colnames(life)=="pcb028amt_F")    #34
which(colnames(life)=="APAP_m")     # 301


######################################################### #
## \C1.2 How many negative counts per chemical? 

    # How the netgative count range?
    # Answer: ranged from a few to ~400

######################################################### #

negnub <- function(x) {
    table(sign(x))["-1"]  # sign(x), then table the results, then select "-1" entries
}


## Extract those variables with negative values
neg1 <- apply(life[, 34:301], 2, negnub)
neg2 <- is.na(neg1)
neg3 <- neg1[!neg2]  

    # Question: how many chemicals out of 
    # Answer: 83 out of 268 chemicals has neg. numbers in the life df (male+female)

## Plot the number of negative count
par (mar= c(8, 3, 1, 1))  # general parameters for base grpahic
barplot(neg3, las=2, ylim = c(0, 520), cex.names=0.5) 


######################################################### #
# \\ C2. NA count in chems ------------------------------------------------------------

######################################################### #
## \\ C2.1 Question: how many NA counts per chemical?

    # Question: how many chemcials has NA in the life df (male + female)
    # Answer: all 268 chemicals has NA
    
    # Question: hows the NA count range?
    # Answer: from ~10 to ~65
######################################################### #

nanub <- function(x) {
    ab <- table(sign(x), useNA = "always")  
    ab[is.na(names(ab))]
}

na1 <- apply(life[, 34:301], 2, nanub)
na2 <- is.na(na1)
na3 <- na1[!na2]  

# quartz(height=6, width=12)

par (mar= c(6, 3, 1, 1))
barplot(na3, las=2, ylim = c(0, 80), cex.names=0.25) 


# ******** -----
# D. X Y freq imputation on LIFE -------------------------------------------------------------
## equivalent to socialdemographic and lifestyle info in Table 1 of a paper for understanding


######################################################### #
# \\ D1. unique var in life vs imputed ----

diff1 <- setdiff(colnames(life), colnames(impute_f))
diff2 <- setdiff(diff1, colnames(impute_m))

# [1] "Age_m"               "Age_f"               "RaceEthnicity_m"     "DMEDU"              
# [5] "DMHSINCO"            "RaceEthnicity_f"     "DFEDU"               "DFHSINCO"           
# [9] "parity_f"            "parity_m"            "LFSMKNOW"            "LMSMKNOW"           
# [13] "catbmi_m"            "catbmi_f"            "conditionalParity_f" "conditionalParity_m"
# [17] "conceptiondelay"     "INFERTILE"           "LFEXERCS"            "LFSMKCGD"           
# [21] "LFAL12DR"            "LFETOHTY"            "MFHIGHBP"            "MFHGHCHL"           
# [25] "MFDIABTS"            "LMEXERCS"            "LMSKCIGD"            "LMETOHBE"           
# [29] "LMETOHTY"            "MMHIGHBP"            "MMHIGHCH"            "MMDIABET"           
# [33] "amy_mean"            "cort_mean"           "SEVOLUME"            "SPCOUNT"            
# [37] "TOTCNT"              "PERMOT"              "STCRIT"              "WHONORM"            
# [41] "catsevolume"         "catspcount"          "cattotcnt"           "catstcrit"          
# [45] "catwhonorm"          "SCSADFI"             "SCSAHDS"             "catscsadfi"         
# [49] "catscsahds"   

    # 1) basic demo
    # 2) preg outcome and female stress
    # 3) male semen paramters

######################################################### #
# \\ D1. Y freq dist ----

# Answer:
    #                 0   1 NA
    # conceptiondelay 294 137 70
    # INFERTILE       345  57 99
    # catsevolume      50 423 28
    # catspcount       45 428 28
    # cattotcnt        47 426 28
    # catstcrit        18 418 65
    # catwhonorm      208 228 65
    # catscsadfi      425  34 42
    # catscsahds      459  0  42  ** NO "1", leave it out from analysis

# No MI, no simple impute to NA. Use complete.case analysis
######################################################### #


# From the codebook, all Ys (2+7) are binary
    # conceptiondelay
    # INFERTILE
    # catsevolume
    # catspcount
    # cattotcnt
    # catstcrit
    # catwhonorm
    # catscsadfi
    # catscsahds

# category distribution of conceptiondelay
table(life$conceptiondelay, useNA = "always") 
    # 0: 294, 1: 137, NA: 70; total 501

## Extract Ys after checking with colnames(life)
ynames <- c("conceptiondelay", "INFERTILE", "catsevolume", "catspcount", "cattotcnt", "catstcrit", "catwhonorm", "catscsadfi", "catscsahds")
        # ynames <- dput(colnames(life[,c(18, 19, 308:312, 315, 316)]))

# apply table function to each of the extracted columns
yout <- apply(life[,ynames], 2, table, useNA = "always")

# Prep condensed table: Combine list to df
indxy <- sapply(yout, length) # length of each sublist

ytable <- as.data.frame(do.call(rbind, lapply(yout, "length<-", max(indxy))))



######################################################### #
# \\ D2. X adj. covariates & imputation ----

    ## impute covariates here with the most freq values (categorical) or mean (continous), as collaborator suggested

######################################################### #


######################################################### #
# ++ D2.1 Female ----
######################################################### #


######################################################### #
# Continous x1:
    # Age_f
    
    ## Answer: No NA, no zero, no negative, OK
######################################################### #

summary(life$Age_f)
which(is.na((life$Age_f)) %in% TRUE)
plot(density(life$Age_f, na.rm = TRUE)) 


######################################################### #
# Count x2: 
    # LFSMKCGD  ## daily cigarettes at enrollment
    # LFETOHTY  ## # usual drinks per occasion at enrollment
        
        # Parity should be a count based on def and the data (integer, 0-5, 1-5). Be aware if you found sth

    ## Answer: No NA, no negative, OK
    ## As continuous: highly skewed for LFSMKCGD; a bit weird dist for LFETOHTY
######################################################### #

table(life$LFSMKCGD, useNA = "always")
plot(density(life$LFSMKCGD, na.rm = TRUE)) 

table(subset(life, LFSMKCGD >0)$LFSMKCGD)
plot(density(subset(life, LFSMKCGD >0)$LFSMKCGD, na.rm = TRUE)) 

table(life$LFETOHTY, useNA = "always")
plot(density(life$LFETOHTY, na.rm = TRUE)) 



######################################################### #
# Categorical x12:
    
    # "catbmi_f", 
    # "RaceEthnicity_f", 
    # "DFEDU", 
    # "DFHSINCO", 
    # "conditionalParity_f", 
    # "parity_f", 
    # "LFSMKNOW", 
    # "LFEXERCS", 
    # "LFAL12DR", 
    # "MFHIGHBP", 
    # "MFHGHCHL", 
    # "MFDIABTS", 

    # Answer: about half has NA <5, only 1 is NA = 10
        # imputed with most freq 
        # bmi (x3), conditional parity (x3), parity (x6) are not binary level
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
# ++ D2.2 Male ----
######################################################### #

######################################################### #
# Continous x1:
# Age_m

    ## Answer: No NA, no zero, no negative, OK
######################################################### #

summary(life$Age_m)
which(is.na((life$Age_m)) %in% TRUE)
plot(density(life$Age_m, na.rm = TRUE)) 


######################################################### #
# Count x2: 
    # LMSKCIGD
    # LMETOHTY
            
        # Parity should be a count based on def and the data (integer, 0-5, 1-5). Be aware if you found sth

    ## Answer: No NA, no negative, OK
    ## As continuous: highly skewed for LMSKCIGD; a bit weird dist for LMETOHTY
######################################################### #

table(life$LMSKCIGD, useNA = "always")
plot(density(life$LMSKCIGD, na.rm = TRUE)) 

table(subset(life, LMSKCIGD >0)$LMSKCIGD)
plot(density(subset(life, LMSKCIGD >0)$LMSKCIGD, na.rm = TRUE)) 

table(life$LMETOHTY)
plot(density(life$LMETOHTY, na.rm = TRUE)) 



######################################################### #
# Categorical x12:

    # "catbmi_m", 
    # "RaceEthnicity_m", 
    # "DMEDU", 
    # "DMHSINCO", 
    # "conditionalParity_m", 
    # "parity_m", 
    # "LMSMKNOW", 
    # "LMEXERCS", 
    # "LMETOHBE", 
    # "MMDIABET", 
    # "MMHIGHBP", 
    # "MMHIGHCH", 
    
    # Answer: about half has NA <5, only 1 is NA = 10
    # imputed with most freq 
    # bmi (x3), conditional parity (x3), parity (x6) are not binary level
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


# To numeric
life <- sapply(life, as.numeric)
life <- as.data.frame(life, stringsAsFactors = FALSE)


# ******** -----
# ***** LIFE df ready for below ----


# ******** -----
# D. X Y characteristics -------------------------------------------------------------


######################################################### #
# \\ D3. X (chem) dist ----

    # Question: How's the distribution of chemicals?
    # Answer: nearly all are heavily skewed to the right (right tail), for both male and female
            # The imputed one has lots of negative

######################################################### #


quartz(height=6, width=16)
    # quartz(height=6, width=8)

## header + female
boxplot(scale(life[, c(1:33, 34:168)]), pars = list(boxwex = 0.6, outcol="black"))
    # 1:33 is the range for basic demo, female and male
par(cex.axis = 1.2, cex = 0.4, las = 2, fg="black", col="black", mar = c(15, 4, 4, 2) + 0.1)

## header + male
par(cex.axis = 1.2, cex = 0.4, las = 2, fg="black", col="black", mar = c(15, 4, 4, 2) + 0.1)
boxplot(scale(life[, c(1:33, 169:301)]), pars = list(boxwex = 0.6, outcol="black"))



######################################################### #
## \\ D4. X selected chem in table ----
######################################################### #

## Females
## cotinine

continine_f <- list()

for(i in 1:10) {
    tmpf <- subset(impute_f, X_Imputation_ == i)
    cotf <- tmpf$cotamt_F
    continine_f[[i]] <- cotf
}

for(i in 1:10) {
    tmpf <- continine_f[[i]]
    cotf <- log(tmpf+1)
    continine_f[[i]] <- cotf
}

cot_mean_vect_f <- sapply(continine_f, mean)
cot_sd_vect_f <- sapply(continine_f, SD)

#rubin mean (log+1 transformed; log = base 2 in r)
cot_mean_f <- mean(cot_mean_vect_f)
cot_sd_f <- mean(cot_sd_vect_f)

# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
cot_within_imp_var_f <- mean((sapply(continine_f, SE))^2)


# between immputation var (using mean imputation)
cot_bet_imp_var_f <- sum((cot_mean_vect_f - cot_mean_f)^2)/(length(continine_f)-1)

# SE of rubin mean (log+1 transformed)
cot_SE_f <- sqrt(cot_within_imp_var_f + (1 + 1/length(continine_f))*cot_bet_imp_var_f)



## creatinine
creatinine_f <- list()

for(i in 1:10) {
    tmpf <- subset(impute_f, X_Imputation_ == i)
    cretf <- tmpf$CREAMOUNT_F
    creatinine_f[[i]] <- cretf
}

for(i in 1:10) {
    tmpf <- creatinine_f[[i]]
    cretf <- log(tmpf+1)
    creatinine_f[[i]] <- cretf
}

cret_mean_vect_f <- sapply(creatinine_f, mean)
cret_sd_vect_f <- sapply(creatinine_f, sd)

#rubin mean (log+1 transformed; log = base 2 in r)
cret_mean_f <- mean(cret_mean_vect_f)
cret_sd_f <- mean(cret_sd_vect_f)

# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
cret_within_imp_var_f <- mean((sapply(creatinine_f, SE))^2)


# between immputation var (using mean imputation)
cret_bet_imp_var_f <- sum((cret_mean_vect_f - cret_mean_f)^2)/(length(creatinine_f)-1)

# SE of rubin mean (log+1 transformed)
cret_SE_f <- sqrt(cret_within_imp_var_f + (1 + 1/length(creatinine_f))*cret_bet_imp_var_f)


## total lipids

## lipids

lipids_f <- list()

for(i in 1:10) {
    tmpf <- subset(impute_f, X_Imputation_ == i)
    lipids__f <- (1.494*tmpf$cholamt_F)+tmpf$trigamt_F+tmpf$phosamt_F
    lipids_f[[i]] <- lipids__f
}


for(i in 1:10) {
    tmpf <- lipids_f[[i]]
    lipids__f <- log(tmpf+1)
    lipids_f[[i]] <- lipids__f  
}

lip_mean_vect_f <- sapply(lipids_f, mean)
lip_sd_vect_f <- sapply(lipids_f, sd)

#rubin mean (log+1 transformed; log = base 2 in r)
lip_mean_f <- mean(lip_mean_vect_f)
lip_sd_f <- mean(lip_sd_vect_f)


# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
lip_within_imp_var_f <- mean((sapply(lipids_f, SE))^2)


# between immputation var (using mean imputation)
lip_bet_imp_var_f <- sum((lip_mean_vect_f - lip_mean_f)^2)/(length(lipids_f)-1)

# SE of rubin mean (log+1 transformed)
lip_SE_f <- sqrt(lip_within_imp_var_f + (1 + 1/length(lipids_f))*lip_bet_imp_var_f)



## Males
## cotinine

continine_m <- list()

for(i in 1:10) {
    tmpm <- subset(impute_m, X_Imputation_ == i)
    cotm <- tmpm$cotamt_M
    continine_m[[i]] <- cotm
}

for(i in 1:10) {
    tmpm <- continine_m[[i]]
    cotm <- log(tmpm+1)
    continine_m[[i]] <- cotm
}

cot_mean_vect_m <- sapply(continine_m, mean)
cot_sd_vect_m <- sapply(continine_m, sd)


#rubin mean (log+1 transformed; log = base 2 in r)
cot_mean_m <- mean(cot_mean_vect_m)
cot_sd_m <- mean(cot_sd_vect_m)

# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
cot_within_imp_var_m <- mean((sapply(continine_m, SE))^2)


# between immputation var (using mean imputation)
cot_bet_imp_var_m <- sum((cot_mean_vect_m - cot_mean_m)^2)/(length(continine_m)-1)

# SE of rubin mean (log+1 transformed)
cot_SE_m <- sqrt(cot_within_imp_var_m + (1 + 1/length(continine_m))*cot_bet_imp_var_m)


## creatinine

creatinine_m <- list()

for(i in 1:10) {
    tmpm <- subset(impute_m, X_Imputation_ == i)
    cretm <- tmpm$CREAMOUNT_M
    creatinine_m[[i]] <- cretm
}


for(i in 1:10) {
    tmpm <- creatinine_m[[i]]
    cretm <- log(tmpm+1)
    creatinine_m[[i]] <- cretm
}

cret_mean_vect_m <- sapply(creatinine_m, mean)
cret_sd_vect_m <- sapply(creatinine_m, sd)


#rubin mean (log+1 transformed; log = base 2 in r)
cret_mean_m <- mean(cret_mean_vect_m)
cret_sd_m <- mean(cret_sd_vect_m)


# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
cret_within_imp_var_m <- mean((sapply(creatinine_m, SE))^2)


# between immputation var (using mean imputation)
cret_bet_imp_var_m <- sum((cret_mean_vect_m - cret_mean_m)^2)/(length(creatinine_m)-1)

# SE of rubin mean (log+1 transformed)
cret_SE_m <- sqrt(cret_within_imp_var_m + (1 + 1/length(creatinine_m))*cret_bet_imp_var_m)


## total lipids

## lipids

lipids_m <- list()

for(i in 1:10) {
    tmpm <- subset(impute_m, X_Imputation_ == i)
    lipids__m <- (1.494*tmpm$cholamt_M)+tmpm$trigamt_M+tmpm$phosamt_M
    lipids_m[[i]] <- lipids__m
}


for(i in 1:10) {
    tmpm <- lipids_m[[i]]
    lipids__m <- log(tmpm+1)
    lipids_m[[i]] <- lipids__m  
}

lip_mean_vect_m <- sapply(lipids_m, mean)
lip_sd_vect_m <- sapply(lipids_m, sd)

#rubin mean (log+1 transformed; log = base 2 in r)
lip_mean_m <- mean(lip_mean_vect_m)
lip_sd_m <- mean(lip_sd_vect_m)

# within imputation var (sq. SE of each imputation)
SE <- function(x) sd(x)/sqrt(length(x))
lip_within_imp_var_m <- mean((sapply(lipids_m, SE))^2)


# between immputation var (using mean imputation)
lip_bet_imp_var_m <- sum((lip_mean_vect_m - lip_mean_m)^2)/(length(lipids_m)-1)

# SE of rubin mean (log+1 transformed)
lip_SE_m <- sqrt(lip_within_imp_var_m + (1 + 1/length(lipids_m))*lip_bet_imp_var_m)



######################################################### #
## \\ D5. X selected chem p value in table ----
######################################################### #

## data frame prep

x_table_f <- c("Age_f", "LFSMKCGD", "LFETOHTY", "catbmi_f", "RaceEthnicity_f", "DFEDU", "DFHSINCO", 
            "conditionalParity_f", "parity_f", "LFSMKNOW", "LFEXERCS", 
            "LFAL12DR", "MFHIGHBP", "MFHGHCHL", "MFDIABTS")

x_table_m <- c("Age_m", "LMSKCIGD", "LMETOHTY", "catbmi_m", "RaceEthnicity_m", "DMEDU", "DMHSINCO", 
            "conditionalParity_m", "parity_m", "LMSMKNOW", "LMEXERCS", 
            "LMETOHBE", "MMDIABET", "MMHIGHBP", "MMHIGHCH")

dput(which(colnames(life) %in% x_table_f))
dput(which(colnames(life) %in% x_table_m))

tabledf_f <- life[, c(3L, 7L, 8L, 9L, 10L, 12L, 15L, 16L, 20L, 21L, 22L, 23L, 24L, 
               25L, 26L)]
tabledf_f$gender <- "F"

tabledf_m <- life[, c(2L, 4L, 5L, 6L, 11L, 13L, 14L, 17L, 27L, 28L, 29L, 30L, 31L, 
                      32L, 33L)]
tabledf_m$gender <- "M"
colnames(tabledf_m) <- colnames(tabledf_f)

tabledf <- rbind(tabledf_f, tabledf_m)

for(i in 1:15){
    tabledf[, i+1] <- as.factor(tabledf[, i+1])     
}
str(tabledf)
summary(tabledf)


## now get summary p value from LIFE df

t.test(Age_f ~ gender, data = tabledf)
# t = -6.1954, df = 971.8, p-value = 8.575e-10

chiresult <- chisq.test(table(tabledf$gender, tabledf$RaceEthnicity_f))
# X-squared = 0, df = 1, p-value = 1


chisq.test(table(tabledf$gender, tabledf$DFEDU))
#X-squared = 3.8806, df = 1, p-value = 0.04885

chisq.test(table(tabledf$gender, tabledf$DFHSINCO))
#X-squared = 0, df = 1, p-value = 1

chisq.test(table(tabledf$gender, tabledf$LFEXERCS))
#X-squared = 0.41251, df = 1, p-value = 0.5207

chisq.test(table(tabledf$gender, tabledf$LFSMKNOW))
#X-squared = 2.5545, df = 1, p-value = 0.11

chisq.test(table(tabledf$gender, tabledf$LFAL12DR))
#X-squared = 17.547, df = 1, p-value = 2.802e-05

chisq.test(table(tabledf$gender, tabledf$MFDIABTS))
#X-squared = 2.4999, df = 1, p-value = 0.1139

chisq.test(table(tabledf$gender, tabledf$MFHIGHBP))
#X-squared = 14.381, df = 1, p-value = 0.0001493

chisq.test(table(tabledf$gender, tabledf$MFHGHCHL))
#X-squared = 12.358, df = 1, p-value = 0.000439

tabledf$catbmi_f <- as.numeric(tabledf$catbmi_f)
with(tabledf, table(catbmi_f, gender))
wilcox.test(catbmi_f ~ gender, data = tabledf)
# W = 88484, p-value < 2.2e-16

tabledf$LFETOHTY <- as.numeric(tabledf$LFETOHTY)
with(tabledf, table(LFETOHTY, gender))
wilcox.test(LFETOHTY ~ gender, data = tabledf)
#W = 86601, p-value < 2.2e-16


## prep for LFSMKCGD
table(tabledf$LFSMKCGD)
table(tabledf$LFSMKCGD, tabledf$gender)


tabledf$LFSMKCGD <- as.numeric(levels(tabledf$LFSMKCGD))[tabledf$LFSMKCGD]
tabledf_smoke <- subset(tabledf, LFSMKCGD > 0)  # total = 117, match the table in paper

tabledf_smoke$test <- as.numeric(100)

tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >=1 & tabledf_smoke$LFSMKCGD <=3)] <- 1
tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >=4 & tabledf_smoke$LFSMKCGD <=6)] <- 2
tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >=7 & tabledf_smoke$LFSMKCGD <=10)] <- 3
tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >=11 & tabledf_smoke$LFSMKCGD <=15)] <- 4
tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >=16 & tabledf_smoke$LFSMKCGD <=25)] <- 5
tabledf_smoke$test[which(tabledf_smoke$LFSMKCGD >25)] <- 6

with(tabledf_smoke, table(test, gender))
wilcox.test(test ~ gender, data = tabledf_smoke)
# W = 1828, p-value = 0.4988



######################################################### #
## For continous

## 1. the df with missing (mice needs this)
length(which(colnames(life) %in% colnames(impute_f)))
dput(which(colnames(life) %in% colnames(impute_f)))

frm_f <- cbind(impute_f[1:501,1], life[, c(1L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 
                45L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 53L, 54L, 55L, 56L, 57L, 
                58L, 59L, 60L, 61L, 62L, 63L, 64L, 65L, 66L, 67L, 68L, 69L, 70L, 
                71L, 72L, 73L, 74L, 75L, 76L, 77L, 78L, 79L, 80L, 81L, 82L, 83L, 
                84L, 85L, 86L, 87L, 88L, 89L, 90L, 91L, 92L, 93L, 94L, 95L, 96L, 
                97L, 98L, 99L, 100L, 101L, 102L, 103L, 104L, 105L, 106L, 107L, 
                108L, 109L, 110L, 111L, 112L, 113L, 114L, 115L, 116L, 117L, 118L, 
                119L, 120L, 121L, 122L, 123L, 124L, 125L, 126L, 127L, 128L, 129L, 
                130L, 131L, 132L, 133L, 134L, 135L, 136L, 137L, 138L, 139L, 140L, 
                141L, 142L, 143L, 144L, 145L, 146L, 147L, 148L, 149L, 150L, 151L, 
                152L, 153L, 154L, 155L, 156L, 157L, 158L, 159L, 160L, 161L, 162L, 
                163L, 164L, 165L, 166L)])

frm_f$gender <- "F"
frm_f[,1] <- 0
colnames(frm_f) <- c(colnames(impute_f), "gender")


length(which(colnames(life) %in% colnames(impute_m)))
dput(which(colnames(life) %in% colnames(impute_m)))

frm_m <- cbind(impute_m[1:501,1], life[, c(1L, 169L, 170L, 171L, 172L, 173L, 174L, 175L, 176L, 177L, 178L, 
                                      179L, 180L, 181L, 182L, 183L, 184L, 185L, 186L, 187L, 188L, 189L, 
                                      190L, 191L, 192L, 193L, 194L, 195L, 196L, 197L, 198L, 199L, 200L, 
                                      201L, 202L, 203L, 204L, 205L, 206L, 207L, 208L, 209L, 210L, 211L, 
                                      212L, 213L, 214L, 215L, 216L, 217L, 218L, 219L, 220L, 221L, 222L, 
                                      223L, 224L, 225L, 226L, 227L, 228L, 229L, 230L, 231L, 232L, 233L, 
                                      234L, 235L, 236L, 237L, 238L, 239L, 240L, 241L, 242L, 243L, 244L, 
                                      245L, 246L, 247L, 248L, 249L, 250L, 251L, 252L, 253L, 254L, 255L, 
                                      256L, 257L, 258L, 259L, 260L, 261L, 262L, 263L, 264L, 265L, 266L, 
                                      267L, 268L, 269L, 270L, 271L, 272L, 273L, 274L, 275L, 276L, 277L, 
                                      278L, 279L, 280L, 281L, 282L, 283L, 284L, 285L, 286L, 287L, 288L, 
                                      289L, 290L, 291L, 292L, 293L, 294L, 295L, 296L, 297L, 298L, 299L, 
                                      300L, 301L)])

frm_m$gender <- "M"
frm_m[,1] <- 0
colnames(frm_m) <- c(colnames(impute_f), "gender")


## 2. the impute_f and impute_m prep
tabledfimpute_f <- impute_f
# tabledfimpute_f$ID <- NULL
tabledfimpute_f$gender <- "F"

tabledfimpute_m <- impute_m
# tabledfimpute_m$ID <- NULL
tabledfimpute_m$gender <- "M"

colnames(tabledfimpute_m) <- colnames(tabledfimpute_f)

tabledfimputed <- rbind(frm_f, frm_m, tabledfimpute_f, tabledfimpute_m)

str(tabledfimputed)
summary(tabledfimputed)

tabledfimputed$gender <- as.factor(tabledfimputed$gender)



library(dplyr) 
tabledfimputed <- arrange(tabledfimputed, X_Imputation_, ID) 


library(mice)

tabledfimpute_m$ID <- NULL # non unique in each imputation

# mice df ready

tablemice <- as.mids(data = tabledfimputed, .imp=1, .id = NULL)



# mice_cot
mice_cot <- with(data = tablemice, glm(I(log(cotamt_F+1)) ~ gender, family = gaussian))

pool_mice_cot <- pool(mice_cot)
summary(pool_mice_cot)

# est         se        t       df     Pr(>|t|)
# (Intercept) 0.6179708 0.08439865 7.322047 976.9062 5.098144e-13
# gender2     0.6211040 0.12060279 5.149998 879.8501 3.215468e-07

# mice_creat
mice_creat <- with(data = tablemice, glm(I(log(CREAMOUNT_F+1)) ~ gender, family = gaussian))

pool_mice_creat <- pool(mice_creat)
summary(pool_mice_creat)

# est         se         t       df Pr(>|t|)
# (Intercept) 4.2239354 0.03754439 112.50509 442.1328        0
# gender2     0.5313962 0.05283152  10.05832 486.3060        0

## mice_lipids
mice_lip <- with(data = tablemice, glm(I(log((1.494*cholamt_F) + trigamt_F + phosamt_F)+1) ~ gender, family = gaussian))

pool_mice_lip <- pool(mice_lip)
summary(pool_mice_lip)

# est        se          t       df Pr(>|t|)
# (Intercept) 7.4102335 0.0103758 714.184550 963.5174        0
# gender2     0.1462654 0.0146563   9.979697 971.9625        0


