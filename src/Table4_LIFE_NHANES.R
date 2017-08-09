# May 2 2017
# Goal: concordance between NHANES and LIFE


rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=4500) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings


# ******** -----
## A. Load Rdata ----
load("/Users/jake/Documents/2. LIFE/data_chirag_NHANES/corrdata6/entire_spearman_2003-2004.Rdata")
load("/Users/jake/Documents/2. LIFE/data_chirag_NHANES/corrdata6/bigtable_correlation_nhanes.Rdata") # load var index
load("results/corr_concordance_spearman_lipids_creat_adj.Rdata") 

# ******** -----
## B. NHANES ----
library(dplyr)
library(magrittr)

# summary table
filter(metaTable, series == "2003-2004") %>% count(var_desc_ewas) %>% print(n=50) # print(tbl_df(df), n=40) 
filter(metaTable, series == "2003-2004") %>% count(var_desc_ewas_sub) %>% print(n=50) # table(metaTable$var_desc_ewas_sub)
sum(is.na(corrDataGather$corr))/nrow(corrDataGather)*100 # % NA in corr

# export the list
N0304 <- filter(metaTable, series == "2003-2004")  

write.csv(N0304, file = "results/n0304.csv")

# import matching csv
n0304_2 <- read.csv("results/n0304_LIFE.csv", stringsAsFactors = FALSE, header = TRUE, na.strings=c(""," ","NA"))
n0304_2 %>% filter(!is.na(LIFE_desc)) %>% dim() # row number correct, 101 match out of 128
n0304_2 <- n0304_2 %>% filter(!is.na(LIFE_desc)) %>% select(series:LIFE_desc)

n0304_mergev1 <- n0304_2 %>% select(var:LIFE_desc)
n0304_mergev2 <- n0304_2 %>% select(var:LIFE_desc)

# var1 in corrdata joining n0304_merg
colnames(n0304_mergev1) <- sprintf("%s_v1", colnames(n0304_mergev1))
names(n0304_mergev1)[names(n0304_mergev1) == "var_v1"] <- "var1"  # rename var

# var2 in corrdata joining n0304_merg
colnames(n0304_mergev2) <- sprintf("%s_v2", colnames(n0304_mergev2))
names(n0304_mergev2)[names(n0304_mergev2) == "var_v2"] <- "var2"


# left join corrdata n0304_2
merge0304 <- left_join(corrDataGather, n0304_mergev1, by = "var1")
    # merge0304 %>% count(LIFE_desc_v1) %>% print(n=50) # print(tbl_df(df), n=40) 
merge0304 <- left_join(merge0304, n0304_mergev2, by = "var2")

# extract rows with both LIFE_desc_v1 & LIFE_desc_v2
final_0304 <- merge0304 %>% filter(!is.na(LIFE_desc_v1) & !is.na(LIFE_desc_v2))
final_0304 <- final_0304 %>% select(corr, LIFE_desc_v1, LIFE_desc_v2)

# create and rename LIFE pairs
final_0304$LIFE_desc_v1 <- gsub("_f$", "", final_0304$LIFE_desc_v1, ignore.case = T)
final_0304$LIFE_desc_v2 <- gsub("_f$", "", final_0304$LIFE_desc_v2, ignore.case = T)
final_0304$pairs <- sprintf("%s-%s", final_0304$LIFE_desc_v1, final_0304$LIFE_desc_v2)

final_0304 <- final_0304 %>% filter(!is.na(corr)) # remove NA
    # final_0304 %>% count(LIFE_desc_v1) %>% print(n=200) # how many chemicals

# ******** -----
## C. LIFE ----

## Females and Males
# rename LIFE matching var
names(spearF)[names(spearF) == "col_row_pair_2"] <- "pairs"  # rename var
names(spearM)[names(spearM) == "col_row_pair_2"] <- "pairs"

# Left join
spearF <- left_join(spearF, final_0304, by = "pairs")
spearF <- spearF %>% filter(!is.na(corr)) # remove NA

spearM <- left_join(spearM, final_0304, by = "pairs")
spearM <- spearM %>% filter(!is.na(corr)) # remove NA


# ******** -----
## D. concordance ----

# dataframe ready for concordance
spearF$mval <- spearM$val  # val is female spearman; corr is NHANES; mval is male spearman

# concordance: female to NHANES
spearF %>% group_by(categlab) %>% summarise(pearson = cor(corr, val), n = n()) # n() is a function in dplyr
    # spearF %>% group_by(categlab) %>% do(mod1 = lm(corr ~ val, data = .)) %>% mutate(r2 = summary(mod1)$r.squared) %>% select(-mod1)  # same
spearF %>% summarise(pearson = cor(corr, val))
    # spearF %>% filter(!(categlab == "Heterogenous")) %>% summarise(pearson = cor(corr, val), n = n())

# A tibble: 10 × 2
# categlab    pearson
# <chr>      <dbl>
# 1         Blood_metals 0.06735253
# 3                 OCPs 0.68849865
# 4                 PCBs 0.87681741
# 5                PFASs 0.31195057
# 6           Phthalates 0.89606761
# 7        Phytoestogens 0.97971074
# 8  Polybrominated_cpds 0.75898311
# 9     Urine_metalloids 0.40921243
# 10        Urine_metals 0.82110286
# total: 0.877945
# total without hetero: 0.8407948

# concordance: male to NHANES
spearF %>% group_by(categlab) %>% summarise(pearson = cor(corr, mval), n = n()) # n() is a function in dplyr)
spearF %>% summarise(pearson = cor(corr, mval))
    # spearF %>% filter(!(categlab == "Heterogenous")) %>% summarise(pearson = cor(corr, mval), n = n())

# categlab    pearson
# <chr>      <dbl>
# 1         Blood_metals  0.1329664
# 3                 OCPs  0.5088518
# 4                 PCBs  0.8821268
# 5                PFASs  0.2925937
# 6           Phthalates  0.8881209
# 7        Phytoestogens  0.9709152
# 8  Polybrominated_cpds  0.7363123
# 9     Urine_metalloids -0.7773426
# 10        Urine_metals  0.7831219
# total: 0.865469
# total without hetero: 0.8356977


## For Couple
load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata") 

library(tidyr)
spearMF <- as.data.frame(corrubinmf) %>% cbind(row = rownames(corrubinmf)) %>% gather(col, corr, PCB_28_f:tellurium_f) # "gather" dataframe only; "melt" matrix OK
    # melt(corrubinmf) # one step

# Create pair ID
spearMF$col2 <- gsub("_f$", "", spearMF$col, ignore.case = T)
spearMF$row2 <- gsub("_m$", "", spearMF$row, ignore.case = T)

spearMF$pairs <- sprintf("%s-%s", spearMF$col2, spearMF$row2)

names(spearMF)[names(spearMF) == "corr"] <- "MF_corr"  # rename var


# Left join
spearMF <- left_join(spearMF, final_0304, by = "pairs")  
    # Cp took a full matrix, likely wo diagnoal, and thus I have the rev pair match alredy for couple
    # final_0304 %>% filter(pairs == "PCB_28-PCB_44")
    # final_0304 %>% filter(pairs == "PCB_44-PCB_28")

spearMF <- spearMF %>% filter(!is.na(corr)) # remove NA

# category list
cl_f2 = list(
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


# remove the _f in the cl_f2
for(i in 1:length(cl_f2)){
    cl_f2[[i]] <- gsub("_f$", "", cl_f2[[i]], ignore.case = T)
    # print(cl_m[[i]])
}


# create category
spearMF$colcateg <- "n____"
spearMF$rowcateg <- "n____"

for(i in 1:length(cl_f2)) {
    mtchC <- which(spearMF$col2 %in% cl_f2[[i]])
    spearMF$colcateg[mtchC] <- names(cl_f2)[i]
    mtchR <- which(spearMF$row2 %in% cl_f2[[i]])
    spearMF$rowcateg[mtchR] <- names(cl_f2)[i]
}

spearMF$categlab <- "Heterogenous"
replistMF <- spearMF$colcateg == spearMF$rowcateg 
spearMF$categlab[replistMF] <- spearMF$colcateg[replistMF]


# concordance: couples to NHANES
spearMF %>% group_by(categlab) %>% summarise(pearson = cor(MF_corr, corr), n = n()) # n() is a function in dplyr)
spearMF %>% summarise(pearson = cor(MF_corr, corr))
    # spearMF %>% filter(!(categlab == "Heterogenous")) %>% summarise(pearson = cor(MF_corr, corr), n = n())

# # A tibble: 10 × 2
# categlab     pearson
# <chr>       <dbl>
#     1         Blood_metals -0.04317300
# 3                 OCPs  0.24545134
# 4                 PCBs  0.77475019
# 5                PFASs  0.32322229
# 6           Phthalates  0.34217362
# 7       Phytoestrogens  0.85753903
# 8  Polybrominated_cpds  0.76613302
# 9     Urine_metalloids  0.33848620
# 10        Urine_metals -0.01490723
# total: 0.7256615
# total wo hetergroup 0.6665564


# ******** -----
## E. concordance plots ----

# females vs NHANES; all groups
library(ggplot2)
quartz(height=6, width=8)

p <- ggplot(data = spearF, aes(x = val, y = corr))
p <- p + geom_point(aes(color = categlab), shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. Within female spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_f_all.png", scale=1, dpi=400)

# by group
quartz(height=10, width=8)

p <- ggplot(data = spearF, aes(x = val, y = corr))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. Within female spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
# p <- p + facet_grid( ~ categlab)
p <- p + facet_wrap( ~ categlab, nrow = 4)
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm) + ylim(-1, 1)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_f_groups.png", scale=1, dpi=400)



# males vs NHANES; all groups
quartz(height=6, width=8)

library(ggplot2)
p <- ggplot(data = spearF, aes(x = mval, y = corr))
p <- p + geom_point(aes(color = categlab), shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. Within male spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_m_all.png", scale=1, dpi=400)

# by group
quartz(height=10, width=8)

p <- ggplot(data = spearF, aes(x = mval, y = corr))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. Within male spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
# p <- p + facet_grid( ~ categlab)
p <- p + facet_wrap( ~ categlab, nrow = 4)
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm) + ylim(-1, 1)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_m_groups.png", scale=1, dpi=400)



# couples vs NHANES; all groups
quartz(height=6, width=8)

library(ggplot2)
p <- ggplot(data = spearMF, aes(x = MF_corr, y = corr))
p <- p + geom_point(aes(color = categlab), shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. between couple spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_mf_all.png", scale=1, dpi=400)

# by group
quartz(height=10, width=8)

p <- ggplot(data = spearMF, aes(x = MF_corr, y = corr))
p <- p + geom_point(shape = 1, alpha = 1/2)
p <- p + xlab("LIFE. between couple spearman corr coef") + 
    ylab("NHANES spearman corr coef") +
    ggtitle("Concordance")
# p <- p + facet_grid( ~ categlab)
p <- p + facet_wrap( ~ categlab, nrow = 4)
p <- p + geom_abline(slope = 1, intercept = 0, color = "red")
p <- p + geom_smooth(method = lm) + ylim(-1, 1)
p

ggsave("results/scatter_r_concordance_NHANES_LIFE_mf_groups.png", scale=1, dpi=400)

