## Jan 25 2017
## Goal: Find p value for spearman r through permutation test

## Initialize the workspace
rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L


# load
load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata") 


# ******** -----
# A. Perm. p for spearman R -------------------------------------------------------------

######################################################### #
## \\ A1. within females ----
######################################################### #

## now my code ## 
nc <- ncol(resid_impute_f)-1

pertime <- 1000 # variable

arrayref <- array(0, dim = c(nc, nc, pertime))
colnames(arrayref) <- colnames(resid_impute_f)[-1]
rownames(arrayref) <- colnames(resid_impute_f)[-1]

corr_per_imp <- list()

# On the lower half matrix
for(h in 1:10){
    for(i in 1:(nc-1)) {            # col, i
        for(j in seq(i+1, nc)){     # row, j
            for (k in 1:pertime){
                tmpfr <- subset(resid_impute_f, X_Imputation_ == h)[2:ncol(resid_impute_f)]
                y <- tmpfr[, i]
                x <- tmpfr[, j]
                x2 <- sample(x, replace = F)
                rs <- cor(x2, y, method = "spearman")
                arrayref[j, i, k] <- rs
                corr_per_imp[[h]] <- arrayref
            }
        } 
    }        
}


## combine array in List to a single array
arrayperm <- array(unlist(corr_per_imp), dim = c(nc, nc, pertime*h))


# rubin spearman r
rubin_r <- corrubinf

p_perm_spear_r <- matrix(NA, nrow = nc, ncol = nc)
colnames(p_perm_spear_r) <- colnames(resid_impute_f)[-1]
rownames(p_perm_spear_r) <- colnames(resid_impute_f)[-1]

# loop to get permuted p value
for(i in 1:(nc-1)){            # col, i
    for(j in seq(i+1, nc)){     # row, j
        rub_r <- rubin_r[j, i]  # mind it's j i and i j
        perm_r <- arrayperm[j, i, ]
        p_spearman_r <- mean(append(abs(perm_r) > abs(rub_r), 1))  
        p_perm_spear_r[j, i] <- p_spearman_r
    }
}

p_perm_spear_f <- p_perm_spear_r

save(file = file.path("results", "perm_p_value_f_resid_lipid_creat_adj_v2.Rdata"), 
     life, impute_f, impute_m,
     resid_impute_f, resid_impute_m,
     p_perm_spear_f,
     corrubinf, corrubinm, corrubinmf)



######################################################### #
## \\ A2. within males ----
######################################################### #


## now my code ## 
nc <- ncol(resid_impute_m)-1

pertime <- 1000 # variable

arrayref <- array(0, dim = c(nc, nc, pertime))
colnames(arrayref) <- colnames(resid_impute_m)[-1]
rownames(arrayref) <- colnames(resid_impute_m)[-1]

corr_per_imp <- list()

# On the lower half matrix
for(h in 1:10){
    for(i in 1:(nc-1)) {            # col, i
        for(j in seq(i+1, nc)){     # row, j
            for (k in 1:pertime){
                tmpfr <- subset(resid_impute_m, X_Imputation_ == h)[2:ncol(resid_impute_m)]
                y <- tmpfr[, i]
                x <- tmpfr[, j]
                x2 <- sample(x, replace = F)
                rs <- cor(x2, y, method = "spearman")
                arrayref[j, i, k] <- rs
                corr_per_imp[[h]] <- arrayref
            }
        } 
    }        
}


## combine array in List to a single array
arrayperm <- array(unlist(corr_per_imp), dim = c(nc, nc, pertime*h))


# rubin spearman r
rubin_r <- corrubinm

p_perm_spear_r <- matrix(NA, nrow = nc, ncol = nc)
colnames(p_perm_spear_r) <- colnames(resid_impute_m)[-1]
rownames(p_perm_spear_r) <- colnames(resid_impute_m)[-1]

# loop to get permuted p value
for(i in 1:(nc-1)){            # col, i
    for(j in seq(i+1, nc)){     # row, j
        rub_r <- rubin_r[j, i]  # mind it's j i and i j
        perm_r <- arrayperm[j, i, ]
        p_spearman_r <- mean(append(abs(perm_r) > abs(rub_r), 1))
        p_perm_spear_r[j, i] <- p_spearman_r
    }
}

p_perm_spear_m <- p_perm_spear_r

save(file = file.path("results", "perm_p_value_m_resid_lipid_creat_adj_v2.Rdata"), 
     life, impute_f, impute_m,
     resid_impute_m, resid_impute_m,
     p_perm_spear_m,
     corrubinf, corrubinm, corrubinmf)


######################################################### #
## \\ A3. within couples ----
######################################################### #

## now my code ## 
nc <- ncol(resid_impute_f)-1

pertime <- 1000 # variable

arrayref <- array(0, dim = c(nc, nc, pertime))
colnames(arrayref) <- colnames(resid_impute_f)[-1]
rownames(arrayref) <- colnames(resid_impute_m)[-1]

corr_per_imp <- list()

# On the whole matrix
for(h in 1:10){
    for(i in 1:nc) {            # col, i
        for(j in 1:nc){     # row, j
            for (k in 1:pertime){
                tmpfr1 <- subset(resid_impute_f, X_Imputation_ == h)[2:ncol(resid_impute_f)]
                tmpfr2 <- subset(resid_impute_m, X_Imputation_ == h)[2:ncol(resid_impute_m)]
                y <- tmpfr1[, i]
                x <- tmpfr2[, j]
                x2 <- sample(x, replace = F)
                rs <- cor(x2, y, method = "spearman")
                arrayref[j, i, k] <- rs
                corr_per_imp[[h]] <- arrayref
            }
        } 
    }        
}


## combine array in List to a single array
arrayperm <- array(unlist(corr_per_imp), dim = c(nc, nc, pertime*h))


# rubin spearman r (rowname is males for couples, same for corr_per_imp. A match)
rubin_r <- corrubinmf

p_perm_spear_r <- matrix(NA, nrow = nc, ncol = nc)
colnames(p_perm_spear_r) <- colnames(resid_impute_f)[-1]
rownames(p_perm_spear_r) <- colnames(resid_impute_m)[-1]

# loop to get permuted p value
for(i in 1:nc){            # col, i
    for(j in 1:nc){     # row, j
        rub_r <- rubin_r[j, i]  # mind it's j i and i j
        perm_r <- arrayperm[j, i, ]
        p_spearman_r <- mean(append(abs(perm_r) > abs(rub_r), 1))
        p_perm_spear_r[j, i] <- p_spearman_r
    }
}

p_perm_spear_mf <- p_perm_spear_r

save(file = file.path("results", "perm_p_value_mf_resid_lipid_creat_adj_v2.Rdata"), 
     life, impute_f, impute_m,
     resid_impute_f, resid_impute_m,
     p_perm_spear_mf,
     corrubinf, corrubinm, corrubinmf)


