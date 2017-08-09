## Jan 4 2017
## Goal: create exposome globe

## Initialize the workspace
rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
# quartz(height=6, width=8)

# Load data (Total lipid and creatinine adjusted spearman r)
load("results/corr_chems_heatmap__resid_lipid_creat_adj_v2.Rdata")

# ******** -----
# A. Create lists  -------------------------------------------

## Create a full chem category and chemical list

cl = list(
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
        "tellurium_f"),
    #PCBs
    PCBs. = c(
        "PCB_28_m", 
        "PCB_44_m", 
        "PCB_49_m", 
        "PCB_52_m", 
        "PCB_66_m", 
        "PCB_74_m", 
        "PCB_87_m", 
        "PCB_99_m", 
        "PCB_101_m", 
        "PCB_105_m", 
        "PCB_110_m", 
        "PCB_118_m", 
        "PCB_114_m", 
        "PCB_128_m", 
        "PCB_138_m", 
        "PCB_146_m", 
        "PCB_149_m", 
        "PCB_151_m", 
        "PCB_153_m", 
        "PCB_156_m", 
        "PCB_157_m", 
        "PCB_167_m", 
        "PCB_170_m", 
        "PCB_172_m", 
        "PCB_177_m", 
        "PCB_178_m", 
        "PCB_180_m", 
        "PCB_183_m", 
        "PCB_187_m", 
        "PCB_189_m", 
        "PCB_194_m", 
        "PCB_195_m", 
        "PCB_196_m", 
        "PCB_201_m", 
        "PCB_206_m", 
        "PCB_209_m" 
    ),
    #OCPs
    OCPs. = c(
        "HCB_m",
        "b_HCB_m",
        "g_HCB_m",
        "op_DDT_m",
        "pp_DDE_m",
        "pp_DDT_m",
        "oxychlordane_m",
        "tr_nonachlor_m",
        "mirex_m"),
    #PBC
    Polybrominated_cpds. = c(
        "BB_153_m", 
        "BDE_17_m", 
        "BDE_28_m", 
        "BDE_47_m", 
        "BDE_66_m", 
        "BDE_85_m", 
        "BDE_99_m", 
        "BDE_100_m", 
        "BDE_153_m", 
        "BDE_154_m", 
        "BDE_183_m"),
    #PFASs
    PFASs. = c(
        "Et_PFOSA_AcOH_m", 
        "Me_PFOSA_AcOH_m", 
        "PFDeA_m", 
        "PFNA_m", 
        "PFOSA_m", 
        "PFOS_m", 
        "PFOA_m"),
    #blood metals
    Blood_metals. = c(
        "blood_Cd_m", 
        "blood_Pb_m", 
        "blood_Hg_m"),
    #cotinine
    Cotinine. = c(
        "cotinine_m"),
    #phytoestrogens
    Phytoestrogens. = c(
        "genistein_m", 
        "daidzein_m", 
        "O_DMA_m", 
        "equol_m", 
        "enterodiol_m", 
        "enterolactone_m"),
    #phthalates
    Phthalates. = c(
        "mMP_m", 
        "mEP_m", 
        "mCPP_m", 
        "mBP_m", 
        "miBP_m", 
        "mECPP_m", 
        "mCMHP_m", 
        "mEHHP_m", 
        "mEOHP_m", 
        "mCHP_m", 
        "mBzP_m", 
        "mEHP_m", 
        "mOP_m", 
        "mNP_m"),
    #phenols
    Phenols. = c(
        "BPA_m", 
        "2_OH_4MeO_BP_m", 
        "4_OH_BP_m", 
        "24_OH_BP_m", 
        "22_OH_4MeO_BP_m", 
        "2244_OH_BP_m"),
    #anti microbial
    Anti_microbial_cpds. = c(
        "MP_m", 
        "EP_m", 
        "PP_m", 
        "BP_m", 
        "BzP_m", 
        "HP_m", 
        "4_HB_m", 
        "34_DHB_m", 
        "OH_Me_P_m", 
        "OH_Et_P_m", 
        "TCS_m", 
        "TCC_m"),
    #paracetamol
    Paracetamols. = c(
        "paracetamol_m", 
        "4_aminophenol_m"),
    #urine metals
    Urine_metals. = c(
        "manganese_m", 
        "chromium_m", 
        "beryllium_m", 
        "cobalt_m", 
        "molybdenum_m", 
        "cadmium_m", 
        "tin_m", 
        "caesium_m", 
        "barium_m", 
        "nickel_m", 
        "copper_m", 
        "zinc_m", 
        "tungsten_m", 
        "platinum_m", 
        "thallium_m", 
        "lead_m", 
        "uranium_m"),
    #urine metalloids
    Urine_metalloids. = c(
        "selenium_m", 
        "arsenic_m", 
        "antimony_m", 
        "tellurium_m")
)



# Make male and female symmetric in the plot: 
# reverse male chem cateogry order
    
ordertemp <- names(cl)
ordertemp[14:length(ordertemp)] <- rev(ordertemp[14:length(ordertemp)])
cl <- cl[ordertemp]


for(i in 14:26) {
    cl[[i]] <- rev(cl[[i]])
}



## create a chemical and chemical category directory vector. unlist(cl, use.names = F) can avoid extract funny name from cd
cd = structure(rep(names(cl), times = sapply(cl, length)),
               names = unlist(cl))


# number of chem groups
n_group = length(cl)


# color x13
colorCodes <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1")

# B. Chord diagram  -------------------------------------------

library(circlize)
circos.clear()
circos.par(
    gap.degree = c(rep(2, n_group/2-1), 10, rep(2, n_group/2-1), 10), 
    cell.padding = c(0, 0, 0, 0), 
    start.degree = 84,
    canvas.xlim = c(-1.5, 1.5),  # 1.2
    canvas.ylim = c(-1.5, 1.5)
)

circos.initialize(factors = names(cl), xlim = cbind(c(rep(0, n_group)), sapply(cl, length)))

circos.trackPlotRegion(
    ylim = c(0, 1), 
    bg.border = NA,
    bg.col = c(colorCodes, rev(colorCodes)),
    track.height = 0.04,
    panel.fun = function(x, y) {
        nm = get.cell.meta.data("sector.index")
        r = cl[[nm]]
        n = length(r)
        circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(1, n), lwd = 0.5)
        circos.text(n/2, 3.7, nm, adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = 0.65) 
    }
)




## color
col_fun = colorRamp2(c(-0.75, 0, 0.75), c("darkgreen", "white", "red"))  # 0.75 range

# C. links for within couples  -------------------------------------------
mat_mf = as.matrix(corrubinmf)

mat_mf[mat_mf < 0.25 & mat_mf > -0.25] <- NA  # set NA threshold

n = nrow(mat_mf)
rn = rownames(mat_mf)
cn = colnames(mat_mf)


for(i in 1:n) {
    for(j in 1:n) {
        g1 = cd[rn[i]]
        g2 = cd[cn[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rn[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == cn[j])
        if (is.na(mat_mf[i, j])) {
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_mf[i, j]*2, h = 0.6, # h = 0.6,
                        col = col_fun(mat_mf[i, j]) )#, rou1 = 0.8, rou2 = 0.8) # rou1 = 0.5, rou2 = 0.5 are optional to shrink
        }
    }
}

# circos.clear()

# D. links for within females  -------------------------------------------

mat_f = as.matrix(corrubinf) 

mat_f[mat_f < 0.25 & mat_f > -0.25] <- NA  # set NA threshold

nf = nrow(mat_f)
rnf = rownames(mat_f)

for(i in 1:(nf-1)) {
    for(j in seq(i+1, nf)) {  # rev(seq(i+1, nf))
        g1 = cd[rnf[i]]
        g2 = cd[rnf[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rnf[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == rnf[j])
        if (is.na(mat_f[i, j])) {
        } else if (g1 == g2) {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*2, h = (abs(j-i)/(nf-1)*-0.35), # h = (abs(j-i)/(nf-1)*0.68 + 0.02)
                        col = col_fun(mat_f[i, j]), rou1 = 1.01, rou2 = 1.01)
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*2, h = (abs(j-i)/(nf-1)*0.45 + 0.01), # *0.7 is default # h = (abs(j-i)/(nf-1)*0.68 + 0.02)
                        col = col_fun(mat_f[i, j]))
        }
    }
    
}
# circos.clear()

# E. links for within males  -------------------------------------------

mat_m = as.matrix(corrubinm) 

mat_m[mat_m < 0.25 & mat_m > -0.25] <- NA  # set NA threshold

nm = nrow(mat_m)
rnm = rownames(mat_m)

for(i in 1:(nm-1)) {
    for(j in seq(i+1, nm)) {
        g1 = cd[rnm[i]]
        g2 = cd[rnm[j]]
        r1 = cd[cd == g1]
        k1 = which(names(r1) == rnm[i])
        r2 = cd[cd == g2]
        k2 = which(names(r2) == rnm[j])
        if (is.na(mat_m[i, j])) {
        } else if (g1 == g2) {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*-0.35), # h = (abs(j-i)/(nm-1)*0.68 + 0.02)
                        col = col_fun(mat_m[i, j]), rou1 = 1.01, rou2 = 1.01)
        } else {
            circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*2, h = (abs(j-i)/(nm-1)*0.45 + 0.01),
                        col = col_fun(mat_m[i, j]))
        }
    }
}
circos.clear()

quartz.save("results/chord_r_mf_chems_lipids_creat_adj_v2.png", type = "png", device = dev.cur(), dpi = 400)


