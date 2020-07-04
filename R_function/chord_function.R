


chord_1_sample <- function(tracksample1, correlationsample1, trackColor,
                           # Track option
                           trackFontSize = 0.65,
                           textDistance = 3.7,
                           # Correlation matrix options  
                           corrLineColor = c("darkgreen", "white", "red"), 
                           corrLineColorRampLv = c(-0.75, 0, 0.75),
                           corrNALv = 0.25,
                           # Line options
                           LineWidthMultipler = 2,
                           WithinCorrHeightMultiplier_inclass = -0.35, 
                           WithinCorrHeightMultiplier_betweenclass = 0.45, 
                           WithinCorrHeightOffset_betweenclass = 0.01) {
    
    # ******** -----
    # A. Track Vector  -------------------------------------------
    
    ## Format vector for track
    cl <- tracksample1
    
    ## create chemical class and chemical names vector
    cd = structure(rep(names(cl), times = sapply(cl, length)),
                   names = unlist(cl))
    
    # number of chem groups
    n_group = length(cl)
    
    
    
    # B. Chord diagram  -------------------------------------------
    
    ## Prepare track figure
    library(circlize)
    circos.clear()
    circos.par(
        gap.degree = 3, 
        cell.padding = c(0, 0, 0, 0), 
        start.degree = 84,
        canvas.xlim = c(-1.5, 1.5),
        canvas.ylim = c(-1.5, 1.5)
    )
    
    circos.initialize(factors = names(cl), xlim = cbind(c(rep(0, n_group)), sapply(cl, length)))
    
    circos.trackPlotRegion(
        ylim = c(0, 1), 
        bg.border = NA,
        bg.col = c(trackColor),
        track.height = 0.04,
        panel.fun = function(x, y) {
            nm = get.cell.meta.data("sector.index")
            r = cl[[nm]]
            n = length(r)
            circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(1, n), lwd = 0.5)
            circos.text(n/2, textDistance, nm, adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = trackFontSize) 
        }
    )
    
    
    
    ## Draw correlation lines
    
    ## Color for lines
    col_fun = colorRamp2(corrLineColorRampLv, corrLineColor)
    
    
    # D. links for within females  -------------------------------------------
    
    mat_f = as.matrix(correlationsample1) 
    
    mat_f[mat_f < corrNALv & mat_f > -corrNALv] <- NA 
    
    nf = nrow(mat_f)
    rnf = rownames(mat_f)
    
    for(i in 1:(nf-1)) {
        for(j in seq(i+1, nf)) {  
            g1 = cd[rnf[i]]
            g2 = cd[rnf[j]]
            r1 = cd[cd == g1]
            k1 = which(names(r1) == rnf[i])
            r2 = cd[cd == g2]
            k2 = which(names(r2) == rnf[j])
            if (is.na(mat_f[i, j])) {
            } else if (g1 == g2) {
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*LineWidthMultipler, h = (abs(j-i)/(nf-1)*WithinCorrHeightMultiplier_inclass),
                            col = col_fun(mat_f[i, j]), rou1 = 1.01, rou2 = 1.01)
            } else {
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*LineWidthMultipler, h = (abs(j-i)/(nf-1)*WithinCorrHeightMultiplier_betweenclass + WithinCorrHeightOffset_betweenclass),
                            col = col_fun(mat_f[i, j]))
            }
        }
        
    }
}


chord_2_samples <- function(tracksample1, tracksample2, correlationsample1, correlationsample2, correlationsample1_2, trackColor,
                          # Track option
                          trackFontSize = 0.65,
                          textDistance = 3.7,
                          # Correlation matrix options  
                          corrLineColor = c("darkgreen", "white", "red"), 
                          corrLineColorRampLv = c(-0.75, 0, 0.75),
                          corrNALv = 0.25,
                          # Line options
                          LineWidthMultipler = 2,
                          betweenCorrLineHeight = 0.6,
                          WithinCorrHeightMultiplier_inclass = -0.35, 
                          WithinCorrHeightMultiplier_betweenclass = 0.45, 
                          WithinCorrHeightOffset_betweenclass = 0.01) {

    # ******** -----
    # A. Track Vector  -------------------------------------------
    
    ## Format vector for track
    #' Track reads in clockwise 
    #' reverse male's list
    
    cl <- c(tracksample1, tracksample2)
    ordertemp <- names(cl)
    ordertemp[(length(tracksample1)+1):length(ordertemp)] <- rev(ordertemp[(length(tracksample1)+1):length(ordertemp)])
    cl <- cl[ordertemp]
    
    #  male chemical order within each cateogry
    for(i in (length(tracksample1)+1):length(cl)) {
        cl[[i]] <- rev(cl[[i]])
    }
    
    
    ## create chemical class and chemical names vector
    cd = structure(rep(names(cl), times = sapply(cl, length)),
                   names = unlist(cl))
    
    # number of chem groups
    n_group = length(cl)
    
    
    
    # B. Chord diagram  -------------------------------------------
    
    ## Prepare track figure
    library(circlize)
    circos.clear()
    circos.par(
        gap.degree = c(rep(2, n_group/2-1), 10, rep(2, n_group/2-1), 10), 
        cell.padding = c(0, 0, 0, 0), 
        start.degree = 84,
        canvas.xlim = c(-1.5, 1.5),
        canvas.ylim = c(-1.5, 1.5)
    )
    
    circos.initialize(factors = names(cl), xlim = cbind(c(rep(0, n_group)), sapply(cl, length)))
    
    circos.trackPlotRegion(
        ylim = c(0, 1), 
        bg.border = NA,
        bg.col = c(trackColor, rev(trackColor)),
        track.height = 0.04,
        panel.fun = function(x, y) {
            nm = get.cell.meta.data("sector.index")
            r = cl[[nm]]
            n = length(r)
            circos.rect(seq(0, n-1), rep(0, n), 1:n, rep(1, n), lwd = 0.5)
            circos.text(n/2, textDistance, nm, adj = c(0, 1), facing = c("clockwise"), niceFacing = TRUE, cex = trackFontSize) 
        }
    )
    
    
    
    ## Draw correlation lines
    
    ## Color for lines
    col_fun = colorRamp2(corrLineColorRampLv, corrLineColor)
    
    # C. links for within couples  -------------------------------------------
    mat_mf = as.matrix(correlationsample1_2)
    
    mat_mf[mat_mf < corrNALv & mat_mf > -corrNALv] <- NA
    
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
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_mf[i, j]*LineWidthMultipler, h = betweenCorrLineHeight,
                            col = col_fun(mat_mf[i, j]) )
            }
        }
    }
    
    
    # D. links for within females  -------------------------------------------
    
    mat_f = as.matrix(correlationsample1) 
    
    mat_f[mat_f < corrNALv & mat_f > -corrNALv] <- NA 
    
    nf = nrow(mat_f)
    rnf = rownames(mat_f)
    
    for(i in 1:(nf-1)) {
        for(j in seq(i+1, nf)) {  
            g1 = cd[rnf[i]]
            g2 = cd[rnf[j]]
            r1 = cd[cd == g1]
            k1 = which(names(r1) == rnf[i])
            r2 = cd[cd == g2]
            k2 = which(names(r2) == rnf[j])
            if (is.na(mat_f[i, j])) {
            } else if (g1 == g2) {
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*LineWidthMultipler, h = (abs(j-i)/(nf-1)*WithinCorrHeightMultiplier_inclass),
                            col = col_fun(mat_f[i, j]), rou1 = 1.01, rou2 = 1.01)
            } else {
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_f[i, j]*LineWidthMultipler, h = (abs(j-i)/(nf-1)*WithinCorrHeightMultiplier_betweenclass + WithinCorrHeightOffset_betweenclass),
                            col = col_fun(mat_f[i, j]))
            }
        }
        
    }
    # circos.clear()
    
    # E. links for within males  -------------------------------------------
    
    mat_m = as.matrix(correlationsample2) 
    
    mat_m[mat_m < corrNALv & mat_m > -corrNALv] <- NA
    
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
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*LineWidthMultipler, h = (abs(j-i)/(nm-1)*WithinCorrHeightMultiplier_inclass),
                            col = col_fun(mat_m[i, j]), rou1 = 1.01, rou2 = 1.01)
            } else {
                circos.link(g1, k1 - 0.5, g2, k2 - 0.5, lwd = mat_m[i, j]*LineWidthMultipler, h = (abs(j-i)/(nm-1)*WithinCorrHeightMultiplier_betweenclass + WithinCorrHeightOffset_betweenclass),
                            col = col_fun(mat_m[i, j]))
            }
        }
    }
    
    # circos.clear()
}



