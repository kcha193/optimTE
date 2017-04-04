test.obj.fun.CRD <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
       
        # inforamtion in within runs and tags stratum while the cage information is eliminated
        
        Rep = Z1.rep
        newZ1 = Z1.mat[ind[-length(ind)], ]
        PP = projMat(newZ1) %*% blk.proj
        
        Rep = trt.rep
        newX = X.trt[ind[-length(ind)], ]
        info.mat = t(newX) %*% blk.proj %*% newX
        e.va = eigen(info.mat, only.values = TRUE)$va
        can.eff = e.va[-which(e.va < 1e-07)]/Rep
        len = length(can.eff)
        
        
        return(list(trtAveEff = 100/mean(1/can.eff), resDF = tr(PP) - len))
    }
