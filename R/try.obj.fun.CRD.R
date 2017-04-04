try.obj.fun.CRD <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
        
        Rep = Z1.rep
        newZ1 = Z1.mat[ind[-length(ind)], ]
        info.mat = t(newZ1) %*% blk.proj %*% newZ1
        e.va = eigen(info.mat, only.values = TRUE)$va
        can.eff = e.va[-which(e.va < 1e-07)]/Rep
        
        aveEff = 1/mean(1/can.eff)
        
        #if(aveEff != 1) return(0)
        
        Rep = trt.rep
        newX = X.trt[ind[-length(ind)], ]
        info.mat = t(newX) %*% blk.proj %*% newX
        e.va = eigen(info.mat, only.values = TRUE)$va
        can.eff = e.va[-which(e.va < 1e-07)]/Rep
        len = length(can.eff)
        
        if(len == 0) return(0)
        
        w = (Z1.rep)^2
        aveEff = ((w - 1)/w) * aveEff + 1/w * ((len + 1/mean(1/can.eff))/nTrt)
        
        #aveEff = ((len + 1/mean(1/can.eff))/nTrt)
        
        
        return(aveEff)
    }
