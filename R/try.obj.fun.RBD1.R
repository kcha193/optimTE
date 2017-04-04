try.obj.fun.RBD1 <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
       
    #inforamtion in within runs and tags stratum while the cage information is eliminated
    info.mat = matMulti(blk.proj ,  Z1.mat[ind[-length(ind)],], C.cage)
    PP = matMulti1(blk.proj, ginv(info.mat), C.cage,  Z1.mat[ind[-length(ind)],])
     
    if(any(abs(info.mat) > 1e-7)){     #check for any information in the cage stratum
     # e.va = Re(eigen(info.mat, only.values = TRUE)$va)
     # can.eff = e.va[-which(e.va < 1e-07)]/cage.Rep
     # aveEff1 = 1/mean(1/can.eff)
     # if(!isTRUE(all.equal(aveEff1, 1))) return(0)

      PP = blk.proj - PP 
    } else{       
    #if there is no information in the cage stratum then no decomposition is required. 
    #If decomposition is performed here, 
    # some very werid thing is going to happen. 
    #This is because taking the inverse of a very small number will become huge. 
    # Hence, large information of the within block within runs can be extracted. 
      PP = blk.proj
    }
    
    #Maximise the amount of animal information in within runs and tags
    info.mat = matMulti(PP , Z1.mat[ind[-length(ind)],], C.ani)
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    can.eff = e.va[-which(e.va < 1e-07)]/ani.Rep 

    aveEff1 = 1/mean(1/can.eff)
     
    if(is.na(aveEff1)) aveEff1 = 0
   
    if(!isTRUE(all.equal(aveEff1, 1))) return(0)
     

    PP =  matMulti1(PP, ginv(info.mat), C.ani, Z1.mat[ind[-length(ind)],]) 
            
    newX = X.trt[ind[-length(ind)], ]
    info.mat = t(newX) %*%  PP %*% newX
    e.va = Re(eigen(info.mat, only.values = TRUE)$va)
    #e.va = Re(svd(info.mat)$d)
    can.eff = e.va[-which(e.va < 1e-07)]/trt.rep
    len = length(can.eff)

    if(!isTRUE(all.equal(len/(nTrt-1), 1))) return(0)
    if(!isTRUE(all.equal((tr(PP) - len), newResDF))) return(0)
    #if((tr(PP) - len)< newResDF) return(0)  
   
    aveEff = 1/(mean(1/can.eff))
     
    return(as.numeric(aveEff))
}
