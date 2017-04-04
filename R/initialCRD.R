initialCRD <- function(nTrt, bRep, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
     
    nAni = nTrt * bRep
    
    
    phase1DesignEX1 <- local({
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Ani, Trt))
    })
    
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    
    
    if (tRep %% 2 == 0) {
        nZ1 = nrow(phase1DesignEX1)
        Z1.rep = n/nZ1
        
        fillIn = function(run, tag, phase1.mat, count, ani.limit, trt.limit) {
            
            for (i in 1:length(count)) {
                
                check1 = sum(phase1.mat[run, , 1] == names(count)[i]) < ani.limit[2]
                check2 = sum(phase1.mat[, tag, 1] == names(count)[i]) < ani.limit[1]
                trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == names(count)[i])])
                check3 = sum(phase1.mat[run, , 2] == trt) < trt.limit[2]
                check4 = sum(phase1.mat[, tag, 2] == trt) < trt.limit[1]
                
                # print(c(check1, check2, check3, check4))
                
                if (length(count) == 1) {
                  return(names(count)[i])
                }
                
                if (check1 && check2 && check3 && check4) {
                  return(names(count)[i])
                } else if (run%%2 == 0 && tag%%2 == 1) {
                  return(phase1.mat[run - 1, tag + 1, 1])
                } else if (run%%2 == 0 && tag%%2 == 0) {
                  return(phase1.mat[run - 1, tag - 1, 1])
                } else if (run == nBlk) {
                  return(names(count)[i])
                }
                
            }
            
            return(names(which(count == max(count)))[1])
        }
        
        #browser()
        ### fill-in then check ####
        ani.char = multiLetters(1:nZ1)
        count = rep(Z1.rep, nZ1)
        names(count) = ani.char
        
        phase1.mat = array("-1", c(nBlk, nPlot, 2))
        
        len = (nBlk%/%tRep) * tRep
        
        ani.limit = c(len, nPlot)/nZ1
        trt.limit = c(len, nPlot)/nTrt
        
        if(len != 0){
        for (i in 1:len) {
            for (j in 1:nPlot) {
                
                phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, trt.limit)
                
                phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i, j, 
                  1])])
                
                count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
                
                # count = sort(count, d=TRUE)
                if (any(count == 0)) 
                  count = count[-which(count == 0)]
            }
        }
        }
        
        if (nBlk != len) {
            if (nBlk%%tRep == 1) {
                for (j in 1:nPlot) {
                  phase1.mat[nBlk, j, 1] = fillIn(nBlk, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, trt.rep))
                  
                  phase1.mat[nBlk, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[nBlk, 
                    j, 1])])
                  count[phase1.mat[nBlk, j, 1]] = count[phase1.mat[nBlk, j, 1]] - 1
                  
                  # count = sort(count, d=TRUE)
                  if (any(count == 0)) 
                    count = count[-which(count == 0)]
                }
                
            } else if (nBlk%%tRep == 2) {
                for (i in (len + 1):nBlk) {
                  for (j in 1:nPlot) {
                    
                    phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep/4, 2), c(nBlk/nTrt, nTrt))
                    
                    phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i, 
                      j, 1])])
                    
                    count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
                    
                    # count = sort(count, d=TRUE)
                    if (any(count == 0)) 
                      count = count[-which(count == 0)]
                  }
                }
            } else {
                for (i in (len + 1):nBlk) {
                  for (j in 1:nPlot) {
                    
                    ani.limit = c(length((len + 1):nBlk), nPlot)/length(count)
                    phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, c(nBlk/nTrt, trt.rep/tRep))
                    
                    phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i, 
                      j, 1])])
                    
                    count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1
                    
                    # count = sort(count, d=TRUE)
                    if (any(count == 0)) 
                      count = count[-which(count == 0)]
                  }
                }
                
            }
        }
        # print(phase1.mat)
        
        as.numeric(as.factor(t(phase1.mat[, , 1])))
        
    } else if (tRep == 3) {
        temp = 0:(nPlot - 1)
        reverse1 = rep(c(2, 3, 1, 4), nPlot/4) + rep(seq(from = 0, length = nPlot/4, by = 4), each = 4)
        reverse2 = rep(c(3, 1, 2, 4), nPlot/4) + rep(seq(from = 0, length = nPlot/4, by = 4), each = 4)
        
        
        return(rep(c(temp, temp[reverse1], temp[reverse2]), nBlk/tRep) + rep(0:(nBlk/tRep - 1), each = nPlot * tRep) * 
            nPlot + 1)
        
        
    } else {
        return(rep(1:nAni, each = tRep))
        
    }
    
} 

