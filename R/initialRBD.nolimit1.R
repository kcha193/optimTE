initialRBD.nolimit1 <-
function(nTrt, bRep, nCag, tRep, nPlot) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt

    nAni = nTrt * bRep

    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }


    phase1DesignEX1 <- local({

        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt

        data.frame(cbind(Cag, Ani, Trt))
    })

   phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)

    nBlk = n/nPlot

    # Parameter's of block structure of Phase 1
    nZ1 = nrow(phase1DesignEX1)
    Z1.rep = n/nZ1

    fillIn = function(run, tag, phase1.mat, count, ani.limit, trt.limit) {

        for (i in 1:length(count)) {

            check1 = sum(phase1.mat[run, , 1] == names(count)[i]) < ani.limit[2]
            check2 = sum(phase1.mat[, tag, 1] == names(count)[i]) < ani.limit[1]
            trt = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == names(count)[i])])
            check3 = sum(phase1.mat[run, , 2] == trt) < trt.limit[2]
            check4 = sum(phase1.mat[, tag, 2] == trt) < trt.limit[1]

             #print(c(check1, check2, check3, check4))

            if (length(count) == 1) {
                return(names(count)[i])
            }

            if (check1 && check2 && check3 && check4) {
                return(names(count)[i])
            } else if (run%%2 == 0 && tag%%2 == 1 && tag == nPlot) {
                return(phase1.mat[run - 1, tag, 1])
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

    ### fill-in then check ####
    ani.char = multiLetters(1:nZ1)
    count = rep(Z1.rep, nZ1)
    names(count) = ani.char

    phase1.mat = array("-1", c(nBlk, nPlot, 2))
    len = (nBlk%/%tRep) * tRep
    ani.limit = c(len, nPlot)/nZ1
    trt.limit = c(Inf, Inf)     
    
    #browser()    
    for (i in 1:len) {
        for (j in 1:nPlot) {
            phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, ani.limit, trt.limit)

            phase1.mat[i, j, 2] = as.character(phase1DesignEX1$Trt[which(phase1DesignEX1$Ani == phase1.mat[i,
                j, 1])])

            count[phase1.mat[i, j, 1]] = count[phase1.mat[i, j, 1]] - 1

            # count = sort(count, d=TRUE)
            if (any(count == 0))
                count = count[-which(count == 0)]
        }
    }
	

    
	if(nBlk == len){
		return(as.numeric(as.factor(phase1.mat[,,1])))
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

                  phase1.mat[i, j, 1] = fillIn(i, j, phase1.mat, count, rep(Z1.rep, 2), c(nBlk/nTrt, nTrt))

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
    
    return(as.numeric(as.factor(phase1.mat[,,1])))
   
}
