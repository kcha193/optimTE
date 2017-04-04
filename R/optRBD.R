optRBD <-
function(nTrt, bRep, nCag, tRep, nPlot, resDF = NA, confoundCag = FALSE, upperValue = 1, iter = 10000, second = FALSE) {
     
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
    
    nAni = nTrt * bRep
    
    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }
     
    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:", nAni, "Cag:", nCag, "Tag:", nPlot, "Run:", nBlk, "\n")
    
    phase1DesignEX1 <- local({
        
        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Cag, Ani, Trt))
    })
    
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)
    
    # Parameter's of block
    nBlk = n/nPlot
    
    Zb = matrix(0, ncol = nBlk, nrow = n)
    Z.des = rep(1:nBlk, each = nPlot)
    Zb[cbind(1:n, Z.des)] <- 1
    
    Zt = matrix(0, ncol = nPlot, nrow = n)
    tag.des = rep(1:nPlot, time = nBlk)
    Zt[cbind(1:n, tag.des)] <- 1
    
    Pb = projMat(Zb)
    Pb1 = projMat(Zt)
    
    betRun = Pb - K(n)
    betTag = Pb1 - K(n)
    withBlock = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)
    
    blk.proj = withBlock
    
    # Parameter's of block structure of Phase 1
    nZ1 = nrow(phase1DesignEX1)
    Z1.rep = n/nZ1
	

  #browser() 
	if(confoundCag){
		Z1.des = as.numeric(initialRBD.nolimit1(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, 
                        nPlot = nBlk))
 
  
  } else{
	  	Z1.des = as.numeric(initialRBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep,
                        nPlot = nPlot))
	}
 
 
 
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1
    
    # Parameter's of treatment
    trt.rep = n/nTrt
    
    trt.des = as.numeric(phase1DesignEX1$Trt[match(multiLetters(Z1.des), phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1

    #Contrast matrix for cage and animals 
    C.cage  = (identityMat(nCag) - K(nCag)) %x%  K(nAni/nCag) 
    cage.Rep = n/nrow(C.cage)

    C.ani =   identityMat(nCag)%x%  (identityMat(nAni/nCag) - K(nAni/nCag))
    ani.Rep = n/nAni
	
	if(is.na(resDF)){
		print("Step 1 optimisation:")
		
	 if(confoundCag){
    if((nCag == 2 || nCag == 6 || nCag == 10) && nPlot == 8){
  		firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
  			test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new1, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new1, 
  			Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
  			trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
  			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, resDF = TRUE, upperValue = upperValue)

     }else{
		
 		firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
			test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage2.new, 
			swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage2.new, 
			Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
			trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, resDF = TRUE, upperValue = upperValue)
		
    } 
    
    } else{
			firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
			test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
			swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
			Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
			trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, resDF = TRUE, upperValue = upperValue)
		}

		
			
		temp1 = as.numeric(obj.fun.RBD(c(firstStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
								  trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
								  ani.Rep, newResDF)) 
		temp2 = as.numeric(obj.fun.RBD(c(1:n,1), Z1.mat, X.trt, nPlot, Z1.rep, 
								  trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
								  ani.Rep, newResDF))
		
		if((temp1 - temp2) < 1e-7){
			temp1 = as.numeric(obj.fun.RBD1(c(firstStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
								  trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
								  ani.Rep, newResDF)) 
			temp2 = as.numeric(obj.fun.RBD1(c(1:n,1), Z1.mat, X.trt, nPlot, Z1.rep, 
								  trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
								  ani.Rep, newResDF))		
		}	
		
		firstStep = rbind(firstStep, 1:n)[which(c(temp1,temp2)==max(c(temp1, temp2)))[1],]	
		
		
		newResDF = test.obj.fun.RBD(c(firstStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
								  trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
								  ani.Rep, newResDF)$res   
    } else{
		newResDF = resDF 
		firstStep = 1:n
	}
	
    print("Step 2 optimisation:")
    
  if(confoundCag){
    if((nCag == 2 || nCag == 6 || nCag == 10) && nPlot == 8){
    
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new1, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new1, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
    }else{
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage2.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage2.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
    
     }
    } else{
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
    
    }
    
  if(second){    
    ans <- readline("Further improve the design (Y/N/C)?: ")    
    
	while((ans != "Y") && (ans != "N")&& (ans != "C")) ans <-  readline("Further improve the design (Y/N/C)?: ")    
	
    while((ans == "Y") || ans == "C"){ 
		
		if(ans == "Y") {
			reduce <-  readline("Reduce by?: ") 
			
			newResDF = newResDF - as.numeric(reduce)
		
		  if(confoundCag){
    #if(confoundCag && nPlot == 4){
    secondStep.new = optThreeStage(init = secondStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage2.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage2.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
    } else{
    secondStep.new = optThreeStage(init = secondStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
    
    }
       
        } else {
			secondStep.new = optThreeStage(init = secondStep, iter = iter, obj.fun = obj.fun.RBD1, 
				  test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
				  swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
				  Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
				  trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
				  cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = upperValue)
 		}		
		
        ans <- readline("Further improve the design (Y/N/C/K)?: ")    
 		
		while((ans != "Y") && (ans != "N")&& (ans != "C") && (ans != "K")) ans <-  readline("Further improve the design (Y/N/C/K)?: ")      

        if(ans == "N") { 
		
             temp1 =  as.numeric(test.obj.fun.RBD(c(secondStep.new,1), Z1.mat, X.trt, nPlot, Z1.rep, 
                              trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
                              ani.Rep, newResDF)$aveEff)
             
             temp2 =  as.numeric(test.obj.fun.RBD(c(secondStep,1), Z1.mat, X.trt, nPlot, Z1.rep, 
                              trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, 
                              ani.Rep, newResDF)$aveEff)                
             
            if(isTRUE(all.equal(temp1, temp2))){
				secondStep  = secondStep   
            }else{
				secondStep  = rbind(secondStep, secondStep.new)[which(c(temp2, temp1) == max(c(temp2, temp1)))[1],]
				
            } 
			break
        } else if(ans == "K") {
			secondStep = secondStep.new
			
			break
		}
       
		
        secondStep = secondStep.new
		
    } 
    }
    new.Z1.mat = Z1.mat[secondStep,]      
    
    # Set the design from the search
    
    colnames(new.Z1.mat) = sort(levels(interaction(multiLetters(1:nZ1))))
    new.Z1.des = apply(new.Z1.mat, 1, function(x) colnames(new.Z1.mat)[which(as.logical(x))])
    new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))
    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
    design.df = cbind(design.df, 
              Cag = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
                                    as.character(phase1DesignEX1$Ani)), ]$Cag, 
              Ani = t(new.Z1.des), 
              Trt = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
                                    as.character(phase1DesignEX1$Ani)), ]$Trt)
    levels(design.df$Ani) = multiLetters(rep(1:(nlevels(design.df$Ani)/nCag), nCag) )

    return(design.df)
}
