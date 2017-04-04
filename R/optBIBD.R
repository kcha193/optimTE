optBIBD <-
function(nTrt, bRep, nCag, tRep, nPlot, resDF = NA, confoundCag = FALSE, iter = 10000, speDes = NA) {
     
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
    
    nAni = nTrt * nCag
    
    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:",  nTrt * bRep, "Cag:", nCag, "Tag:", nPlot, "Run:", nBlk, "\n")
    
    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }

    phase1DesignEX1 <- local({
        
        Cag = rep(1:nCag, each = (nTrt * nCag)/nCag)
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Cag, Ani, Trt))
    })
    
    
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)
    
    #browser()
    
    test = matrix(as.numeric(phase1DesignEX1$Trt) - 1, nrow = nCag, byrow = TRUE)
        
    bibd = apply( test, 2,  function(x) multiLetters((x + 0:(nrow(test)-1))%%ncol(test) + 1, FALSE))

    if(all(is.na(speDes))){ 
      newTrt = as.factor(t(bibd[,-(1:(nCag - bRep))]))
    } else{
      newTrt = as.factor(apply(cbind(bibd, speDes), 1, function(x) letters[which(table(x)==1)]))
    }
    
    phase1DesignEX1 <- local({
        
        Cag = rep(1:nCag, each = (nTrt * bRep)/nCag)
        Trt = newTrt
        
        data.frame(cbind(Cag, Trt))
    })
    
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    phase1DesignEX1$Cag = as.factor(phase1DesignEX1$Cag)
    phase1DesignEX1$Ani = as.factor(multiLetters(1:(nrow(phase1DesignEX1))))
 
    nAni = nrow(phase1DesignEX1)
 
  Zc = matrix(0, ncol = nCag, nrow = nAni)  
 Zc[cbind(1:nAni, as.numeric( phase1DesignEX1$Cag))]<-1   
  
  Xt = matrix(0, ncol = nTrt, nrow = nAni)  
  Xt[cbind(1:nAni, as.numeric( phase1DesignEX1$Trt))]<-1   
   
    e.va = Re(eigen( t(Xt) %*% (identityMat(nAni) - projMat(Zc)) %*% Xt, only.values = TRUE)$va)
    
    can.eff = e.va[-which(e.va < 1e-07)]
    
    if(all(e.va == 0) || !all(outer(can.eff, can.eff , "-") <1e-7 )) {
        stop("The Phase 1 design is not a BIBD!")
    }         
    
    n = nAni * tRep
    k = nAni /nCag
    
    aveEff.phase1 = (nTrt * (k-1))/(k * (nTrt-1)) 

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
	
	  # browser() 
 
	  if(confoundCag){
      	Z1.des =  as.numeric(initialBIBD.nolimit(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, 
                        nPlot = nBlk, speDes = speDes))
    } else{
      Z1.des = as.numeric(initialBIBD(nTrt = nTrt, bRep = bRep, nCag = nCag, tRep = tRep, 
                         nPlot = nPlot, speDes = speDes))
     } 
    
    
    
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1
    
    # Parameter's of treatment
    trt.rep = n/nTrt
    
    trt.des = as.numeric(phase1DesignEX1$Trt[match(multiLetters(Z1.des), phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1
    #############################################################################################
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
  			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, 
        resDF = TRUE, upperValue = aveEff.phase1)

     }else{
		
 		firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
			test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage2.new, 
			swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage2.new, 
			Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
			trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, 
      resDF = TRUE, upperValue = aveEff.phase1)
		
    } 
    
    } else{
			firstStep = optThreeStage(init = 1:n, iter = iter/10, obj.fun = obj.fun.RBD, 
			test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
			swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
			Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
			trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
			cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = NA, 
      resDF = TRUE, upperValue = aveEff.phase1)
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
        resDF = FALSE, upperValue = aveEff.phase1)
    }else{
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage2.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage2.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = aveEff.phase1)
    
     }
    } else{
    secondStep = optThreeStage(init = firstStep, iter = iter, obj.fun = obj.fun.RBD1, 
        test.obj.fun = test.obj.fun.RBD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, 
        trt.rep = trt.rep, blk.proj = blk.proj, nTrt = nTrt, C.cage =  C.cage, 
        cage.Rep =  cage.Rep, C.ani = C.ani, ani.Rep = ani.Rep, newResDF = newResDF, 
        resDF = FALSE, upperValue = aveEff.phase1)
    
    }        

    new.Z1.mat = Z1.mat[secondStep,]      
    
    ###############################################################################
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
