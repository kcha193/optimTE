optCRD <-
function(nTrt, bRep, tRep, nPlot, iter = 10000) {
    n = nTrt * bRep * tRep
    nBlk = n/nPlot
    trt.rep = n/nTrt
    
    nAni = nTrt * bRep
    
    if (!is.wholenumber(n/nPlot)) {
        stop("The samples from Phase 1 experiment can not be fitted to the Phase 2 experiment!")
    }
 
    cat("Design parameters:\n")
    cat("Trt:", nTrt, "Ani:", nAni, "Tag:", nPlot, "Run:", nBlk, "\n")
     
    phase1DesignEX1 <- local({
        Ani = 1:nAni
        Trt = 1:nTrt
        
        data.frame(cbind(Ani, Trt))
    })
    phase1DesignEX1$Ani = as.factor(multiLetters(phase1DesignEX1$Ani))
    phase1DesignEX1$Trt = as.factor(multiLetters(phase1DesignEX1$Trt, FALSE))
    
   #browser()
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
    
    
    Z1.des = as.numeric(initialCRD(nTrt = nTrt, bRep = bRep, tRep = tRep, nPlot = nPlot))
    
     
    Z1.mat = matrix(0, ncol = nZ1, nrow = nZ1 * Z1.rep)
    Z1.mat[cbind(1:(nZ1 * Z1.rep), Z1.des)] = 1
    
    
    # Parameter's of treatment
    trt.rep = n/nTrt
    
    trt.des = as.numeric(phase1DesignEX1$Trt[match(multiLetters(Z1.des), phase1DesignEX1$Ani)])
    X.trt = matrix(0, ncol = nTrt, nrow = nTrt * trt.rep)
    X.trt[cbind(1:(nTrt * trt.rep), trt.des)] = 1
   
     
    #Finding the optimal design by SA 
    #if(tRep ==2){
        new.Z1.mat = Z1.mat[optThreeStage(init = 1:n, iter = iter, obj.fun = obj.fun.CRD, 
        test.obj.fun = test.obj.fun.CRD, swap.stage1.new = swap.stage1.new, 
        swap.stage2.new = swap.stage2.new, swap.stage3.new = swap.stage3.new, 
        Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, trt.rep = trt.rep,
        blk.proj = blk.proj, nTrt = nTrt, C.cage = NA, cage.Rep = NA, C.ani = NA, 
        ani.Rep = NA, newResDF = NA, resDF = FALSE, upperValue = 1), 
        ]
    #} else {
    #    new.Z1.mat = Z1.mat[optThreeStage(init = 1:n, iter = iter, obj.fun = obj.fun.CRD, 
     #   test.obj.fun = test.obj.fun.CRD, swap.stage1.new = swap.stage1, 
    #    swap.stage2.new = swap.stage2, swap.stage3.new = swap.stage3, 
    #    Z1.mat = Z1.mat, X.trt = X.trt, nPlot = nPlot, Z1.rep = Z1.rep, trt.rep = trt.rep,
    #    blk.proj = blk.proj, nTrt = nTrt, C.cage = NA, cage.Rep = NA, C.ani = NA, 
    #    ani.Rep = NA, newResDF = NA, resDF = FALSE, upperValue = 1), 
   #     ]
    #}
   

   
    # Set new optimal design     
    colnames(new.Z1.mat) = sort(levels(interaction(multiLetters(1:nZ1))))
    new.Z1.des = apply(new.Z1.mat, 1, function(x) colnames(new.Z1.mat)[which(as.logical(x))])
    new.Z1.des = as.data.frame(t(sapply(strsplit(new.Z1.des, "\\."), function(x) x)))
    
    design.df = data.frame(Run = as.factor(Z.des), Tag = factor(1:nPlot))
    
    design.df = cbind(design.df, Ani = t(new.Z1.des), Trt = phase1DesignEX1[match(as.character(t(new.Z1.des)), 
        as.character(phase1DesignEX1$Ani)), ]$Trt)
     
    return(design.df)
}
