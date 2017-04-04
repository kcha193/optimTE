optThreeStage <- function(init, iter, obj.fun, test.obj.fun, swap.stage1.new, swap.stage2.new, 
    swap.stage3.new, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, 
    cage.Rep, C.ani, ani.Rep, newResDF, resDF, upperValue) {
    
    newInit = c(init, init[1])
    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, 
        C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
        nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
    
    cat("Level: 1, Finding the temperatures for both stages\n")
    
    
    temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, nPlot, 
        Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low
    
    temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun, Z1.mat, X.trt, nPlot, 
        Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
    res.stage1 = temp$newInit
    y2.stage1 = temp$y2
    y2low.stage1 = temp$y2low
    
    temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun, Z1.mat, X.trt, nPlot, 
        Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
    res = temp$newInit
    y2 = temp$y2
    y2low = temp$y2low
    
    cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
    cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
    
    test.newInit = apply(rbind(res, res.stage1, res.stage2), 1, obj.fun, Z1.mat, 
        X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
        newResDF)
    newInit = rbind(res, res.stage1, res.stage2)[which(test.newInit == max(test.newInit))[1], 
        ]
    
    cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, 
        C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
        nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
    
    if (!resDF && isTRUE(all.equal(cur, upperValue))) 
        return(newInit[-length(newInit)])
    
    
    while ((cur - old) > 1e-07) {
        
        temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, 
            nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
            newResDF)
        
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
        
        
        temp = initialTemp(iter, newInit, swap.stage1.new, obj.fun, Z1.mat, X.trt, 
            nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
            newResDF)
        
        res.stage1 = temp$newInit
        y2.stage1 = temp$y2
        y2low.stage1 = temp$y2low
        
        
        temp = initialTemp(iter, newInit, swap.stage3.new, obj.fun, Z1.mat, X.trt, 
            nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
            newResDF)
        
        res = temp$newInit
        y2 = temp$y2
        y2low = temp$y2low
        
        test.newInit = apply(rbind(res, res.stage1, res.stage2), 1, obj.fun, Z1.mat, 
            X.trt, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
            newResDF)
        newInit = rbind(res, res.stage1, res.stage2)[which(test.newInit == max(test.newInit))[1], 
            ]
        
        old = cur
        
        
        cat("Stage 2 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
        cat("Stage 1 temperature range: ", y2.stage1, " to ", y2low.stage1, "\n")
        cat("Stage 0 temperature range: ", y2, " to ", y2low, "\n")
        
        cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, 
            C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
        
        if (!resDF && isTRUE(all.equal(cur, upperValue))) 
            return(newInit[-length(newInit)])
    }
    
    
    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, 
        C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
    cat("Level: 2, Current temp: ", y2.stage2, ",", y2.stage1, ",", y2, "\n")
    
    res <- optim(newInit, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, 
        temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
        Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
        C.ani, ani.Rep, newResDF)
    
    res1 <- optim(res$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, 
        temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
        Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
        C.ani, ani.Rep, newResDF)
    
    res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, 
        temp = y2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, 
        X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, 
        newResDF)
    
    cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, 
        C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
    print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
        nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
    
    if (!resDF && isTRUE(all.equal(cur, upperValue))) 
        return(res1$par[-length(res1$par)])
    
    while ((cur - old) > 1e-07) {
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, 
            temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, 
            temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, 
            temp = y2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        old = cur
        cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        
        print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
        
        
        if (!resDF && isTRUE(all.equal(cur, upperValue))) 
            return(res1$par[-length(res1$par)])
    }
    
    
    inter.stage1 = exp(log(y2.stage1/y2low.stage1)/8)
    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)
    inter = exp(log(y2/y2low)/8)
    
    i = 1
    
    unstable = FALSE
    
    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE
        y2.stage2 = y2.stage2/inter.stage2
        y2.stage1 = y2.stage1/inter.stage1
        y2 = y2/inter
        
        old = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        
        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, ",", y2.stage1, ",", 
            y2, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, 
            temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, 
            temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, 
            temp = y2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
            Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
            C.ani, ani.Rep, newResDF)
        
        cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        
        print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
            nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
        
        if (!resDF && isTRUE(all.equal(cur, upperValue))) 
            return(res1$par[-length(res1$par)])
        
        while ((cur - old) > 1e-07) {
            unstable = TRUE
            
            old = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
                nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
            
            res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, 
                temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
                C.ani, ani.Rep, newResDF)
            
            res1 <- optim(res1$par, obj.fun, swap.stage1.new, method = "SANN", control = list(maxit = iter, 
                temp = y2.stage1, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
                C.ani, ani.Rep, newResDF)
            
            res1 <- optim(res1$par, obj.fun, swap.stage3.new, method = "SANN", control = list(maxit = iter, 
                temp = y2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, 
                C.ani, ani.Rep, newResDF)
            
            
            cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
                nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
            print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, 
                nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
            
            if (!resDF && isTRUE(all.equal(cur, upperValue))) 
                return(res1$par[-length(res1$par)])
            
        }
        
        i = i + 1
        
    }
    
    
    return(res1$par[-length(res1$par)])
} 
