optOneStage <-
function(init, iter, obj.fun, test.obj.fun, swap.stage1.new, swap.stage2.new, 
    swap.stage3.new, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, 
    cage.Rep, C.ani, ani.Rep, newResDF, resDF, upperValue) {
    
   newInit = c(init, init[1])
    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))

    cat("Level: 1, Finding the temperatures for stages\n")

    temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

    res.stage2 = temp$newInit
    y2.stage2 = temp$y2
    y2low.stage2 = temp$y2low
    
    cat("Stage 0 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
 
    newInit =  temp$newInit
    
    cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
      print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
	  
    if(!resDF && isTRUE(all.equal(cur, upperValue))) return(newInit[-length(newInit)])

    
    while ((cur - old) > 1e-07) {
        
        temp = initialTemp(iter, newInit, swap.stage2.new, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    
        res.stage2 = temp$newInit
        y2.stage2 = temp$y2
        y2low.stage2 = temp$y2low
    
               
       newInit =  temp$newInit
        
        old = cur
    
        cat("Stage 0 temperature range: ", y2.stage2, " to ", y2low.stage2, "\n")
    
        cur = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        print(test.obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))  

        if(!resDF && isTRUE(all.equal(cur, upperValue))) return(newInit[-length(newInit)])
     }
      

    old = obj.fun(newInit, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
  
    cat("Level: 2, Current temp: ", y2.stage2, "\n")

    res <- optim(newInit, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
        tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

    res1 <- res
  
    cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

    print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))
    
    if(!resDF && isTRUE(all.equal(cur, upperValue))) return(res1$par[-length(res1$par)])

    while ((cur - old) > 1e-07) {
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

        old = cur
        cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

         print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))


        if(!resDF && isTRUE(all.equal(cur, upperValue))) return(res1$par[-length(res1$par)])
    }


    inter.stage2 = exp(log(y2.stage2/y2low.stage2)/8)

    i = 1
    
    unstable = FALSE

    # for (i in 1:8) {
    while (i < 9 || unstable) {
        unstable = FALSE
        y2.stage2 = y2.stage2/inter.stage2
 		
		old = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
		
        cat("Level: ", i + 2, ", Current temp: ", y2.stage2, "\n")
        res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter, temp = y2.stage2,
            tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
				
        cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
        
		 print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))

        if(!resDF && isTRUE(all.equal(cur, upperValue))) return(res1$par[-length(res1$par)])
        
        while ((cur - old) > 1e-07) {
            unstable = TRUE
			
			old = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)

            res1 <- optim(res1$par, obj.fun, swap.stage2.new, method = "SANN", control = list(maxit = iter,
                temp = y2.stage2, tmax = iter/10, trace = TRUE, REPORT = 10, fnscale = -1), 
                Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
			  
         
            cur = obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
                print(test.obj.fun(res1$par, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF))

            if(!resDF && isTRUE(all.equal(cur, upperValue))) return(res1$par[-length(res1$par)])

        }

        i = i + 1
          
    }

 
    return(res1$par[-length(res1$par)])
}
