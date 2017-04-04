initialTemp <-
function(iter, newInit, swap, obj.fun, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, 
            blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF){
  
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  U = numeric(iter)
  temp = matrix(0, nrow = iter, ncol = length(newInit))
  temp[1,] = newInit
  U[1] =   obj.fun(temp[1,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
  
  tol = 10
  check.counter = 0
  check.U = numeric(tol) 
  #browser()
  for(i in 2:iter){
    setTxtProgressBar(pb, i)
    temp[i,] = swap(temp[i-1,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
    U[i] = obj.fun(temp[i,], Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF)
	check.U = c(check.U, U[i])
    check.U = check.U[-1]
	#print(check.U)
	check.counter = check.counter + 1
    if((check.counter > tol) && all(outer(check.U,check.U, "-") < 1e-7)) {
	  temp[i,] = newInit
      check.U = numeric(tol)
	  check.counter = 0
    }		
  }
  close(pb)
  
  if(all(U==0) || all(diff(U)<1e-8)){
    return(list(y2 = 1, y2low = 1e-7, newInit = newInit))
  
  }else{
    xxx = max(U, na.rm = TRUE) - range(U[which(U < (max(U, na.rm = TRUE)-1e-10))], na.rm = TRUE)
    
    if(isTRUE(all.equal(xxx[1], xxx[2]))) xxx[2] = 1e-7
     
    newInit = temp[which(U == max(U, na.rm = TRUE))[1],]
    return(list(y2 = max(xxx), y2low = min(xxx), newInit =newInit))
  }
}
