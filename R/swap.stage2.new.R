swap.stage2.new <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
        
    idx <- seq(1, length(ind) - 1)
    changepoint1 <- sample(idx, size = 1)
   
    s = unique(c(seq(changepoint1, max(idx), nPlot), seq(changepoint1, 1, -nPlot)))

    if(length(s[-which(s == changepoint1)]) == 1){
      changepoint2 = s[-which(s == changepoint1)]
    } else {
      changepoint2 = sample(s[-which(s == changepoint1)], size=1)
    }
      
    if(Z1.rep < 5){  
      changepoint1 = which(Z1.mat[,which(Z1.mat[ind[changepoint1],] == 1)]==1)  
      changepoint2 = which(Z1.mat[,which(Z1.mat[ind[changepoint2],] == 1)]==1)  
    }
                        
    tmp <- ind[changepoint1]
    ind[changepoint1] <- ind[changepoint2]
    ind[changepoint2] <- tmp
    return(ind)
}
