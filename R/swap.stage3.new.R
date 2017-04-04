swap.stage3.new <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
        
        idx <- seq(1, length(ind) - 1)
        changepoints <- sample(idx, size = 2, replace = FALSE)
        # Two checks for omitting the swap of two indtical animals and animals in the same block
        
        check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
         check2 = diff((changepoints - 1)%/%nPlot) == 0
        check3 = diff(changepoints%%nPlot) == 0
        
        while (check1 || check2 || check3) {
            changepoints <- sample(idx, size = 2, replace = FALSE)
            check1 = all(Z1.mat[ind[changepoints[1]], ] == Z1.mat[ind[changepoints[2]], ])
            check2 = diff((changepoints - 1)%/%nPlot) == 0
            check3 = diff(changepoints%%nPlot) == 0
        }
        
      
    if(Z1.rep < 5){  
      changepoint1 = which(Z1.mat[,which(Z1.mat[ind[changepoints[1]],] == 1)]==1)  
      changepoint2 = which(Z1.mat[,which(Z1.mat[ind[changepoints[2]],] == 1)]==1)  
    }else {
        changepoint1 = changepoints[1]
        changepoint2 = changepoints[2]    
    }
                        
    tmp <- ind[changepoint1]
    ind[changepoint1] <- ind[changepoint2]
    ind[changepoint2] <- tmp
    return(ind)
}
