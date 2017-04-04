obj.fun.CRD <-
function(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF) {
        
        ans = try(try.obj.fun.CRD(ind, Z1.mat, X.trt, nPlot, Z1.rep, trt.rep, blk.proj, nTrt, C.cage, cage.Rep, C.ani, ani.Rep, newResDF), silent = TRUE)
        
        ifelse(class(ans) == "try-error", 0, ans)
    }
