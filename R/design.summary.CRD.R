design.summary.CRD <-
function(design.df, simple = TRUE, ...) {

    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, Ani)))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (identityMat(n) - Pb)

    cat("Design parameters:\n")
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Tag:", ncol(Tag.mat), "Run:",
        ncol(Run.mat), "\n")

   
      cat("Animal design:\n")
      print(matrix(design.df$Ani, nrow = nBlk, ncol = nPlot, byrow = TRUE))
     
     if(simple){
      cat("Animal incidence matrix:\n")
      print((N = with(design.df, table(Ani, Run))))
  
      cat("Animal concurrence matrix:\n")
      print(N %*% t(N))
  
      cat("Animal efficiency:\n")
      print(test.CRD(X.trt = Ani.mat, blk.proj, Rep = n/ncol(Ani.mat)))
      }
      
      cat("Treatment design:\n")
      print(matrix(design.df$Trt, nrow = nBlk, ncol = nPlot, byrow = TRUE))
     if(simple){
      cat("Treatment incidence matrix:\n")
      print((N = with(design.df, table(Trt, Run))))
  
      cat("Treatment concurrence matrix:\n")
      print(N %*% t(N))
      }
      cat("Treatment efficiency:\n")
      print(test.CRD(X.trt = Trt.mat, (identityMat(n) - Pb) %*% (identityMat(n) - Pb1), Rep = n/ncol(Trt.mat)))
    

    cat("Phase 1 theoretical ANOVA:\n")
    print(summaryAovOnePhase(design.df, blk.str = "Ani", trt.str = "Trt", ...))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Ani", trt.str = "Tag + Trt", ...))


}
