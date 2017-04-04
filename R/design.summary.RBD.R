design.summary.RBD <-
function(design.df, simple = TRUE, ...) {

    n = nrow(design.df)
    nBlk = nlevels(design.df$Run)
    nPlot = nlevels(design.df$Tag)
    nCag = nlevels(design.df$Cag)

    nAni = nlevels(design.df$Ani)
    Run.mat = with(design.df, as.matrix(table(1:n, Run)))
    Tag.mat = with(design.df, as.matrix(table(1:n, Tag)))
    Cag.mat = with(design.df, as.matrix(table(1:n, Cag)))
    Ani.mat = with(design.df, as.matrix(table(1:n, paste(Cag, Ani, sep = ""))))
    Trt.mat = with(design.df, as.matrix(table(1:n, Trt)))

    C.cage = (identityMat(nCag) - K(nCag)) %x% K(nAni)
    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))
    cage.Rep = n/nCag

    C.ani = identityMat(nCag) %x% (identityMat(nAni) - K(nAni))

    Pb = projMat(Run.mat)
    Pb1 = projMat(Tag.mat)

    blk.proj = (identityMat(n) - Pb) %*% (identityMat(n) - Pb1)

    info.mat = matMulti(identityMat(n) - Pb, Ani.mat, C.cage)

    if (any(abs(info.mat) > 1e-07)) {
        PP = (identityMat(n) - Pb) - matMulti1((identityMat(n) - Pb), ginv(matMulti((identityMat(n) - Pb), Ani.mat,
            C.cage)), C.cage, Ani.mat)

        PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
    } else {
        PP1 = matMulti1((identityMat(n) - Pb), ginv(matMulti((identityMat(n) - Pb), Ani.mat, C.ani)),
            C.ani, Ani.mat)

    }


    cat("Design parameters:\n")
    cat("Trt:", ncol(Trt.mat), "Ani:", ncol(Ani.mat), "Cag:", ncol(Cag.mat), "Tag:", ncol(Tag.mat), "Run:",
        ncol(Run.mat), "\n")

    if(simple){
      
      print((N = with(design.df, table(Cag, Run))))
  
      print(N %*% t(N))
      }
      cat("Cage and Animal design:\n")
      print(matrix(paste(design.df$Cag, design.df$Ani, sep = ""), nrow = nBlk, ncol = nPlot,
          byrow = TRUE))
  
      cat("Cage and Animal efficiency:\n")
      print(test.RBD(X.trt = identityMat(n), blk.proj = PP1, Rep = 1)[c("nCan", "can.eff")])
  
      info.mat = matMulti(blk.proj, Ani.mat, C.cage)
  
      if (any(abs(info.mat) > 1e-07)) {
          PP = blk.proj - matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.cage)),
              C.cage, Ani.mat)
  
          PP1 = matMulti1(PP, ginv(matMulti(PP, Ani.mat, C.ani)), C.ani, Ani.mat)
      } else {
          PP1 = matMulti1(blk.proj, ginv(matMulti(blk.proj, Ani.mat, C.ani)), C.ani,
              Ani.mat)
  
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
  
      print(test.RBD(X.trt = Trt.mat, PP1, Rep = n/ncol(Trt.mat)))
    

    cat("Phase 1 theoretical ANOVA:\n")
    print(summaryAovOnePhase(design.df, blk.str = "Cag/Ani", trt.str = "Trt", ...))

    cat("Phase 2 theoretical ANOVA:\n")
    print(summaryAovTwoPhase(design.df, blk.str2 = "Run", blk.str1 = "Cag/Ani", trt.str = "Tag + Trt", ...))

}
