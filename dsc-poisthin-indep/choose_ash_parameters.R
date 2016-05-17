choose.ash.parameters <- function(input, args) {
  # Assuming that sebetahat/scale are exchangeable
  # For voom transformation, scale is not constant
  scale = rep(1, length(input$sebetahat.voom))
  if (args$transform=="voom"){
    lim = lmFit(input$v)
    scale = lim$stdev.unscaled[,2]
  }else if (args$transform=="RUVvoom"){
    lim = lmFit(input$RUVv)
    scale = lim$stdev.unscaled[,2]
  }else if (args$transform=="SVAvoom"){
    lim = lmFit(input$SVAv)
    scale = lim$stdev.unscaled[,2]
  }
  
  if (args$transform == "voom") {
    sebetahat <- input$sebetahat.voom
    betahat <- input$betahat.voom
    df <- input$df.voom
  }else if (args$transform == "RUVvoom"){
    sebetahat <- input$sebetahat.RUVvoom
    betahat <- input$betahat.RUVvoom
    df <- input$df.RUVvoom
  }else if (args$transform == "SVAvoom"){
    sebetahat <- input$sebetahat.SVAvoom
    betahat <- input$betahat.SVAvoom
    df <- input$df.SVAvoom
  }else if (args$transform == "quasibinom") {
    betahat <- input$betahat.qb
    sebetahat <- input$sebetahat.qb
    df <- input$df.qb
  } else if (args$transform == "Myrna+quasibinom") {
    betahat <- input$betahat.Myrnaqb
    sebetahat <- input$sebetahat.Myrnaqb
    df <- input$df.Myrnaqb
  } else if (args$transform == "Myrnaoff+quasibinom") {
    betahat <- input$betahat.Myrnaoffqb
    sebetahat <- input$sebetahat.Myrnaoffqb
    df <- input$df.Myrnaoffqb
  } else if (args$transform == "RUV+quasibinom") {
    betahat <- input$betahat.RUVqb
    sebetahat <- input$sebetahat.RUVqb
    df <- input$df.RUVqb
  } else if (args$transform == "SVA+quasibinom") {
    betahat <- input$betahat.SVAqb
    sebetahat <- input$sebetahat.SVAqb
    df <- input$df.SVAqb
  } else if (args$transform == "edgeRglm") {
    betahat <- input$betahat.edgeRglm
    sebetahat <- input$sebetahat.edgeRglm
    df <- input$df.edgeRglm
  } else if (args$transform == "DESeqglm") {
    betahat <- input$betahat.DESeqglm
    sebetahat <- input$sebetahat.DESeqglm
    df <- input$df.DESeqglm
  } else if (args$transform == "DESeq2glm") {
    betahat <- input$betahat.DESeq2glm
    sebetahat <- input$sebetahat.DESeq2glm
    df <- input$df.DESeq2glm
  } else if (args$transform == "cate") {
    betahat <- input$betahat.cate
    sebetahat <- input$sebetahat.cate
    df <- input$df.cate
  }
  
  return(list(betahat = betahat, sebetahat = sebetahat, df = df, scale=scale))
} 
