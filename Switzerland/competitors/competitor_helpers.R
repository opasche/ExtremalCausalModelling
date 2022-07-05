library(pcalg)


pc_algorithm <- function(dat, rank_version=FALSE, alpha=0.01, labels=NULL, u2pd="retry",
                         conservative=FALSE, maj.rule=FALSE, skel.method="stable", ...){#0.01 5e-2
  
  n <- NROW(dat)
  p <- NCOL(dat)
  if(is.null(labels)){
    labels <- paste0("X",1:p)
  }
  
  if(rank_version){
    suff_stat <- list(C = 2 * sin(stats::cor(dat, method="spearman") * pi/6), n=n)
  }else{
    suff_stat <- list(C=stats::cor(dat), n=n)
  }
  
  output <- tryCatch({
    pc_fit <- pcalg::pc(suffStat=suff_stat, indepTest=pcalg::gaussCItest,
                        labels=labels, alpha=alpha, u2pd=u2pd,#"retry" "relaxed"
                        skel.method=skel.method, conservative=conservative, maj.rule=maj.rule, ...)#"stable" "stable.fast"
    cpdag <- unname(methods::as(pc_fit@graph, "matrix"))
    effects <- which(cpdag==1, arr.ind=TRUE)
    colnames(effects) <- c("cause", "effect")
    list(cpdag=cpdag, effects=effects, fit=pc_fit)
  },
  error=function(e){
    warning(paste0("Error in 'pc_algorithm' call: ", e, "."))
    return(list(cpdag=NA, effects=NA, fit=NA))
  })
  
  return(output)
}

LiNGAM <- function(dat, labels=NULL, verbose=FALSE){
  
  n <- NROW(dat)
  p <- NCOL(dat)
  if(is.null(labels)){
    labels <- paste0("X",1:p)
  }
  
  output <- tryCatch({
    lingam_fit <- pcalg::lingam(dat, verbose=verbose)
    #lingam_fit$coeffs <- lingam_fit$Bpruned
    lingam_fit$adj <- t(lingam_fit$Bpruned != 0)
    # colnames(lingam_fit$adj) <- labels
    # rownames(lingam_fit$adj) <- labels
    lingam_fit$effects <- which(lingam_fit$adj, arr.ind=TRUE)
    colnames(lingam_fit$effects) <- c("cause", "effect")
    lingam_fit
  },
  error=function(e){
    warning(paste0("Error in 'LiNGAM' call: ", e, "."))
    return(list(Bpruned=NA, stde=NA, ci=NA, adj=NA, effects=NA))
  })
  output$labels <- labels
  return(output)
}

causal_cpdag_helper <- function(cpdag, ind_cause, ind_effect){
  if(cpdag[ind_cause, ind_effect]==0 & cpdag[ind_effect, ind_cause]==0){
    return("cond. independent")
  }
  if(cpdag[ind_cause, ind_effect]==1 & cpdag[ind_effect, ind_cause]==1){
    return("related")
  }
  if(cpdag[ind_cause, ind_effect]==1 & cpdag[ind_effect, ind_cause]==0){
    return("causal")
  }
  if(cpdag[ind_cause, ind_effect]==0 & cpdag[ind_effect, ind_cause]==1){
    return("inverse causal")
  }
}

# dat <- matrix(rnorm(3000),ncol=3)
# X1=rt(1000,4)
# X2=0.5*X1+rt(1000,4)
# dat <- cbind(X1,X2,0.5*X1+0.5*X2+rt(1000,4))
# 
# pcres <- pc_algorithm(dat, rank_version=TRUE, alpha=0.01, labels=NULL)
# pcres$fit
# pcres$cpdag
# pcres$effects
# 
# linres <- LiNGAM(dat)
# linres$adj
# linres$effects
