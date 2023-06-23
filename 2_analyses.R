library(hdp)
library(igraph)
################################
# EXTRACT CONVERGED COMPONENTS #
################################
extractConvergedComponentsNOCI <- function(chain) {
  if (is.null(chain)) return(NULL)
  n_features <- nrow(chain@clust_categ_counts[[1]])
  actual_numclust <- unlist(Map(ncol, clust_categ_counts(chain)))
  convChains_cond <- which(actual_numclust==names(which.max(table(actual_numclust))))
  ccc_conv <- clust_categ_counts(chain)[convChains_cond]
  ccc_conv <- lapply(ccc_conv, function(x) {
    x <- t(x)
    x <- x[rowSums(x)!=0,,drop=FALSE]
    x <- x/rowSums(x)
    return(x)
  })
  ccc_ncomp <- unique(unlist(Map(nrow,ccc_conv)))
  stopifnot(length(ccc_ncomp)==1)
  ccc_convComp <- list()
  for (compIdx in 1:ccc_ncomp) {
    ccc_convComp[[compIdx]] <- sapply(ccc_conv, function(x) {
      return(x[compIdx,])
    })
  }
  FinalComps4Chain <- matrix(0,nrow=length(ccc_convComp),ncol=n_features)
  for (compIdx in 1:ccc_ncomp) {
    for (j in 1:n_features) {
      presence_across_samplings <- ccc_convComp[[compIdx]][j,]
      FinalComps4Chain[compIdx,j]<- mean(presence_across_samplings)
    }
  }
  return(FinalComps4Chain)
}

extractConvergedComponents <- function(chain) {
  if (is.null(chain)) return(NULL)
  n_features <- nrow(chain@clust_categ_counts[[1]])
  actual_numclust <- unlist(Map(ncol, clust_categ_counts(chain)))
  convChains_cond <- which(actual_numclust==names(which.max(table(actual_numclust))))
  ccc_conv <- clust_categ_counts(chain)[convChains_cond]
  ccc_conv <- lapply(ccc_conv, function(x) {
    x <- t(x)
    x <- x[rowSums(x)!=0,,drop=FALSE]
    x <- x/rowSums(x)
    return(x)
  })
  ccc_ncomp <- unique(unlist(Map(nrow,ccc_conv)))
  stopifnot(length(ccc_ncomp)==1)
  ccc_convComp <- list()
  for (compIdx in 1:ccc_ncomp) {
    ccc_convComp[[compIdx]] <- sapply(ccc_conv, function(x) {
      return(x[compIdx,])
    })
  }
  FinalComps4Chain <- matrix(0,nrow=length(ccc_convComp),ncol=n_features)
  for (compIdx in 1:ccc_ncomp) {
    for (j in 1:n_features) {
      presence_across_samplings <- ccc_convComp[[compIdx]][j,]
      FinalComps4Chain[compIdx,j]<- ifelse(coda::HPDinterval(coda::as.mcmc(presence_across_samplings))[1]>0, mean(presence_across_samplings), 0)
    }
  }
  which2keep <- which(rowSums(FinalComps4Chain)!=0)
  FinalComps4Chain <- FinalComps4Chain[rowSums(FinalComps4Chain)!=0,,drop=FALSE]
  return(list("comps"=FinalComps4Chain,"which2keep"=which2keep))
}

#######################################
# FUNCTIONS TO COMPUTE SOFT DISTANCES #
#######################################
AssociationsCos2 <- function(a, b, k1, k2) {
  a <- a/rowSums(a)
  b <- b/rowSums(b)
  cc <- a %*% t(b)
  dd <- sqrt(diag(a %*% t(a))) %o% sqrt(diag(b %*% t(b)))
  cosdist <- cc/dd
  allzero_x <- apply(cosdist, 1, function(x) all(x==0))
  allzero_y <- apply(cosdist, 2, function(x) all(x==0))
  # stopifnot(all(!allzero_x))
  # stopifnot(all(!allzero_y))
  # one <- cbind(paste0(1:nrow(cosdist), "_", k1), paste0(apply(cosdist, 1, which.max), "_", k2))
  # one <- one[!allzero_x,]
  # two <- cbind(paste0(apply(cosdist, 2, which.max), "_", k1), paste0(1:ncol(cosdist), "_", k2))
  # two <- two[!allzero_y,]
  # jt <- rbind(one,two)
  # jt_dup <- jt[duplicated(jt),]
  # jt  <- unique(jt)
  one <- cbind(paste0(1:nrow(cosdist), "_", k1), paste0(apply(cosdist, 1, which.max), "_", k2))
  one <- one[!allzero_x,,drop=FALSE]
  two <- cbind(paste0(1:ncol(cosdist), "_", k2), paste0(apply(cosdist, 2, which.max), "_", k1))
  two <- two[!allzero_y,,drop=FALSE]
  gd <- graph_from_data_frame(rbind(one,two))
  multi_in_nodes <- V(gd)$name[degree(gd, mode="in")>1]
  gu <- graph_from_data_frame(rbind(one,two), directed=FALSE)
  gu_dup <- get.edgelist(gd)[duplicated(get.edgelist(gu)),]
  ml_nodes <- unique(as.vector(gu_dup)[as.vector(gu_dup) %in% multi_in_nodes])
  gu_next <- delete_vertices(gu, ml_nodes)
  leftover_edges <- get.edgelist(gu_next)
  leftover_nodes <- V(gu_next)$name[degree(gu_next)==0]
  leftover_edges_extra <- get.edgelist(gd)[(get.edgelist(gd)[,1] %in% leftover_nodes) | ((get.edgelist(gd)[,2] %in% leftover_nodes)),]
  pairwise_graph <- graph_from_data_frame(rbind(gu_dup, leftover_edges, leftover_edges_extra), vertices=c(paste0(1:nrow(cosdist), "_", k1), paste0(1:ncol(cosdist), "_", k2)), directed=FALSE)
  pairwise_graph <- simplify(pairwise_graph)
  return(pairwise_graph)
}

AssociationsCos <- function(a, b, k1, k2) {
  cc <- a %*% t(b)
  dd <- sqrt(diag(a %*% t(a))) %o% sqrt(diag(b %*% t(b)))
  cosdist   <- cc / dd
  allzero_x <- apply(cosdist, 1, function(x) all(x==0))
  # allzero_y <- apply(cosdist, 2, function(x) all(x==0))
  # stopifnot(all(!allzero_x))
  # stopifnot(all(!allzero_y))
  one   <- cbind(paste0(1:nrow(cosdist), "_", k1), paste0(apply(cosdist, 1, which.max), "_", k2))
  one <- one[!allzero_x,,drop=FALSE]
  onejt <- graph_from_data_frame(one, vertices=c(paste0(1:nrow(cosdist), "_", k1), paste0 (1:ncol(cosdist), "_", k2)), directed=FALSE)
  onejt <- simplify(onejt)
  return(onejt)
}



#########################
# SOFT DISTANCE MEASURE #
#########################
soft_dist <- function(compsObjList) {
  cosgraph <- make_empty_graph(directed=FALSE)
  for (k1 in 1:length(compsObjList)) {
    # message(k1)
    if(k1==length(compsObjList)){break}
    for (k2 in (k1+1):length(compsObjList)) {
      tmpAssocs <- AssociationsCos2(a=compsObjList[[k1]],
                                    b=compsObjList[[k2]],
                                    k1=k1,
                                    k2=k2)
      cosgraph <- cosgraph + tmpAssocs
    }
  }
  cosgraph <- simplify(cosgraph)
  return(cosgraph)
}

rm_toNA <- function(cosgraph, K) {
  names_split <- as.data.frame(stringr::str_split_fixed(V(cosgraph)$name,"_",3))
  whichNA_nodes <- which(names_split[,1]=="NA")
  chs2rm <- as.vector(unlist(apply(names_split[whichNA_nodes,c(2,3)], 1, function(x) (as.integer(x[1]) + K*(as.integer(x[2])-1)))))
  return(chs2rm)
}

rm_unconverged <- function(cosgraph, crit_conv) {
  unconv_chains <- c()
  K <- length(crit_conv)
  for (kn in 1:K) {
    if (kn %in% which(!crit_conv)) {
      unconv_chains <- c(unconv_chains, paste(kn,sep="_"))
    }
  }
  comps_to_remove <- V(cosgraph)$name[stringr::str_split_fixed(V(cosgraph)$name,"_",2)[,2] %in% unconv_chains]
  cosgraph <- delete.vertices(cosgraph, comps_to_remove)
  return(cosgraph)
}

rm_iso <- function(conv_g) {
  while (TRUE) {
    iso_vertices <- which(degree(conv_g)==0)
    if (length(iso_vertices)<1){
      break
    }
    chains_not_explicitly_converged <- unique(stringr::str_split_fixed(names(iso_vertices),"_",2)[,2])
    comps_to_remove <- which(stringr::str_split_fixed(V(conv_g)$name,"_",2)[,2] %in% chains_not_explicitly_converged)
    conv_g <- delete.vertices(conv_g, comps_to_remove)
  }
  return(conv_g)
}


extract_overall_components <- function(cosgraph, nchains) {
  cc <- components(cosgraph)
  collectRobustComps <- make_empty_graph(directed=FALSE)
  for (membIdx in 1:cc$no) {
    cosgraph_cc <- induced_subgraph(cosgraph, names(which(cc$membership==membIdx)))
    fullccs <- igraph::max_cliques(cosgraph_cc, min = nchains, max = nchains)
    for (fullcc in fullccs) {
      collectRobustComps <- collectRobustComps + induced_subgraph(cosgraph, fullcc$name)
    }
  }
  decomposeRobust <- decompose.graph(collectRobustComps)
  return(decomposeRobust)
}


###########################
# CRITERION FOR SPLITTING #
###########################
getNsplits <- function(x,thres=0.01) {
  nret <- as.numeric(names(x))
  while(length(x)>1) {
    test <- XNomial::xmulti(obs = as.vector(x), expr= rep(1/length(x), length(x)), detail=0)
    if (test$pProb<thres) {
      x <- x[names(x) != names(which.min(x))]
      nret <- as.numeric(names(x))
      } else {
      break
    }
  }
  return(nret)
}


check_split <- function(oc, thres=0.01, verbose=FALSE) {
  cc_occurences <- table(table(stringr::str_split_fixed(V(oc)$name,"_",2)[,2]))
  if (length(cc_occurences)>1) {
      subs <- getNsplits(x=cc_occurences)
      if (verbose) {
        cat(paste0("Keeping ", paste(subs,collapse="-"), " sub-components\n"))
      }
      return(subs)
    } else {
      if (verbose) {
        cat("No sub-components found\n")
      }
      return(as.numeric(names(cc_occurences)))
  }
}


##############################
# DETAILED EXTRACT FUNCTIONS #
##############################
extract_detail_components <- function(overall_components, verbose=FALSE) {
  detailed_components <- list()
  mi <- 1
  for (memb in overall_components) {
    detailed_components[[mi]] <- list()
    dosplit <- check_split(oc=memb, verbose=verbose)
    # fullccs <- igraph::max_cliques(memb, min = nchains, max = nchains)
    labels_ccs <- stringr::str_split_fixed(V(memb)$name,"_",2)
    tab_labels_ccs <- table(labels_ccs[,2])
    inccs <- names(tab_labels_ccs[tab_labels_ccs %in% c(dosplit)])
    subfullccs <- decompose.graph(induced_subgraph(memb, V(memb)$name[labels_ccs[,2] %in% inccs]))
    stopifnot(length(subfullccs)==min(dosplit))
    mii <- 1
    for (fullcc in subfullccs) {
      detailed_components[[mi]][[mii]] <- V(fullcc)$name
      mii <- mii + 1
    }
    mi <- mi + 1
  }
  return(detailed_components)
}


extract_final_components <- function(Ccs) {
  Fcs <- list()
  mi <- 1
  for (idx in 1:length(Ccs)) {
    for (Cc in Ccs[[idx]]) {
      Fcs[[mi]] <- Cc
      mi <- mi + 1
    }
  }
  return(Fcs)
}

extract_refined_mostexclusive_components_new <- function(infocc, whichin, rawCompl, fSizCompl, features=NULL, distOpt="fnchd", choice_opt=3) {
  final_weights <- data.frame(matrix(0, nrow=ncol(rawCompl[[1]][[1]]), ncol=1))
  if (is.null(features)) {
    features <- 1:ncol(rawCompl[[1]][[1]])
  }
  iccIdx <- 1
  SumAllDataItems <- c()
  # icc <- infocc[[1]]
  for ( icc in infocc ) {
    DataItems <- 0
    ticc <- stringr::str_split_fixed(icc, "_", 2)
    ticcSum <- matrix(ncol=length(features),nrow=0)
    for ( uticc in unique(ticc[,2]) ) {
      # uticc <- unique(ticc[,2])[1]
      x <- as.numeric(strsplit(uticc,"_")[[1]])

      comps_ticc <- as.numeric(ticc[ticc[,2]==uticc,,drop=FALSE][,1])
      comps_ticc_conv <- whichin[[x[[1]]]][comps_ticc]

      mergedCounts <- Reduce("+", rawCompl[[x[[1]]]][comps_ticc_conv])
      DataItems <- DataItems + sum(mergedCounts)
      mergedSizes <- Reduce("+", fSizCompl[[x[[1]]]][comps_ticc_conv])
      # mus_comptmp <- apply(mergedCounts,2, function(x){t<-table(x)
      #                                                  as.numeric(names(t)[which.max(t)])})
      mus_comptmp <- apply(mergedCounts,2, function(x){mean(x)})
      if (choice_opt==1) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          mstmp <- as.vector(size_of_feats) # + 1 plus1
          mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
          multn_params <- mus_comptmp/sum(mus_comptmp)
          n_comptmp <- sum(mus_comptmp)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
        } else if (choice_opt==2) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          mstmp <- unlist(lapply(size_of_feats, function(x) max(x, max(mus_comptmp)))) # + 1 plus1
          mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
          multn_params <- mus_comptmp/sum(mus_comptmp)
          n_comptmp <- sum(mus_comptmp)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
        } else if (choice_opt==3) {
          size_of_feats <- ceiling(rep(mean(rowSums(mergedSizes>0)), length(features)))
          mstmp <- size_of_feats # + 1 plus1
          stopifnot(!any(mus_comptmp>size_of_feats[1]))
          mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
          multn_params <- mus_comptmp/sum(mus_comptmp)
          n_comptmp <- sum(mus_comptmp)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
        } else if (choice_opt==4) {
          xtots <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          ntot <- ncol(mergedSizes)
          nobs <- mean(rowSums(mergedSizes>0))
          xexps <- xtots*nobs/ntot
          mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          nexps <- (xexps*nobs)/mus_comptmp
          mus_comptmp <- ifelse(nexps==Inf, 1e-20, nexps)
          mus_comptmp <- as.vector(xexps)
          mstmp <- as.vector(nexps) # + 1 plus1
          mstmp[mstmp==Inf] <- max(mstmp[is.finite(mstmp)])
          multn_params <- mus_comptmp/sum(mus_comptmp)
          n_comptmp <- sum(mus_comptmp)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
        } else if (choice_opt==5) {
          multn_params <- mus_comptmp/sum(mus_comptmp)
          fnchd_params  <- multn_params/median(multn_params[multn_params>0])
        } else if (choice_opt==6) {
          xtots <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          ntot <- ncol(mergedSizes)
          nobs <- mean(rowSums(mergedSizes>0))
          xexps <- xtots*nobs/ntot
          mus_comptmp <- mus_comptmp/xexps
          multn_params <- mus_comptmp/sum(mus_comptmp)
          fnchd_params  <- multn_params/median(multn_params[multn_params>0])
        } else if (choice_opt==7) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          mstmp <- as.vector(size_of_feats) # + 1 plus1
          # mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          # mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
          mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
          multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
          n_comptmp <- sum(mus_comptmp_norm)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=mstmp, n=n_comptmp, precision=0.1))
          mus_comptmp <- mus_comptmp_norm
        } else if (choice_opt==8) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          mstmp <- as.vector(size_of_feats) # + 1 plus1
          mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
          multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
          n_comptmp <- sum(mus_comptmp_norm)
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=rep(1,length(mstmp)), n=n_comptmp, precision=0.1))
          mus_comptmp <- mus_comptmp_norm
        } else if (choice_opt==9) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          # whichnotin <- !(1:length(raw_comp2feature[[x[1]]]))%in% comps_ticc
          # mergedCountsOpp <- Reduce("+", raw_comp2feature[[x[1]]][whichnotin])
          # mus_comptmpOpp <- apply(mergedCountsOpp,2, function(x){mean(x)})
          # multn_params_opp-mus_comptmpOpp
          mstmp <- as.vector(size_of_feats) # + 1 plus1
          # mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
          # mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
          mus_comptmpOpp <- mstmp-mus_comptmp
          mus_comptmp <- mus_comptmpOpp
          multn_params_opp <- (mus_comptmpOpp)/sum(mus_comptmpOpp)
          n_comptmp <- sum(mus_comptmpOpp)
          if (n_comptmp<1) {
            fnchd_params <- rep(0,length(mus_comptmpOpp))
            multn_params <- fnchd_params
            } else {
            fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmpOpp, m=mstmp, n=n_comptmp, precision=0.1))
            multn_params <- multn_params_opp
          }
        } else if (choice_opt==10) {
          size_of_feats <- unique(Reduce("+", rawCompl[[x[[1]]]]))
          mstmp <- as.vector(size_of_feats) # + 1 plus1
          mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
          mus_comptmp <- mus_comptmp_norm
          multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
          n_comptmp <- sum(mus_comptmp_norm)
          ntot <- ncol(mergedSizes)
          nobs <- mean(rowSums(mergedSizes>0))
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=ceiling(mstmp*nobs/ntot), n=n_comptmp, precision=0.1))
          multn_params <- mus_comptmp_norm
        } else {
          stop("ERROR!")
      }

      if (sum(mus_comptmp!=0)==1) {
          featuresOdds <- rep(0, length(mus_comptmp))
          featuresOdds[mus_comptmp!=0] <- Inf
          ticcSum <- rbind(ticcSum, featuresOdds)
        } else {
          featuresOdds <- switch(distOpt,
                                 "multinomial" = multn_params,
                                 "fnchd" = fnchd_params,
                                 "other" = stop("ERROR: distribution error (only FNCHD or multinomial)"))
          ticcSum <- rbind(ticcSum, featuresOdds)
      }
    }
    # all(is.finite(ticcSum))
    # all(!is.nan(ticcSum))
    # ticcSum[!is.finite(ticcSum)] <- .Machine$double.xmax
    if ((sum(ticcSum)>0) & (distOpt=="fnchd") & (sum(colSums(ticcSum==1)==nrow(ticcSum))!=1)) {
      featRef <- which.max(colSums(ticcSum==1))
      ticcSum <- ticcSum[(ticcSum[,featRef]>0) & (is.finite(ticcSum[,featRef])),,drop=FALSE]
      ticcSum <- ticcSum*1./ticcSum[,featRef]
    }
    stopifnot(all(!is.nan(ticcSum)))
    # which(colSums(ticcSum==0)!=nrow(ticcSum))
    for (widx in 1:(dim(ticcSum)[2])) {
      final_weights[widx,iccIdx] <- median(ticcSum[,widx])
    }
    # features[rev(order(final_weights[,1]))]
    # for (widx in rev(order(final_weights))) {
    #   hist(ticcSum[,widx][is.finite(ticcSum[,widx])],50, main=features[widx])
    # }
    rownames(final_weights) <- features
    if ((distOpt=="multinomial") & (sum(final_weights[,iccIdx])!=0)) {
      final_weights[,iccIdx] <- final_weights[,iccIdx]/sum(final_weights[,iccIdx])
    }
    SumAllDataItems <- c(SumAllDataItems, DataItems)
    iccIdx <- iccIdx + 1
  }
  final_weights <- final_weights[,rev(order(SumAllDataItems)),drop=FALSE]
  colnames(final_weights) <- 1:ncol(final_weights)
  return(final_weights)
}

# chain <- chlist[[1]]
getCompsFromChain2 <- function(chain) {
  if (is.null(chain)) return(NULL)
  actual_numclust <- unlist(Map(ncol, clust_categ_counts(chain)))
  convChains_cond <- which(actual_numclust==names(which.max(table(actual_numclust))))
  ccc_conv <- clust_categ_counts(chain)[convChains_cond]
  ccc_conv <- lapply(ccc_conv, function(x) {
    x <- t(x)
    return(x)
  })
  ccc_ncomp <- unique(unlist(Map(nrow,ccc_conv)))
  stopifnot(length(ccc_ncomp)==1)
  ccc_convComp <- list()
  for (compIdx in 1:ccc_ncomp) {
    ccc_convComp[[compIdx]] <- t(sapply(ccc_conv, function(x) {
      return(x[compIdx,])
    }))
  }
  return(ccc_convComp)
}

components_structure <- function (actual_feats, features) {
  # nozeros <- actual_feats != 0
  nozeros <- !is.nan(actual_feats)
  feats_nozeros <- features[nozeros]
  actual_feats <- actual_feats[nozeros]
  return(feats_nozeros[rev(order(actual_feats))])
}

driver_list <- function(cs, features) {
  str_ccs <- Map(function(x) components_structure(x,features), cs)
  outDrivers <- data.frame(matrix(NA, ncol=length(str_ccs), nrow=max(unlist(Map(length,str_ccs)))))
  for (idx in 1:length(str_ccs)) {
    str_cc <- str_ccs[[idx]]
    outDrivers[1:length(str_cc),idx] <- str_cc
  }
  outDrivers[is.na(outDrivers)] <- ""
  if ("0" %in% colnames(cs)) {
      colnames(outDrivers) <- paste("Component",(1:ncol(outDrivers))-1,sep="_")
    } else {
      colnames(outDrivers) <- paste("Component",1:ncol(outDrivers),sep="_")
  }
  return(outDrivers)
}

extract_cdcsize <- function(chain, cdc_idxs) {
  if (is.null(chain)) return(NULL)
  n_features <- nrow(chain@clust_categ_counts[[1]])
  actual_numclust <- unlist(Map(ncol, clust_categ_counts(chain)))
  convChains_cond <- which(actual_numclust==names(which.max(table(actual_numclust))))
  cdc_conv <- clust_dp_counts(chain)[convChains_cond]
  cdc_size <- lapply(cdc_conv, function(x) {
    return(as.matrix(x)[cdc_idxs,])
  })
  cdc_ncomp <- unique(unlist(Map(ncol,cdc_conv)))
  cdc_convComp <- list()
  for (compIdx in 1:cdc_ncomp) {
    cdc_convComp[[compIdx]] <- t(sapply(cdc_size, function(x) {
      return(x[,compIdx])
    }))
  }
  return(cdc_convComp)
}

getNormalizingFactors <- function(params, nMax) {
  clust2dup <- as.data.frame(t(params))
  colnames(clust2dup) <- NULL
  rownames(clust2dup) <- NULL

  # outfile <- "/home/PERSONALE/daniele.dallolio3/HARMONY_AML/tmp_Weights.csv"
  outfile <- "/home/PERSONALE/daniele.dallolio3/temp/tmp_Weights.csv"
  data.table::fwrite(clust2dup, outfile)
  # filep0s <- "/home/PERSONALE/daniele.dallolio3/HARMONY_AML/temp_logP0s.csv"
  filep0s <- "/home/PERSONALE/daniele.dallolio3/temp/temp_logP0s.csv"
  cmdline <- paste0("/home/PERSONALE/daniele.dallolio3/FisherNonCentralHyper/P0FisherNonCentral_3.out -w ", outfile," -o ", filep0s," -n ", nMax," -t 24")
  system(cmdline)
  normalizing_constants <- as.matrix(data.table::fread(filep0s, header=FALSE))
  colnames(normalizing_constants) <- NULL
  system(paste0("rm ",outfile))
  system(paste0("rm ",filep0s))
  return(normalizing_constants)
}

getGarbageComponent <- function(whichin, rawCompl, fSizCompl, features=NULL, distOpt="fnchd", choice_opt=3) {
  whichout <- lapply(1:length(rawCompl), function(x) setdiff(1:length(rawCompl[[x]]),whichin[[x]]) )
  ticcSum <- matrix(ncol=length(features),nrow=0)
  for (tabo in 1:length(rawCompl)) {
    comps_ticc <- whichout[[tabo]]
    mergedCounts <- Reduce("+", rawCompl[[tabo]][comps_ticc])
    if (sum(mergedCounts)==0) {
      next
    }
    mergedSizes <- Reduce("+", fSizCompl[[tabo]][comps_ticc])
    # mus_comptmp <- apply(mergedCounts,2, function(x){t<-table(x)
    #                                                  as.numeric(names(t)[which.max(t)])})
    mus_comptmp <- apply(mergedCounts,2, function(x){mean(x)})
    if (choice_opt==1) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        mstmp <- as.vector(size_of_feats) # + 1 plus1
        mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
        multn_params <- mus_comptmp/sum(mus_comptmp)
        n_comptmp <- sum(mus_comptmp)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
      } else if (choice_opt==2) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        mstmp <- unlist(lapply(size_of_feats, function(x) max(x, max(mus_comptmp)))) # + 1 plus1
        mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
        multn_params <- mus_comptmp/sum(mus_comptmp)
        n_comptmp <- sum(mus_comptmp)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
      } else if (choice_opt==3) {
        size_of_feats <- ceiling(rep(mean(rowSums(mergedSizes>0)), length(features)))
        mstmp <- size_of_feats # + 1 plus1
        stopifnot(!any(mus_comptmp>size_of_feats[1]))
        mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
        multn_params <- mus_comptmp/sum(mus_comptmp)
        n_comptmp <- sum(mus_comptmp)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
      } else if (choice_opt==4) {
        xtots <- unique(Reduce("+", rawCompl[[tabo]]))
        ntot <- ncol(mergedSizes)
        nobs <- mean(rowSums(mergedSizes>0))
        xexps <- xtots*nobs/ntot
        mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        nexps <- (xexps*nobs)/mus_comptmp
        mus_comptmp <- ifelse(nexps==Inf, 1e-20, nexps)
        mus_comptmp <- as.vector(xexps)
        mstmp <- as.vector(nexps) # + 1 plus1
        mstmp[mstmp==Inf] <- max(mstmp[is.finite(mstmp)])
        multn_params <- mus_comptmp/sum(mus_comptmp)
        n_comptmp <- sum(mus_comptmp)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp, m=mstmp, n=n_comptmp, precision=0.1))
      } else if (choice_opt==5) {
        multn_params <- mus_comptmp/sum(mus_comptmp)
        fnchd_params  <- multn_params/median(multn_params[multn_params>0])
      } else if (choice_opt==6) {
        xtots <- unique(Reduce("+", rawCompl[[tabo]]))
        ntot <- ncol(mergedSizes)
        nobs <- mean(rowSums(mergedSizes>0))
        xexps <- xtots*nobs/ntot
        mus_comptmp <- mus_comptmp/xexps
        multn_params <- mus_comptmp/sum(mus_comptmp)
        fnchd_params  <- multn_params/median(multn_params[multn_params>0])
      } else if (choice_opt==7) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        mstmp <- as.vector(size_of_feats) # + 1 plus1
        # mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        # mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
        mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
        multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
        n_comptmp <- sum(mus_comptmp_norm)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=mstmp, n=n_comptmp, precision=0.1))
      } else if (choice_opt==8) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        mstmp <- as.vector(size_of_feats) # + 1 plus1
        mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
        multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
        n_comptmp <- sum(mus_comptmp_norm)
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=rep(1,length(mstmp)), n=n_comptmp, precision=0.1))
      } else if (choice_opt==9) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        # whichnotin <- !(1:length(raw_comp2feature[[x[1]]]))%in% comps_ticc
        # mergedCountsOpp <- Reduce("+", raw_comp2feature[[x[1]]][whichnotin])
        # mus_comptmpOpp <- apply(mergedCountsOpp,2, function(x){mean(x)})
        # multn_params_opp-mus_comptmpOpp
        mstmp <- as.vector(size_of_feats) # + 1 plus1
        # mus_comptmp <- ifelse(mus_comptmp>1e-10, mus_comptmp-1e-10, 0)
        # mus_comptmp <- ifelse(mus_comptmp==0, 1e-20, mus_comptmp)
        mus_comptmpOpp <- mstmp-mus_comptmp
        multn_params_opp <- (mus_comptmpOpp)/sum(mus_comptmpOpp)
        n_comptmp <- sum(mus_comptmpOpp)
        if (n_comptmp<1) {
          fnchd_params <- rep(0,length(mus_comptmpOpp))
          multn_params <- fnchd_params
          } else {
          fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmpOpp, m=mstmp, n=n_comptmp, precision=0.1))
          multn_params <- multn_params_opp
        }
        multn_params <- multn_params_opp
      } else if (choice_opt==10) {
        size_of_feats <- unique(Reduce("+", rawCompl[[tabo]]))
        mstmp <- as.vector(size_of_feats) # + 1 plus1
        mus_comptmp_norm <- apply(mergedCounts/rowSums(mergedSizes>0),2, function(x){mean(x)})
        multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
        n_comptmp <- sum(mus_comptmp_norm)
        ntot <- ncol(mergedSizes)
        nobs <- mean(rowSums(mergedSizes>0))
        fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=ceiling(mstmp*nobs/ntot), n=n_comptmp, precision=0.1))
        multn_params <- multn_params_opp
      } else {
        stop("ERROR!")
    }

    if (sum(mus_comptmp!=0)==1) {
        featuresOdds <- rep(0, length(mus_comptmp))
        featuresOdds[mus_comptmp!=0] <- Inf
        ticcSum <- rbind(ticcSum, featuresOdds)
      } else {
        featuresOdds <- switch(distOpt,
                               "multinomial" = multn_params,
                               "fnchd" = fnchd_params,
                               "other" = stop("ERROR: distribution error (only FNCHD or multinomial)"))
        ticcSum <- rbind(ticcSum, featuresOdds)
      }
  }
  if (nrow(ticcSum)<1) {
    return(NULL)
  }
  if ((sum(ticcSum)>0) & (distOpt=="fnchd") & (sum(colSums(ticcSum==1)==nrow(ticcSum))!=1)) {
    featRef <- which.max(colSums(ticcSum==1))
    ticcSum <- ticcSum[(ticcSum[,featRef]>0) & (is.finite(ticcSum[,featRef])),,drop=FALSE]
    ticcSum <- ticcSum*1./ticcSum[,featRef]
  }
  stopifnot(all(!is.nan(ticcSum)))

  tmp <- apply(ticcSum, 2, median)
  final_weights <- data.frame("V1"=tmp)
  colnames(final_weights) <- c(0)
  rownames(final_weights) <- features
  if ((distOpt=="multinomial") & (sum(tmp)!=0)) {
    final_weights[,1] <- final_weights[,1,drop=FALSE]/sum(tmp)
  }
  return(final_weights)
}



extract_fnchd_components <- function(chlist, n_subjects, features=NULL, distOpt="fnchd", verbose=FALSE, choice_opt=3, garbageOpt=FALSE) {
  subjs_dps <- 1:n_subjects+(nrow(chlist[[1]]@clust_dp_counts[[1]])-n_subjects)
  if (is.null(features)) {
    features <- 1:nrow(chlist[[1]]@clust_categ_counts[[1]])
  }
  beforecompslist <- lapply(chlist, extractConvergedComponents)
  compslist <- lapply(beforecompslist, function(x) x[[1]])
  compswhichkept <- lapply(beforecompslist, function(x) x[[2]])
  crit0_conv <- unlist(Map(function(x) nrow(x)>0, compslist))
  if (all(!crit0_conv)) {
    return(NULL)
  }
  # unlist(lapply(chlist, function(x) (max(table(x@numcluster))/1e4)))
  # table(chlist[[10]]@numcluster)
  crit1_conv  <- unlist(lapply(chlist, function(x) (max(table(x@numcluster))/1e4)>0.))
  explained_multinomials <- lapply(compslist, function(x) min(rowSums(x)))
  crit2_thres <- coda::HPDinterval(coda::as.mcmc(unlist(explained_multinomials)))[1]
  crit2_conv <- !unlist(lapply(explained_multinomials, function(x) any(x<crit2_thres) ))
  crit_conv <- crit1_conv & crit2_conv
  linknet <- soft_dist(compsObjList=compslist)
  additional_unconverged <- rm_toNA(cosgraph=linknet, K=length(compslist))
  crit_conv[additional_unconverged] <- FALSE
  linknet <- rm_unconverged(cosgraph=linknet, crit_conv=crit_conv)
  linknet <- rm_iso(linknet)
  nchains <- length(unique(stringr::str_split_fixed(V(linknet)$name, "_", 2)[,2]))
  ocs <- extract_overall_components(cosgraph=linknet, nchains=nchains)
  if (length(ocs)<1) {
    return(NULL)
  }
  # cat("Number of overall components: ", length(ocs), "\n")
  comps_in <- unique(unlist(Map(function(x) V(x)$name, ocs)))
  comps_out <- V(linknet)$name[!V(linknet)$name %in% comps_in]
  prop.ex <- 1 - sum(V(linknet)$name %in% comps_out)/length(V(linknet))
  # cat(paste0("Proportion of components preserved: ", round(prop.ex,3)*100, "%\n"))
  dcs <- extract_detail_components(overall_components=ocs, verbose=verbose)
  focs <- extract_final_components(dcs)
  raw_comp2feature_complete <- lapply(chlist, getCompsFromChain2)
  raw_comp2feature <- lapply(1:length(chlist), function(x) raw_comp2feature_complete[[x]][compswhichkept[[x]]])
  raw_featuresSize_complete <- lapply(chlist, function(ch) extract_cdcsize(ch, subjs_dps))
  raw_featuresSize <- lapply(1:length(chlist), function(x) raw_featuresSize_complete[[x]][compswhichkept[[x]]])
  # featuresSize <- lapply(chlist, function(ch) {
  #   rowSums(clust_categ_counts(ch)[[1]])
  # })
  averages_comps_raw <- extract_refined_mostexclusive_components_new(infocc = focs,
                                                                     whichin = compswhichkept,
                                                                     rawCompl = raw_comp2feature_complete,
                                                                     fSizCompl = raw_featuresSize_complete,
                                                                     features = features,
                                                                     distOpt = distOpt,
                                                                     choice_opt = choice_opt)
  if (garbageOpt) {
    garbage_component <- getGarbageComponent(whichin = compswhichkept,
                                             rawCompl = raw_comp2feature_complete,
                                             fSizCompl = raw_featuresSize_complete,
                                             features = features,
                                             distOpt = distOpt,
                                             choice_opt = choice_opt
                                           )
    if (!is.null(garbage_component)) {
      averages_comps_raw <- cbind(garbage_component, averages_comps_raw)
    }
  }

  if (distOpt=="fnchd") {
      averages_comps <- as.data.frame(apply(averages_comps_raw, 2, function(x) {
        if (sum(x)==0) {
            x
          } else {
            x/median(x[x!=0 & x!=Inf])
        }
      }))
    } else {
      averages_comps <- averages_comps_raw
  }
  components_overview <- driver_list(cs=averages_comps, features=features)
  components_details <- lapply(1:ncol(averages_comps), function(x) {
    xc <- averages_comps[, x, drop=FALSE]
    # xc <- xc[rowSums(xc)>0,, drop=FALSE]
    xc <- xc[rev(order(xc[,1])),,drop=FALSE]
    return(xc)
  })
  outObj <- list("overview"=components_overview, "parameters"=components_details)
  if (distOpt=="fnchd") {
    ExclusiveDefiningLesions <- features[which(rowSums(averages_comps==Inf)>0)]
    # stopifnot(all(unique(as.vector(as.matrix(averages_comps[ExclusiveDefiningLesions,]))) %in% c(0,Inf)))
    maxMut <- max(rowSums(as.matrix(clust_dp_counts(chlist[[1]])[[1]])[subjs_dps,]))
    if (sum(averages_comps)!=0) {
      logP0s <- getNormalizingFactors(averages_comps[!rownames(averages_comps) %in% ExclusiveDefiningLesions,], maxMut)
      logP0s <- as.data.frame(logP0s)
      } else {
      logP0s <- data.frame(matrix(0,nrow=1,ncol=maxMut))
    }
    colnames(logP0s) <- NULL
    outObj$logP0s <- logP0s
  }
  return(outObj)
}


#%%
############
# ANALYSIS #
############
# typePriorGens <- c("1/n", "1")
# typePriorHDPs <- c("1/n", "1")
# idxGen <- 2
# idxHDP <- 2
# res <- data.frame()
# for (idxGen in 1:2) {
#   message(idxGen)
#   for (idxHDP in 1:2) {
#     message(idxHDP)
#     for (n_components in c(1,5,10)) {
#       for (mean_freq_rate in 1:10) {
#         # mean_freq_rate <- 2
#         print(mean_freq_rate)
#         n_people       <- 1000
#         n_features     <- 50
#         n_components   <- as.numeric(n_components)
#         set.seed(1)
#         #######################################################################
#         # SIMULATE MIXTURE OF FISHER NON-CENTRAL HYPERGEOMETRIC DISTRIBUTIONS #
#         #######################################################################
#         freq_rate_all <- rpois(n_people, lambda=mean_freq_rate)
#
#
#         subj2add <- BiasedUrn::rMFNCHypergeo(nran = 1,
#                                              m = rep(1, n_features),
#                                              n = freq_rate_all[idx],
#                                              odds = odds_vec_all[true_clusters[idx],])
#         rev(order(suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=colSums(data[true_clusters==3,]), m=rep(sum(true_clusters==3), ncol(data)), n=sum(colSums(data[true_clusters==3,])), precision=0.1))))
#         rev(order(ws_vec_all[4,]))
#         rev(sort(ws_vec_all[4,]))
#
#         if (typePriorGens[idxGen]=="1/n") {
#             ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
#                                                        alpha=rep(1/n_features,n_features))
#             dirtest <- "alpha_1_over_n_features/"
#           } else {
#             ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
#                                                        alpha=rep(1.,n_features))
#             dirtest <- "alpha_1s/"
#         }
#
#         odds_vec_all  <- ws_vec_all/apply(ws_vec_all, 1, median)
#         # t(apply(odds_vec_all, 1, function(x) head(rev(order(x))) ))
#         composition   <- LaplacesDemon::rdirichlet(n=1, alpha=rep(1,n_components))
#         true_clusters <- rep(0, n_people)
#         true_clusters[freq_rate_all!=0] <- sample(seq(1,n_components),
#                                                   n_people-sum(freq_rate_all==0),
#                                                   replace=TRUE,
#                                                   prob=composition)
#         simdata <- data.frame()
#         for (idx in 1:length(freq_rate_all)) {
#           if (freq_rate_all[idx]==0) {
#             stopifnot(true_clusters[idx]==0)
#             subj2add <- rep(0, n_features)
#             } else {
#             subj2add <- BiasedUrn::rMFNCHypergeo(nran = 1,
#                                                  m = rep(1, n_features),
#                                                  n = freq_rate_all[idx],
#                                                  odds = odds_vec_all[true_clusters[idx],])
#           }
#           simdata <- rbind(simdata, subj2add)
#         }
#         table(true_clusters)
#         median(rowSums(simdata))
#         #################################
#         # FILTER OUT ALL ZEROS SUBJECTS #
#         #################################
#         dim(simdata)
#         data <- simdata[,colSums(simdata)>1]
#         dim(data)
#         data <- data[rowSums(data)>0,]
#         features4components <- as.character(as.vector(which(colSums(simdata)>1)))
#
#         ###########################
#         # DIRICHLET PROCESS SETUP #
#         ###########################
#         n_features <- ncol(data)
#         n_patients <- nrow(data)
#
#         files_dirtest <- list.files(dirtest)
#         files_dirtest <- files_dirtest[files_dirtest!="results"]
#         conf_used <- stringr::str_split_fixed(files_dirtest, "_", 6)[,c(3,4,5,6)]
#         conf_used <- conf_used[conf_used[,2]==ncol(simdata),]
#
#         if (typePriorHDPs[idxHDP]=='1') {
#             simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="dp1s",,drop=FALSE]
#           } else {
#             simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="",,drop=FALSE]
#         }
#         n_people       <- simPrecomputed[,1]
#         n_features     <- simPrecomputed[,2]
#         n_components   <- simPrecomputed[,3]
#         type_dp        <- simPrecomputed[,4]
#
#         pathdir <- paste0(dirtest,"simulated_data_", n_people, "_", n_features, "_", n_components)
#         pathdir <- ifelse(type_dp!="", paste0(pathdir,"_",type_dp), pathdir)
#         logdir <- file.path(pathdir, mean_freq_rate)
#         multichains <- vector("list", 10)
#         for (k in 1:10) {
#           outfile_chain_folder <- file.path(logdir, "chains")
#           allfiles <- list.files(outfile_chain_folder)
#           fl <- allfiles[startsWith(allfiles, paste0("HDMM_it_", k, "_"))]
#           stopifnot(length(fl)==1)
#           multichains[[k]] <- readRDS(file.path(outfile_chain_folder, fl))
#         }
#
#         whc_freqs <- c()
#         dom_freqs <- c()
#         for (k in 1:10) {
#           tabnumcl <- table(multichains[[k]]@numcluster)
#           whc_freqs <- c(whc_freqs, names(which.max(tabnumcl)))
#           dom_freqs <- c(dom_freqs, max(tabnumcl)/length(multichains[[k]]@numcluster))
#         }
#         tmp4res <- data.frame("N"=n_people,
#                               "M"=n_features,
#                               "K"=n_components,
#                               "AverageOccurence"=mean_freq_rate,
#                               "DominantNumber"=whc_freqs,
#                               "DominanceFrequency"=dom_freqs,
#                               "priorHDP"=ifelse(type_dp=="", "1/n", "1"),
#                               "priorGen"=ifelse(typePriorGens[idxGen]=="1/n","1/n", "1")
#                             )
#         res <- rbind(res, tmp4res)
#       }
#     }
#   }
# }
#
# data.table::fwrite(res, "summary_conv.csv")
#
#
#
# #%%
library(ggplot2)
library(dplyr)
options(repr.plot.width=16, repr.plot.height=10)

summary_conv <- data.table::fread("summary_conv.csv")

df <- summary_conv
df$K <- factor(paste0("K==", df$K), levels=c("K==1","K==5","K==10"))
df$priorGen[df$priorGen=="1/n"] <- "1/M"
df$priorHDP[df$priorHDP=="1/n"] <- "1/M"
df$priorGen <- factor(df$priorGen, levels=c("1","1/M"))
df$DominantNumber <- as.numeric(as.character(df$DominantNumber))
df$DominanceFrequency <- as.numeric(as.character(df$DominanceFrequency))
df$AverageOccurence <- factor(df$AverageOccurence, levels=1:10)
df$priorHDP <- paste0("alpha[HDMM]==", df$priorHDP)
df$priorGen <- paste0("alpha[sim]==", df$priorGen)
df$priorHDP <- factor(df$priorHDP, levels=c("alpha[HDMM]==1","alpha[HDMM]==1/M"))
df$priorGen <- factor(df$priorGen, levels=c("alpha[sim]==1","alpha[sim]==1/M"))
domplot <- ggplot(df[df$K %in% c("K==5","K==10")], aes(x=AverageOccurence, y=log10(DominantNumber), colour=priorHDP)) +
           geom_jitter() +
           geom_hline(aes(yintercept=log10(as.numeric(as.character(gsub("K==", "", K))))), color = "gold", size=1.5) +
           scale_color_discrete(labels=c('alpha[HDMM]==1'=expression(alpha[HDMM]==1),'alpha[HDMM]==1/M'=expression(alpha[HDMM]==1/M))) +
           xlab("Average number of events per sample") +
           ylab("Estimated number of components (log scale)") +
           facet_wrap(~K+priorGen, nrow = 3, scales = "free_y", labeller=label_parsed) +
           theme(legend.position="top", text = element_text(size = 20), legend.title=element_blank())
ggsave(file.path("PAPER/", paste0("domplot.eps")), print(domplot), width=400,height=280,units="mm",dpi=600)

freplot <- ggplot(df[df$K %in% c("K==5","K==10")], aes(x=AverageOccurence, y=DominanceFrequency*100, colour=priorHDP)) +
           geom_jitter() +
           scale_color_discrete(labels=c('alpha[HDMM]==1'=expression(alpha[HDMM]==1),'alpha[HDMM]==1/M'=expression(alpha[HDMM]==1/M))) +
           xlab("Average number of events per sample") +
           ylab("Frequency of estimated number of cluster (%)") +
           facet_wrap(~K+priorGen, nrow = 3, scales = "free_y", labeller=label_parsed) +
           theme(legend.position="top", text = element_text(size = 20), legend.title=element_blank())
ggsave(file.path("PAPER/", paste0("freplot.eps")), print(freplot), width=400,height=280,units="mm",dpi=600)


#%%
###########
# GET P0s #
###########
typePriorGens <- c("1/n", "1")
typePriorHDPs <- c("1/n", "1")
choice_opt <- 7
gOpt <- TRUE
for (idxGen in 1:2) {
  message(idxGen)
  for (idxHDP in 1:2) {
    message(idxHDP)
    for (n_components in c(1,5,10)) {
      for (mean_freq_rate in 1:10) {
        # mean_freq_rate <- 4
        # n_components <- 5
        print(mean_freq_rate)
        n_people       <- 1000
        n_features     <- 50
        n_components   <- as.numeric(n_components)
        set.seed(1)
        #######################################################################
        # SIMULATE MIXTURE OF FISHER NON-CENTRAL HYPERGEOMETRIC DISTRIBUTIONS #
        #######################################################################
        freq_rate_all <- rpois(n_people, lambda=mean_freq_rate)
        if (typePriorGens[idxGen]=="1/n") {
            ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                       alpha=rep(1/n_features,n_features))
            dirtest <- "alpha_1_over_n_features/"
          } else {
            ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                       alpha=rep(1.,n_features))
            dirtest <- "alpha_1s/"
        }

        odds_vec_all  <- ws_vec_all/apply(ws_vec_all, 1, median)
        # t(apply(odds_vec_all, 1, function(x) head(rev(order(x))) ))
        composition   <- LaplacesDemon::rdirichlet(n=1, alpha=rep(1,n_components))
        true_clusters <- rep(0, n_people)
        true_clusters[freq_rate_all!=0] <- sample(seq(1,n_components),
                                                  n_people-sum(freq_rate_all==0),
                                                  replace=TRUE,
                                                  prob=composition)
        simdata <- data.frame()
        for (idx in 1:length(freq_rate_all)) {
          if (freq_rate_all[idx]==0) {
            stopifnot(true_clusters[idx]==0)
            subj2add <- rep(0, n_features)
            } else {
            subj2add <- BiasedUrn::rMFNCHypergeo(nran = 1,
                                                 m = rep(1, n_features),
                                                 n = freq_rate_all[idx],
                                                 odds = odds_vec_all[true_clusters[idx],])
          }
          simdata <- rbind(simdata, subj2add)
        }
        table(true_clusters)
        median(rowSums(simdata))
        #################################
        # FILTER OUT ALL ZEROS SUBJECTS #
        #################################
        dim(simdata)
        data <- simdata[,colSums(simdata)>1]
        dim(data)
        data <- data[rowSums(data)>0,]
        features4components <- as.character(as.vector(which(colSums(simdata)>1)))

        ###########################
        # DIRICHLET PROCESS SETUP #
        ###########################
        n_features <- ncol(data)
        n_patients <- nrow(data)

        files_dirtest <- list.files(dirtest)
        files_dirtest <- files_dirtest[files_dirtest!="results"]
        conf_used <- stringr::str_split_fixed(files_dirtest, "_", 6)[,c(3,4,5,6)]
        conf_used <- conf_used[conf_used[,2]==ncol(simdata),]

        if (typePriorHDPs[idxHDP]=='1') {
            simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="dp1s",,drop=FALSE]
          } else {
            simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="",,drop=FALSE]
        }
        n_people       <- simPrecomputed[,1]
        n_features     <- simPrecomputed[,2]
        n_components   <- simPrecomputed[,3]
        type_dp        <- simPrecomputed[,4]

        pathdir <- paste0(dirtest,"simulated_data_", n_people, "_", n_features, "_", n_components)
        pathdir <- ifelse(type_dp!="", paste0(pathdir,"_",type_dp), pathdir)
        logdir <- file.path(pathdir, mean_freq_rate)

        if (gOpt) {
            outfile_multn_folder <- file.path(logdir, paste0("multinomials_garbage_", choice_opt))
            outfile_fnchd_folder <- file.path(logdir, paste0("fnchds_garbage_", choice_opt))
          } else {
            outfile_multn_folder <- file.path(logdir, paste0("multinomials_", choice_opt))
            outfile_fnchd_folder <- file.path(logdir, paste0("fnchds_", choice_opt))
        }
        if (file.exists(file.path(outfile_multn_folder, "finalresult_multinomial.rds")) & file.exists(file.path(outfile_fnchd_folder, "finalresult_fnchd.rds"))) {
            next
          } else {
            multichains <- vector("list", 10)
            for (k in 1:10) {
              outfile_chain_folder <- file.path(logdir, "chains")
              allfiles <- list.files(outfile_chain_folder)
              fl <- allfiles[startsWith(allfiles, paste0("HDMM_it_", k, "_"))]
              stopifnot(length(fl)==1)
              multichains[[k]] <- readRDS(file.path(outfile_chain_folder, fl))
            }
            # Compute components
            finalresult_multinomial <- extract_fnchd_components(chlist=multichains,
                                                                n_subjects=nrow(data),
                                                                distOpt="multinomial",
                                                                features=features4components,
                                                                choice_opt=choice_opt,
                                                                garbageOpt=gOpt
                                                              )
            if (is.null(finalresult_multinomial)) {
              next
            }
            finalresult_fnchd       <- extract_fnchd_components(chlist=multichains,
                                                                n_subjects=nrow(data),
                                                                distOpt="fnchd",
                                                                features=features4components,
                                                                choice_opt=choice_opt,
                                                                garbageOpt=gOpt
                                                              )
            if (is.null(finalresult_fnchd)) {
              next
            }
            dir.create(outfile_multn_folder, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            saveRDS(finalresult_multinomial, file.path(outfile_multn_folder, "finalresult_multinomial.rds"))
            dir.create(outfile_fnchd_folder, showWarnings = FALSE, recursive = TRUE, mode = "0777")
            saveRDS(finalresult_fnchd, file.path(outfile_fnchd_folder, "finalresult_fnchd.rds"))
        }
      }
    }
  }
}


#%%
###############################
# COMPONENTS ACCURACY METRICS #
###############################

# MultinomialAssignment <- function(ps, data) {
#   ps <- t(ps)
#   Pmultinomial <- data.frame(matrix(nrow=nrow(data),ncol=nrow(ps)))
#   colnames(Pmultinomial) <- 1:ncol(Pmultinomial)
#   for (pidx in 1:nrow(Pmultinomial)) {
#     for (cidx in 1:nrow(ps)) {
#       Pmultinomial[pidx,cidx] <- stats::dmultinom(as.matrix(data[pidx,]), prob=as.matrix(ps[cidx,]))
#     }
#   }
#   ambiguous <- apply(Pmultinomial,1, function(x) sum(x=max(x))>1)
#   stopifnot(any(!ambiguous))
#   assignments <- apply(Pmultinomial,1, which.max)
#   return(assignments)
# }
MultinomialAssignmentOpp <- function(ps, data) {
  Pmultinomial <- data.frame(matrix(nrow=nrow(data),ncol=ncol(ps)))
  colnames(Pmultinomial) <- 1:ncol(Pmultinomial)
  for (cidx in 1:ncol(Pmultinomial)) {
    partial_num <- log(t(t(data)*ps[,cidx]))
    cond_to0 <- (data==1) & (partial_num==-Inf)
    # stopifnot(!any((data==1) & (partial_num==-Inf)))
    partial_num[!is.finite(partial_num)] <- 0
    partial_num[cond_to0] <- -Inf
    lognums <- rowSums(partial_num)
    Pmultinomial[,cidx] <- 1-exp(as.vector(lognums))
  }
  ambiguous <- apply(Pmultinomial,1, function(x) sum(x=max(x))>1)
  stopifnot(any(!ambiguous))
  assignments <- apply(Pmultinomial,1, function(x) as.numeric(colnames(ps)[which.max(x)]))
  return(assignments)
}


MultinomialAssignment <- function(ps, data) {
  Pmultinomial <- data.frame(matrix(nrow=nrow(data),ncol=ncol(ps)))
  colnames(Pmultinomial) <- 1:ncol(Pmultinomial)
  for (cidx in 1:ncol(Pmultinomial)) {
    partial_num <- log(t(t(data)*ps[,cidx]))
    cond_to0 <- (data==1) & (partial_num==-Inf)
    # stopifnot(!any((data==1) & (partial_num==-Inf)))
    partial_num[!is.finite(partial_num)] <- 0
    partial_num[cond_to0] <- -Inf
    lognums <- rowSums(partial_num)
    Pmultinomial[,cidx] <- as.vector(lognums)
  }
  ambiguous <- apply(Pmultinomial,1, function(x) sum(x==max(x))>1)
  stopifnot(any(!ambiguous))
  # assignments <- apply(Pmultinomial,1, function(x) as.numeric(colnames(ps)[which.max(x)]))
  assignments <- apply(Pmultinomial,1, which.max)
  return(assignments)
}

FNCHDAssignmentOpp <- function(ps, data, logP0s) {
  numMuts <- as.vector(rowSums(data))
  ALM <- data.frame(matrix(NA, nrow=nrow(data), ncol=nrow(logP0s)))
  # idx <- 4
  for (idx in 1:ncol(ALM)) {
    logP0s_idx <- as.vector(unlist(logP0s[idx,]))
    partial_num <- log(t(t(data)*ps[,idx]))
    cond_to0 <- (data==1) & (partial_num==-Inf)
    # stopifnot(!any((data==1) & (partial_num==-Inf)))
    partial_num[!is.finite(partial_num)] <- 0
    partial_num[cond_to0] <- -Inf
    lognums <- rowSums(partial_num)
    ALM[,idx] <- 1-exp(as.vector(lognums-ifelse(any(ps[,idx]==Inf), .Machine$double.xmax, logP0s_idx[numMuts])))
  }
  ambiguous <- apply(ALM,1, function(x) sum(x==max(x))>1)
  stopifnot(any(!ambiguous))
  assignments <- apply(ALM,1, function(x) as.numeric(colnames(ps)[which.max(x)]))
  return(assignments)
}


FNCHDAssignment <- function(ps, data, logP0s) {
  numMuts <- as.vector(rowSums(data))
  ALM <- data.frame(matrix(NA, nrow=nrow(data), ncol=nrow(logP0s)))
  # idx <- 2
  for (idx in 1:ncol(ALM)) {
    logP0s_idx <- as.vector(unlist(logP0s[idx,]))
    partial_num <- log(t(t(data)*ps[,idx]))
    cond_to0 <- (data==1) & (partial_num==-Inf)
    cond_Inf <- (data==1) & (partial_num==Inf)
    # stopifnot(!any((data==1) & (partial_num==-Inf)))
    partial_num[!is.finite(partial_num)] <- 0
    partial_num[cond_to0] <- -Inf
    partial_num[cond_Inf] <- Inf
    lognums <- rowSums(partial_num)
    if (any(ps[,idx]==Inf)) {
        ALM[,idx] <- as.vector(lognums-.Machine$double.xmax)
      } else {
        ALM[,idx] <- as.vector(lognums-logP0s_idx[numMuts])
    }
  }
  ambiguous <- apply(ALM,1, function(x) sum(x==max(x))>1)
  stopifnot(any(!ambiguous))
  # assignments <- apply(ALM,1, function(x) as.numeric(colnames(ps)[which.max(x)]))
  assignments <- apply(ALM,1, which.max)
  return(assignments)
}

getContingencyTable <- function(grph, onevec, twovec, onelab, twolab) {
  cc_assess <- components(grph)
  ref_vec <- rep(NA, length(onevec))
  oth_vec <- rep(NA, length(onevec))
  for (cc_tmp in sort(unique(cc_assess$membership))) {
    cc_elements <- names(cc_assess$membership)[cc_assess$membership==cc_tmp]
    cc_els_one <- as.numeric(stringr::str_split_fixed(cc_elements[endsWith(cc_elements,onelab)], "_",2)[,1])
    cc_els_two <- as.numeric(stringr::str_split_fixed(cc_elements[endsWith(cc_elements,twolab)], "_",2)[,1])
    ref_vec[onevec %in% cc_els_one] <- cc_tmp
    oth_vec[twovec %in% cc_els_two] <- cc_tmp
  }
  return(table(ref_vec, oth_vec))
}

# cm<-contMultn

custom_mcc <- function(cm){
  if (all(dim(cm)==c(1,1))) {
    return(1.)
  }
  samples <- sum(cm)
  correct <- sum(diag(cm))
  y <- colSums(cm)
  x <- rowSums(cm)
  cov_x_y <- correct * samples - x%*%y
  cov_y_y <- samples * samples - y%*%y
  cov_x_x <- samples * samples - x%*%x
  denom <- sqrt(cov_x_x * cov_y_y)
  denom <- ifelse(denom != 0.0, denom, 1.0)
  return(as.vector(cov_x_y / denom))
}


res <- data.frame()
# gOpt <- FALSE
for (choice_opt in c(1,2,3,5,6,7,8)) {
  message(choice_opt)
  for (idxGen in 1:2) {
    message(idxGen)
    for (idxHDP in 1:2) {
      message(idxHDP)
      for (n_components in c(1,5,10)) {
        for (mean_freq_rate in 1:10) {

          # n_components=10
          # mean_freq_rate=8
          # idxGen=1
          # idxHDP=1

          print(mean_freq_rate)
          n_people       <- 1000
          n_features     <- 50
          n_components   <- as.numeric(n_components)
          set.seed(1)
          #######################################################################
          # SIMULATE MIXTURE OF FISHER NON-CENTRAL HYPERGEOMETRIC DISTRIBUTIONS #
          #######################################################################
          freq_rate_all <- rpois(n_people, lambda=mean_freq_rate)

          if (typePriorGens[idxGen]=="1/n") {
              ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                         alpha=rep(1/n_features,n_features))
              dirtest <- "alpha_1_over_n_features/"
            } else {
              ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                         alpha=rep(1.,n_features))
              dirtest <- "alpha_1s/"
          }

          odds_vec_all  <- ws_vec_all/apply(ws_vec_all, 1, median)
          # t(apply(odds_vec_all, 1, function(x) head(rev(order(x))) ))
          composition   <- LaplacesDemon::rdirichlet(n=1, alpha=rep(1,n_components))
          true_clusters <- rep(0, n_people)
          true_clusters[freq_rate_all!=0] <- sample(seq(1,n_components),
                                                    n_people-sum(freq_rate_all==0),
                                                    replace=TRUE,
                                                    prob=composition)
          simdata <- data.frame()
          for (idx in 1:length(freq_rate_all)) {
            if (freq_rate_all[idx]==0) {
              stopifnot(true_clusters[idx]==0)
              subj2add <- rep(0, n_features)
              } else {
              subj2add <- BiasedUrn::rMFNCHypergeo(nran = 1,
                                                   m = rep(1, n_features),
                                                   n = freq_rate_all[idx],
                                                   odds = odds_vec_all[true_clusters[idx],])
            }
            simdata <- rbind(simdata, subj2add)
          }
          table(true_clusters)
          median(rowSums(simdata))
          #################################
          # FILTER OUT ALL ZEROS SUBJECTS #
          #################################
          dim(simdata)
          data <- simdata[,colSums(simdata)>1]
          dim(data)
          true_clusters <- true_clusters[rowSums(data)>0]
          data <- data[rowSums(data)>0,]
          features4components <- as.character(as.vector(which(colSums(simdata)>1)))

          ###########################
          # DIRICHLET PROCESS SETUP #
          ###########################
          n_features <- ncol(data)
          n_patients <- nrow(data)

          files_dirtest <- list.files(dirtest)
          files_dirtest <- files_dirtest[files_dirtest!="results"]
          conf_used <- stringr::str_split_fixed(files_dirtest, "_", 6)[,c(3,4,5,6)]
          conf_used <- conf_used[conf_used[,2]==ncol(simdata),]

          if (typePriorHDPs[idxHDP]=='1') {
              simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="dp1s",,drop=FALSE]
            } else {
              simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="",,drop=FALSE]
          }
          n_people       <- simPrecomputed[,1]
          n_features     <- simPrecomputed[,2]
          n_components   <- simPrecomputed[,3]
          type_dp        <- simPrecomputed[,4]

          pathdir <- paste0(dirtest,"simulated_data_", n_people, "_", n_features, "_", n_components)
          pathdir <- ifelse(type_dp!="", paste0(pathdir,"_",type_dp), pathdir)
          logdir <- file.path(pathdir, mean_freq_rate)

          if (gOpt) {
              outfile_multn_folder <- file.path(logdir, paste0("multinomials_garbage_", choice_opt))
              outfile_fnchd_folder <- file.path(logdir, paste0("fnchds_garbage_", choice_opt))
            } else {
              outfile_multn_folder <- file.path(logdir, paste0("multinomials_", choice_opt))
              outfile_fnchd_folder <- file.path(logdir, paste0("fnchds_", choice_opt))
          }
          if (!file.exists(file.path(outfile_multn_folder, "finalresult_multinomial.rds")) & !file.exists(file.path(outfile_fnchd_folder, "finalresult_fnchd.rds"))) {
            next
          }
          finalresult_multinomial <- readRDS(file.path(outfile_multn_folder, "finalresult_multinomial.rds"))
          finalresult_fnchd <- readRDS(file.path(outfile_fnchd_folder, "finalresult_fnchd.rds"))

          compMultinomial <- do.call("cbind",lapply(finalresult_multinomial$parameters, function(x) x[rownames(x)[order(as.numeric(rownames(x)))],,drop=FALSE]))
          compFNCHD <- do.call("cbind",lapply(finalresult_fnchd$parameters, function(x) x[rownames(x)[order(as.numeric(rownames(x)))],,drop=FALSE]))
          true_clusts <- as.data.frame(t(ws_vec_all[,as.numeric(features4components),drop=FALSE]))
          rownames(true_clusts) <- features4components
          colnames(true_clusts) <- 1:ncol(true_clusts)
          true_clusts_odds <- as.data.frame(t(odds_vec_all[,as.numeric(features4components),drop=FALSE]))
          rownames(true_clusts_odds) <- features4components
          colnames(true_clusts_odds) <- 1:ncol(true_clusts_odds)

          # mean(as.matrix((true_clusts-compMultinomial)**2))
          # mean(as.matrix((t(t(true_clusts_odds)/colSums(true_clusts_odds))-t(t(compFNCHD)/colSums(compFNCHD)))**2))

          # compFNCHD[compFNCHD==Inf] <- 1e6
          # ass_true_multn <- AssociationsCos2(a=t(true_clusts),
          #                                    b=t(compMultinomial),
          #                                    k1="true",
          #                                    k2="multinomial")
          # ass_true_fnchd <- AssociationsCos2(a=t(true_clusts),
          #                                    b=t(compFNCHD),
          #                                    k1="true",
          #                                    k2="fnchd")
          # edges_true_multn <- as_edgelist(ass_true_multn)
          # edges_true_multn <- edges_true_multn[order(as.numeric(stringr::str_split_fixed(edges_true_multn[,1],"_",2)[,1])),,drop=FALSE]
          # edges_true_fnchd <- as_edgelist(ass_true_fnchd)
          # edges_true_fnchd <- edges_true_fnchd[order(as.numeric(stringr::str_split_fixed(edges_true_fnchd[,1],"_",2)[,1])),,drop=FALSE]
          # # stopifnot(all(stringr::str_split_fixed(edges_true_multn[,2],"_", 2)[,1]==stringr::str_split_fixed(edges_true_fnchd[,2],"_", 2)[,1]))
          #
          # # plot(ass_true_multn)
          # # plot(ass_true_fnchd)
          #
          # outDrivers <- data.frame(matrix(NA, ncol=ncol(true_clusts), nrow=nrow(true_clusts)))
          # for (idx in 1:ncol(outDrivers)) {
          #   outDrivers[,idx] <- rownames(true_clusts)[rev(order(true_clusts[,idx]))]
          # }
          # outDrivers[is.na(outDrivers)] <- ""
          # colnames(outDrivers) <- paste("Component",1:ncol(outDrivers),sep="_")
          #
          # # stopifnot(all(stringr::str_split_fixed(edges_true_multn[,2],"_", 2)[,1]==stringr::str_split_fixed(edges_true_fnchd[,2],"_", 2)[,1]))
          #
          # a <- outDrivers[,as.numeric(stringr::str_split_fixed(edges_true_multn[,1],"_",2)[,1]),drop=FALSE]
          # b <- finalresult_multinomial$overview[,as.numeric(stringr::str_split_fixed(edges_true_multn[,2],"_",2)[,1]),drop=FALSE]
          # taus_multinomial <- unlist(lapply(1:ncol(a), function(idx) cor(as.numeric(a[,idx]), as.numeric(b[,idx]), method="kendall")))
          #
          # a <- outDrivers[,as.numeric(stringr::str_split_fixed(edges_true_fnchd[,1],"_",2)[,1]),drop=FALSE]
          # b <- finalresult_fnchd$overview[,as.numeric(stringr::str_split_fixed(edges_true_fnchd[,2],"_",2)[,1]),drop=FALSE]
          # taus_fnchd <- unlist(lapply(1:ncol(a), function(idx) cor(as.numeric(a[,idx]), as.numeric(b[,idx]), method="kendall")))

          # assignment multinomial
          if (choice_opt==9) {
              multinomial_assignments <- MultinomialAssignmentOpp(compMultinomial, data)
              fnchd_assignments <- FNCHDAssignmentOpp(ps=compFNCHD, data=data, logP0s=finalresult_fnchd$logP0s)
              both_assignments <- paste(multinomial_assignments, fnchd_assignments, sep="_")
            } else {
              multinomial_assignments <- MultinomialAssignment(compMultinomial, data)
              fnchd_assignments <- FNCHDAssignment(ps=compFNCHD, data=data, logP0s=finalresult_fnchd$logP0s)
              both_assignments <- paste(multinomial_assignments, fnchd_assignments, sep="_")
          }

          if (gOpt) {
              if (length(unique(true_clusters))==1 & length(unique(fnchd_assignments))==1) {
                ari <- c(1,1,1)
              } else {
                ari <- c(pdfCluster::adj.rand.index(true_clusters[multinomial_assignments!=0], multinomial_assignments[multinomial_assignments!=0]), pdfCluster::adj.rand.index(true_clusters[fnchd_assignments!=0], fnchd_assignments[fnchd_assignments!=0]), pdfCluster::adj.rand.index(true_clusters, both_assignments))
              }
            } else {
              if (length(unique(true_clusters))==1 & length(unique(fnchd_assignments))==1) {
                ari <- c(1,1,1)
              } else {
                ari <- c(pdfCluster::adj.rand.index(true_clusters, multinomial_assignments), pdfCluster::adj.rand.index(true_clusters, fnchd_assignments), pdfCluster::adj.rand.index(true_clusters, both_assignments))
              }
          }
          print(ari)

          which(multinomial_assignments!=fnchd_assignments)
          table(true_clusters, fnchd_assignments)

          # contMultn <- getContingencyTable(grph=ass_true_multn, onevec=true_clusters, twovec=multinomial_assignments, onelab="true", twolab="multinomial")
          # contFNCHD <- getContingencyTable(grph=ass_true_multn, onevec=true_clusters, twovec=fnchd_assignments, onelab="true", twolab="multinomial")
          tmp4res <- data.frame("N"=n_people,
                                "M"=n_features,
                                "K"=n_components,
                                "choice_opt"=choice_opt,
                                "AverageOccurence"=mean_freq_rate,
                                # "kendall_rank"=c(mean(taus_multinomial), mean(taus_fnchd)),
                                "distribution"=c("multinomial", "fnchd", "both"),
                                # "accuracy"=c(sum(diag(contMultn))/sum(contMultn), sum(diag(contFNCHD))/sum(contFNCHD)),
                                # "mcc"=c(custom_mcc(contMultn), custom_mcc(contFNCHD)),
                                "ARI" = ari,
                                "nchanges" = sum(multinomial_assignments!=fnchd_assignments),
                                "priorHDP"=ifelse(type_dp=="", "1/n", "1"),
                                "priorGen"=ifelse(typePriorGens[idxGen]=="1/n","1/n", "1"),
                                "zero_component" = c(sum(multinomial_assignments==0), sum(fnchd_assignments==0), sum(both_assignments==0))
                              )
          res <- rbind(res, tmp4res)
        }
      }
    }
  }
}
res_filt <- res[res$priorGen=="1/n" & res$priorHDP=="1/n",]
res_multn <- res[res$distribution=="multinomial",]
res_fnchd <- res[res$distribution=="fnchd",]
res_both <- res[res$distribution=="both",]
prova <- res_filt %>% group_by(K,AverageOccurence, priorHDP, priorGen) %>% summarize(best=paste(unique(choice_opt[ARI==max(ARI)]), collapse="-"),
                                                                                     best_d=paste(unique(distribution[ARI==max(ARI)]), collapse="-"),
                                                                                     max_ARI = max(ARI),
                                                                                     nchg = nchanges[ARI==max(ARI)]
                                                                                    )
prova[order(prova$nchg),]
rev(sort(table(prova$best)))
rev(sort(table(prova$best_d)))
prova[prova$best_d=="multinomial",]



#%%
for (choice_opt_tmp in c(1,2,3,5,6,7,8)) {
  res4plot <- res[res$choice_opt==choice_opt_tmp,]
  res4plot$Setting <- apply(res4plot[,c("priorGen","priorHDP")], 1, function(x) paste(x,collapse=","))
  res4plot$Setting <- factor(res4plot$Setting, levels=c("1,1", "1/n,1", "1,1/n", "1/n,1/n"))
  res4plot$distribution <- factor(res4plot$distribution, levels=c("multinomial", "fnchd", "both"))
  resplot <- ggplot(res4plot[res4plot$K%in%c(5,10),], aes(x=Setting, y=ARI, fill=distribution)) +
             geom_boxplot() +
             scale_fill_manual(values=c("multinomial"="#ffcc00", "fnchd"="#009900", "both"="#cc0099")) +
             facet_wrap(~K, nrow = 3, scales = "free_y") +
             xlab("Setting: alpha generating, alpha HDP") +
             ggtitle(paste0("Approach index: ", choice_opt_tmp)) +
             ylim(c(0,1))
  # ggsave(file.path("Results", paste0("paperplot_choice",choice_opt_tmp,".png")), print(resplot), width=400,height=280,units="mm",dpi=300)
  # nchplot <- ggplot(res4plot, aes(x=Setting, y=nchanges, fill=distribution)) +
  #            geom_boxplot() +
  #            facet_wrap(~K, nrow = 3, scales = "free_y") +
  #            xlab("Setting: alpha generating, alpha HDP") +
  #            ggtitle(paste0("Approach index: ", choice_opt_tmp))
  # ggsave(file.path("Results", paste0("paperplot_nchplot_choice",choice_opt_tmp,".png")), print(nchplot), width=400,height=280,units="mm",dpi=300)
}


#%%
###########
# 4 PAPER #
###########
choice_opt_tmp <- 7
res4plot <- res[res$choice_opt==choice_opt_tmp,]
res4plot$priorGen <- as.character(res4plot$priorGen)
res4plot$priorHDP <- as.character(res4plot$priorHDP)
res4plot$priorGen[res4plot$priorGen=="1/n"] <- "1/M"
res4plot$priorHDP[res4plot$priorHDP=="1/n"] <- "1/M"
res4plot$priorGen <- factor(res4plot$priorGen, levels=c("1", "1/M"))
res4plot$priorHDP <- factor(res4plot$priorHDP, levels=c("1", "1/M"))
res4plot$Setting <- apply(res4plot[,c("priorGen","priorHDP")], 1, function(x) paste(x,collapse=","))
res4plot$Setting <- factor(res4plot$Setting, levels=c("1,1", "1/M,1", "1,1/M", "1/M,1/M"))
res4plot$distribution <- as.character(res4plot$distribution)
res4plot$distribution[res4plot$distribution!="fnchd" & res4plot$distribution!="both"] <- "Multi"
res4plot$distribution[res4plot$distribution=="fnchd" & res4plot$distribution!="Multi"] <- "MFNCH"
res4plot$distribution <- factor(res4plot$distribution, levels=c("Multi", "MFNCH", "both"))
res4plot$K <- factor(paste0("K==", as.character(res4plot$K)), levels=c("K==1","K==5","K==10"))

res4plot[res4plot$distribution=="Multi" & (res4plot$K %in% c("K==5","K==10")),] %>% group_by(Setting, K) %>% summarize(x=median(ARI))

resplot <- ggplot(res4plot[res4plot$distribution=="Multi" & (res4plot$K %in% c("K==5","K==10")),], aes(x=Setting, y=ARI)) +
           geom_boxplot(fill="#009900") +
           # scale_fill_manual(values=c("Multi"="#009900")) +
           facet_wrap(~K, nrow = 3, scales = "free_y", labeller=label_parsed) +
           xlab(expression(paste(alpha[sim],",",alpha[HDMM]))) +
           theme(text = element_text(size = 20))+
           # ggtitle(paste0("Approach index: ", choice_opt_tmp)) +
           ylim(c(0,1))
ggsave(file.path("PAPER/", paste0("paperplot_choice",choice_opt_tmp,".eps")), print(resplot), width=400,height=280,units="mm",dpi=600)


res4plotpaper2 <- res4plot[res4plot$distribution!="Multi",]
res4plotpaper2$Difference <- res4plot[res4plot$distribution!="Multi","ARI"]-res4plot[res4plot$distribution=="Multi","ARI"]

res4plotpaper2 %>% group_by(Setting, K) %>% summarize(x=max(Difference))


resplot2 <- ggplot(res4plotpaper2[res4plotpaper2$K%in%c("K==5","K==10"),], aes(x=Setting, y=Difference)) +
            geom_boxplot(fill="#3399ff") +
            facet_wrap(~K, nrow = 2, scales = "free_y", labeller=label_parsed) +
            xlab(expression(paste(alpha[sim],",",alpha[HDMM]))) +
            ylab("MFNCH impact on ARI w.r.t. Multi") +
            theme(text = element_text(size = 20))
            # ggtitle(paste0("Approach index: ", choice_opt_tmp))
ggsave(file.path("PAPER/", paste0("resplot2_",choice_opt_tmp,".eps")), print(resplot2), width=400,height=280,units="mm",dpi=600)




#%%
library(ggplot2)
library(dplyr)
options(repr.plot.width=16, repr.plot.height=10)

summary_conv <- as.data.frame(data.table::fread("summary_conv.csv"))


df <- res[res$choice_opt==6,]
df$K <- factor(df$K, levels=c("1","5","10"))
df$priorGen <- factor(df$priorGen, levels=c("1","1/n"))
df$distribution <- factor(df$distribution, levels=c("multinomial","fnchd"))
df$AverageOccurence <- factor(df$AverageOccurence, levels=1:10)
# df$kendall_rank <- as.numeric(df$kendall_rank)

df$Condition <- NA
for (ridx in 1:nrow(df)) {
  condridx <- summary_conv$N==df[ridx,"N"]
  condridx <- condridx & summary_conv$M==df[ridx,"M"]
  condridx <- condridx & summary_conv$K==df[ridx,"K"]
  condridx <- condridx & summary_conv$AverageOccurence==df[ridx,"AverageOccurence"]
  condridx <- condridx & summary_conv$priorHDP==df[ridx,"priorHDP"]
  condridx <- condridx & summary_conv$priorGen==df[ridx,"priorGen"]
  df[ridx,"Condition"] <- median(summary_conv[ridx,"DominanceFrequency"])>0.5
}
stopifnot(!anyNA(df$Condition))

which.max(df[df$distribution == "multinomial","ARI"]-df[df$distribution != "multinomial","ARI"])
df[df$distribution == "multinomial",][72,]
df[df$distribution != "multinomial",][72,]


dir.create(file.path(dirtest, "results"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
# kenplot <- ggplot(df, aes(x=AverageOccurence, y=kendall_rank, fill=distribution)) +
#            geom_boxplot() +
#            facet_wrap(~K, nrow = 3, scales = "free_y")
# # ggsave(file.path(dirtest, "results", paste0("kenplot_",basename(dirtest),".png")), print(kenplot), width=400,height=280,units="mm",dpi=300)
# accplot <- ggplot(df, aes(x=AverageOccurence, y=accuracy, fill=distribution)) +
#            geom_boxplot() +
#            facet_wrap(~K, nrow = 3, scales = "free_y")
#[df$priorHDP=="1",]
# mccplot <- ggplot(df_prova2[df_prova2$priorHDP=="1" & df_prova2$priorGen=="1",], aes(x=AverageOccurence, y=mcc, fill=distribution)) +
#            geom_boxplot() +
#            facet_wrap(~K, nrow = 3, scales = "free_y")
#
ariplot <- ggplot(df[df$priorHDP=="1/n" & df$priorGen=="1/n",], aes(x=AverageOccurence, y=ARI, fill=distribution)) +
           geom_boxplot() +
           facet_wrap(~K, nrow = 3, scales = "free_y")
# ggsave(file.path(dirtest, "results", paste0("ariplot_",basename(dirtest),".png")), print(ariplot), width=400,height=280,units="mm",dpi=300)


#%%
#############
# REAL DATA #
#############
realdata <- as.data.frame(data.table::fread("realdata.csv", head=TRUE))
rownames(realdata) <- realdata$Gene
realdata["Gene"] <- NULL
for (rd in colnames(realdata)) {
  realdata[,rd] <- as.numeric(realdata[,rd])
}
realdata <- realdata[rowSums(realdata)>1,]

size_of_feats <- rowSums(realdata)
mstmp <- as.vector(size_of_feats)
RealComponents <- data.frame(matrix(0, nrow=length(size_of_feats),ncol=ncol(realdata)))
RealComponentsMultn <- data.frame(matrix(0, nrow=length(size_of_feats),ncol=ncol(realdata)))
for (idx in 1:ncol(realdata)) {
  mus_comptmp_norm <- realdata[,idx]
  mus_comptmp_norm <- ifelse(mus_comptmp_norm>1e-10, mus_comptmp_norm-1e-10, 0)
  mus_comptmp_norm <- ifelse(mus_comptmp_norm==0, 1e-20, mus_comptmp_norm)
  multn_params <- mus_comptmp_norm/sum(mus_comptmp_norm)
  n_comptmp <- sum(mus_comptmp_norm)
  fnchd_params <- suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=mus_comptmp_norm, m=mstmp, n=n_comptmp, precision=0.1))
  RealComponents[,idx] <- fnchd_params
  RealComponentsMultn[,idx] <- multn_params
}
rownames(RealComponents) <- rownames(realdata)
rownames(RealComponentsMultn) <- rownames(realdata)
realcomponents_overview <- driver_list(cs=RealComponents, features=rownames(realdata))
realcomponents_overviewmultn <- driver_list(cs=RealComponentsMultn, features=rownames(realdata))
head(realcomponents_overviewmultn[,2:11])

realcomponents_details <- lapply(1:ncol(RealComponents), function(x) {
  xc <- RealComponents[, x, drop=FALSE]
  # xc <- xc[rowSums(xc)>0,, drop=FALSE]
  xc <- xc[rev(order(xc[,1])),,drop=FALSE]
  return(xc)
})
head(realcomponents_details[[2]])

knitr::kable(head(realcomponents_overview[,2:11]),"latex")


realdata










#%%


















































































#%%
for (idxGen in 1:2) {
  message(idxGen)
  for (idxHDP in 1:2) {
    message(idxHDP)
    for (n_components in c(1,5,10)) {
      for (mean_freq_rate in 1:10) {
        # mean_freq_rate <- 10
        # n_components <- 5
        print(mean_freq_rate)
        n_people       <- 1000
        n_features     <- 50
        n_components   <- as.numeric(n_components)
        set.seed(1)
        #######################################################################
        # SIMULATE MIXTURE OF FISHER NON-CENTRAL HYPERGEOMETRIC DISTRIBUTIONS #
        #######################################################################
        freq_rate_all <- rpois(n_people, lambda=mean_freq_rate)
        if (typePriorGens[idxGen]=="1/n") {
            ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                       alpha=rep(1/n_features,n_features))
            dirtest <- "alpha_1_over_n_features/"
          } else {
            ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                                       alpha=rep(1.,n_features))
            dirtest <- "alpha_1s/"
        }

        odds_vec_all  <- ws_vec_all/apply(ws_vec_all, 1, median)
        # t(apply(odds_vec_all, 1, function(x) head(rev(order(x))) ))
        composition   <- LaplacesDemon::rdirichlet(n=1, alpha=rep(1,n_components))
        true_clusters <- rep(0, n_people)
        true_clusters[freq_rate_all!=0] <- sample(seq(1,n_components),
                                                  n_people-sum(freq_rate_all==0),
                                                  replace=TRUE,
                                                  prob=composition)
        simdata <- data.frame()
        for (idx in 1:length(freq_rate_all)) {
          if (freq_rate_all[idx]==0) {
            stopifnot(true_clusters[idx]==0)
            subj2add <- rep(0, n_features)
            } else {
            subj2add <- BiasedUrn::rMFNCHypergeo(nran = 1,
                                                 m = rep(1, n_features),
                                                 n = freq_rate_all[idx],
                                                 odds = odds_vec_all[true_clusters[idx],])
          }
          simdata <- rbind(simdata, subj2add)
        }
        table(true_clusters)
        median(rowSums(simdata))
        #################################
        # FILTER OUT ALL ZEROS SUBJECTS #
        #################################
        estimated_ms <- do.call("rbind", lapply(1:max(true_clusters), function(x) colSums(simdata[true_clusters==x,])))
        estimated_nums <- do.call("rbind", lapply(1:max(true_clusters), function(x) sum(true_clusters==x)))


        # multinomials_estimates <- estimated_ms/rowSums(estimated_ms)
        # fnchd_estimates <- do.call("rbind",lapply(1:max(true_clusters), function(x) suppressWarnings(BiasedUrn::oddsMFNCHypergeo(mu=estimated_ms[x,], m=rep(estimated_nums[x,], ncol(simdata)), n=sum(estimated_ms[x,]), precision=0.1))))

        multinomials_estimates <- estimated_ms/sum(estimated_ms)
        fnchd_estimates  <- multinomials_estimates/median(multinomials_estimates[multinomials_estimates>0])

        # xtots <- unique(Reduce("+", raw_comp2feature[[x[[1]]]]))
        # ntot <- ncol(mergedSizes)
        # nobs <- mean(rowSums(mergedSizes>0))
        # xexps <- xtots*nobs/ntot
        # mus_comptmp <- mus_comptmp/xexps
        # multn_params <- mus_comptmp/sum(mus_comptmp)
        # fnchd_params  <- multn_params/median(multn_params[multn_params>0])

        compFNCHD <- as.data.frame(t(fnchd_estimates))
        ExclusiveDefiningLesions <- (1:ncol(simdata))[which(rowSums(compFNCHD==Inf)>0)]
        stopifnot(all(unique(as.vector(as.matrix(compFNCHD[ExclusiveDefiningLesions,]))) %in% c(0,Inf)))
        logP0s <- getNormalizingFactors(compFNCHD[!rownames(compFNCHD) %in% ExclusiveDefiningLesions,], max(rowSums(simdata,na.rm=TRUE)))
        logP0s <- as.data.frame(logP0s)
        colnames(logP0s) <- NULL


        multinomial_assignments <- MultinomialAssignment(t(multinomials_estimates), simdata)
        fnchd_assignments <- FNCHDAssignment(ps=t(fnchd_estimates), data=simdata, logP0s=logP0s)
        if (length(unique(true_clusters))==1 & length(unique(fnchd_assignments))==1) {
            ari <- c(1,1)
          } else {
            ari <- c(pdfCluster::adj.rand.index(true_clusters, multinomial_assignments), pdfCluster::adj.rand.index(true_clusters, fnchd_assignments))
        }
        print(ari)






        # files_dirtest <- list.files(dirtest)
        # files_dirtest <- files_dirtest[files_dirtest!="results"]
        # conf_used <- stringr::str_split_fixed(files_dirtest, "_", 6)[,c(3,4,5,6)]
        # conf_used <- conf_used[conf_used[,2]==ncol(simdata),]
        #
        # if (typePriorHDPs[idxHDP]=='1') {
        #     simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="dp1s",,drop=FALSE]
        #   } else {
        #     simPrecomputed <- conf_used[conf_used[,3]==n_components & conf_used[,4]=="",,drop=FALSE]
        # }
        # n_people       <- simPrecomputed[,1]
        # n_features     <- simPrecomputed[,2]
        # n_components   <- simPrecomputed[,3]
        # type_dp        <- simPrecomputed[,4]
        #
        # pathdir <- paste0(dirtest,"simulated_data_", n_people, "_", n_features, "_", n_components)
        # pathdir <- ifelse(type_dp!="", paste0(pathdir,"_",type_dp), pathdir)
        # logdir <- file.path(pathdir, mean_freq_rate)
        # dir.create(file.path(logdir, "simulated"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
        # saveRDS(multinomials_estimates, file.path(logdir, "simulated", "est_mltns.rds"), )
        # saveRDS(fnchd_estimates, file.path(logdir, "simulated", "est_fnchds.rds"))
        # saveRDS(logP0s, file.path(logdir, "simulated", "est_fnchds_logP0s.rds"))
      }
    }
  }
}


#%%
###########
# THE END #
###########
