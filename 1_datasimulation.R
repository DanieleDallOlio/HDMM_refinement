library(hdp)
# mean_freq_rate <- 3
for (mean_freq_rate in 1:10) {
  message(mean_freq_rate)
  n_people       <- 1000
  n_features     <- 50
  n_components   <- 10
  # pathdir <- paste0("alpha_1_over_n_features/simulated_data_", n_people, "_", n_features, "_", n_components)
  # pathdir <- paste0("alpha_1s/simulated_data_", n_people, "_", n_features, "_", n_components)
  # pathdir <- paste0("alpha_1s/simulated_data_", n_people, "_", n_features, "_", n_components, "_dp1s")
  pathdir <- paste0("alpha_1_over_n_features/simulated_data_", n_people, "_", n_features, "_", n_components, "_dp1s")
  set.seed(1)
  #######################################################################
  # SIMULATE MIXTURE OF FISHER NON-CENTRAL HYPERGEOMETRIC DISTRIBUTIONS #
  #######################################################################
  freq_rate_all <- rpois(n_people, lambda=mean_freq_rate)
  ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
                                             alpha=rep(1/n_features,n_features))
  # ws_vec_all    <- LaplacesDemon::rdirichlet(n=n_components,
  #                                            alpha=rep(1.,n_features))
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
  dim(data)


  #################
  # LOGGING SETUP #
  #################
  logdir <- file.path(pathdir, mean_freq_rate)
  dir.create(file.path(logdir, "chains"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
  dp_log <- file.path(logdir, paste0("hdp_log_", mean_freq_rate, ".txt"))
  printToLog <- function(s) {
    cat(s, file = dp_log, append = T, sep = "\n")
  }



  ###########################
  # DIRICHLET PROCESS SETUP #
  ###########################
  n_features <- ncol(data)
  n_patients <- nrow(data)

  burnin      <- 5000
  postsamples <- 1e4
  spacebw     <- 20
  cpsamples   <- 10
  dom_thres   <- 0.6
  # hhOpt <- rep(1/n_features, n_features)
  hhOpt <- rep(1., n_features)


  #################
  # HDP FIRST RUN #
  #################
  seeds_activated <- sample(1:1e3,10)
  seeds_posterior <- sample(1:1e2,10)
  firstRunHDP <- function(outfile_chain_fn, seed_k_act, seed_k_pos) {
    hdp <- hdp_init(ppindex=c(0, rep(1,n_patients)),
                    cpindex=c(1, rep(1,n_patients)),
                    hh=hhOpt,
                    alphaa=c(1),
                    alphab=c(1))
    hdp <- hdp_setdata(hdp=hdp,
                       dpindex=2:numdp(hdp),
                       data=data)
    hdp_activated <- dp_activate(hdp=hdp,
                                 dpindex=1:numdp(hdp),
                                 initcc=25,
                                 seed=seed_k_act
                                )
    chlist <- hdp_posterior(hdp=hdp_activated,
                            burnin=burnin,
                            n=postsamples,
                            space=spacebw,
                            cpiter=cpsamples,
                            seed=seed_k_pos
                            )
    saveRDS(chlist, outfile_chain_fn)
  }


  ###############################
  # DIRICHLET PROCESS EXECUTION #
  ###############################
  # k<-1
  for (k in 1:10) {
    seed_act <- seeds_activated[[k]]
    seed_post <- seeds_posterior[[k]]
    outfile_chain_fn <- file.path(logdir, "chains", paste0("HDMM_it_", k, "_seedact_", seed_act, "_seedpost_", seed_post, ".rds"))
    if (!file.exists(outfile_chain_fn)) {
      firstRunHDP(outfile_chain_fn, seed_k_act=seed_act, seed_k_pos=seed_post)
      # outfile_chain <- readRDS(outfile_chain_fn)
      # saveRDS(outfile_chain, outfile_chain_fn)
      rm(outfile_chain)
      gc()
    }
  }
}


#%%
###########
# THE END #
###########
