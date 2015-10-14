#' Item Focused Trees for the Identification of Items in Differential Item Functioning 
#' 
#' @description
#' A function to estimate Item focussed trees for simultaneous selection of items and variables 
#' that induce DIF (Differential Item Functioning) in Item Response models. 
#' The method of item focussed recursive partitioning is described in Tutz and Berger (2015).
#' 
#' @param Y Matrix or Data.frame of binary 0/1 response (rows correspond to persons, columns correspond to items)
#' @param X Data.frame of (not scaled) covariates (rows correspond to persons, columns correspond to covariates)
#' @param alpha global significance level for the permutation tests
#' @param nperm number of permutations for the permutation tests
#' @param trace If true, information about the estimation progress is printed
#' @param penalize If true, a small ridge penalty is added to ensure existence of model parameters 
#' @param x,object Object of class \code{"DIFtree"}
#' @param ... further arguments passed to or from other methods
#' 
#' @details The method assumes a Rasch model where DIF can appear in some of the items. 
#' Items with DIF are gradually identified by recursive partitioning. 
#' In each iteration of the algorithm one item, variable and threshold ist chosen for splitting. 
#' Significance of the split is verified by permutation tests. 
#' The result of the permutation test can strongly depend on the number of permutations \code{nperm}.
#' 
#' In the case of pure terminal nodes estimates of the model do not exist. If \code{penalize=TRUE} 
#' a small ridge penalty is added during estimation to ensure existence of all parameters. 
#'
#' 
#' @return object of class \code{"DIFtree"}. 
#' An object of class \code{"DIFtree"} is a list containing the following components:
#' 
#' \item{splits}{Matrix with detailed information about all executed splits during the estimation process}
#' \item{theta_hat}{Estimated person abilities}
#' \item{beta_hat_nodif}{Estimated item difficulties for items without DIF} 
#' \item{beta_hat_dif}{List of estimated item difficulties for items with DIF. 
#' Each element contains diffculties for one item.}
#' \item{pvalues}{P-values of each permutation test during the estimation process}
#' \item{Y}{Response matrix used in the estimation process}
#' \item{X}{Model matrix used in the estimation process}
#' \item{persons}{Number of persons} 
#' \item{items}{Number of items} 
#' 
#' @author Moritz Berger <moritz.berger@@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Department of Statistics, LMU Munich
#' 
#' @seealso \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}
#' 
#' @examples 
#' 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,alpha=0.05,nperm=100,trace=TRUE)
#'  
#' print(mod)
#' summary(mod)
#' } 
#'
#
#' @exportClass DIFtree
#' @export
#' @importFrom penalized penalized 
#' @importFrom stats binomial coef deviance formula glm predict quantile 

DIFtree <-
function(Y, X, alpha, nperm, trace=FALSE, penalize=FALSE){
  
  # check input 
  if(!all(as.vector(t(Y)) %in% c(0,1))){
    stop("Y must be a binary 0/1 matrix")
  }
  if(!is.data.frame(X)){
    stop("X must be of class 'data.frame'")
  }
  if(nrow(X)!=nrow(Y)){
    stop("Dimensions of X and Y don't match")
  }
  if(any(sapply(1:ncol(X),function(j) class(X[,j]))=="character")){
    stop("variable of class 'character' is not useful")
  }
  if(any(sapply(1:ncol(X),function(j) class(X[,j]))=="logical")){
    stop("variable of class 'logical' is not useful")
  }
  if(any(sapply(1:ncol(X), function(j) {
    is.numeric(X[,j]) && (mean(X[,j])<10e-6 | var(X[,j])==1)}))){
    stop("Don't use scaled covariates")
  }
  
  npersons <- nrow(Y)           # Anzahl an Beobachtungen 
  nitems   <- ncol(Y)           # Anzahl an Items
  nvar     <- ncol(X)           # Anzahl an Variablen 
  y        <- as.vector(t(Y))   # Erzeuge aus Zielmatrix einen Vektor 
  
  # modify DM_kov
  DM_kov <- X 
  for(i in 1:nvar){
    if(is.factor(DM_kov[,i])){
      if(is.ordered(DM_kov[,i])){
        DM_kov[,i] <- as.numeric(DM_kov[,i])
      } else{
        if(nlevels(DM_kov[,i])==2){
          DM_kov[,i] <- as.numeric(DM_kov[,i])-1
        } else{
          labs <- levels(DM_kov[,i])
          DM_kov[,i] <- as.numeric(DM_kov[,i])
          for(j in 2:length(labs)){
            DM_kov[,labs[j]] <- ifelse(DM_kov[,i]==j,1,0)
          }
          DM_kov <- DM_kov[,-i]
        }
      }
    }
  }
  
  # Designmatrix des Rasch-Modells 
  pp_design <- diag(npersons)                   # Personenparameter, Person P ist Referenz 
  pp_design <- pp_design[rep(1:nrow(pp_design),each=nitems),]
  pp_design <- pp_design[,-npersons]
  
  ip_design <- -1*diag(nitems)                  # Itemparameter 
  ip_design <- ip_design[rep(1:nrow(ip_design),times=npersons),]
  
  dm_rasch <- cbind(pp_design,ip_design)               # Designmatrix fuer Rasch-Modell zusammenfuegen 
  
  names_rasch    <- c(paste("theta",1:(npersons-1),sep=""),paste("beta",1:nitems,sep=""))
  colnames(dm_rasch) <- names_rasch      # Variablennamen vergeben 
  
  # Bestimme geordnete Werte
  ordered_values <- lapply(1:nvar, function(j){
    if(!all((DM_kov[,j] - round(DM_kov[,j])) == 0)){
      quantile(DM_kov[,j],seq(0.05,1,by=0.05))
    } else{
      unique(sort(DM_kov[,j]))
    }
  })
  names(ordered_values) <- names(DM_kov)
  
  n_levels <- sapply(1:nvar, function(j) length(ordered_values[[j]]))
  n_s      <- n_levels-1
  
  designlists <- function(DM_kov){
  
    # Generiere Designmatrizen
    design_help_upper <- lapply(1:nvar, function(j){
      sapply(1:(n_levels[j]-1),function(k) { ifelse(DM_kov[,j] > ordered_values[[j]][k],1,0)})
    })
    # design_help_upper <- do.call(cbind,design_help_upper)
    
    design_help_lower <- lapply(1:nvar, function(j){
      sapply(1:(n_levels[j]-1),function(k) { ifelse(DM_kov[,j] > ordered_values[[j]][k],0,1)})
    })
    # design_help_lower <- do.call(cbind,design_help_lower)
    
    v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1)) 
    w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))
    
    design_upper <- lapply(1:nitems, function(j) {
      lapply(1:nvar, function(var){
        design_tree <- matrix(0,nrow=nitems*npersons,ncol=ncol(design_help_upper[[var]]))
        rows        <- seq(j,(nitems*npersons),by=nitems)
        design_tree[rows,] <- design_help_upper[[var]]
        z <- rep(paste0("_u",j),ncol(design_help_upper[[var]]))
        colnames(design_tree) <- paste0(w[[var]],v[[var]],z)
        return(design_tree)
      })
    })
    
    v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1)) 
    w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))
    
    design_lower <- lapply(1:nitems, function(j) {
      lapply(1:nvar, function(var){
        design_tree <- matrix(0,nrow=nitems*npersons,ncol=ncol(design_help_lower[[var]]))
        rows        <- seq(j,(nitems*npersons),by=nitems)
        design_tree[rows,] <- design_help_lower[[var]]
        z <- rep(paste0("_l",j),ncol(design_help_lower[[var]]))
        colnames(design_tree) <- paste0(w[[var]],v[[var]],z)
        return(design_tree)
      })
    }) 
    return(list(design_upper,design_lower))
  }
  
  #########################################################################################
  
  mod_potential <- list()
  devs          <- c()
  crit          <- c()
  splits        <- c()
  pvalues       <- c() 
  ip            <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  numbers       <- list()
  count <- 1
  
  pp     <- paste("theta",1:(npersons-1),sep="")
  help_p <- paste0(pp,collapse="+")
  
  numbers[[1]]     <- lapply(1:nitems,function(j) 1)
  which_obs[[1]]   <- lapply(1:nitems,function(j) matrix(1:npersons,nrow=1))
  splits_evtl[[1]] <- lapply(1:nitems,function(j) lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1)))
  vars_evtl[[1]]   <- lapply(1:nitems,function(j) nvar)
  ip[[1]]          <- lapply(1:nitems,function(j) paste0("beta",j))
  help0            <- formula(paste("y~",help_p,"+",paste0(unlist(ip[[1]]),collapse="+"),"-1"))
  dat0             <- data.frame(y,dm_rasch)
  mod0             <- mod_potential[[1]] <- glm(help0,family=binomial(link="logit"),data=dat0)
  start            <- predict(mod0)
  
  # function to compute all models in one knot
  
  allmodels <- function(i,var,kn,design_lower,design_upper){
    
    deviances <- rep(0,n_s[var])
    help_kn <- ip[[count]][[i]][kn]
    help1   <- paste0(unlist(ip[[count]])[-which(unlist(ip[[count]])==help_kn)],collapse="+")
    splits_aktuell <- splits_evtl[[count]][[i]][[var]][kn,] 
    splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]
    
    param_knot <- coef(mod0)[help_kn]
    obs_aktuell <- which_obs[[count]][[i]][kn,]
    obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
    rows_obs   <- seq(i,(nitems*npersons),by=nitems)[obs_aktuell]
    start[rows_obs] <- start[rows_obs]+param_knot
        
    if(length(splits_aktuell)>0){
        
      for(j in splits_aktuell){
        dat   <- data.frame(dat0,design_lower[[i]][[var]][,j,drop=FALSE],design_upper[[i]][[var]][,j,drop=FALSE])
        help2 <- paste(ip[[count]][[i]][kn],c(colnames(design_lower[[i]][[var]])[j],colnames(design_upper[[i]][[var]])[j]),sep=":")
        help3 <- paste(help2,collapse="+")
        help4 <- formula(paste("y~",help3,"-1"))
        suppressWarnings(
        mod <- glm(help4,family=binomial(link="logit"),data=dat,offset=start)
        )
        deviances[j] <-  deviance(mod0)-deviance(mod)
      }
    }
    return(deviances)
  }  
  
  design_upper <- designlists(DM_kov)[[1]]
  design_lower <- designlists(DM_kov)[[2]]
  sig   <- TRUE
  anysplit <- TRUE
  
  while(sig & anysplit){
    
    # compute all models
    dv <- lapply(1:nvar,function(var) {
            lapply(1:nitems,function(i) {
              n_knots   <- length(ip[[count]][[i]])
              deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
              for(kn in 1:n_knots){
                deviances[,kn] <- allmodels(i,var,kn,design_lower,design_upper)
              }
              return(deviances)
            })
          })
                      
    # select optimum
    variable <- which.max(lapply(1:nvar,function(j) max(unlist(dv[[j]]))))
    item     <- which.max(lapply(1:nitems, function(j) max(dv[[variable]][[j]])))
    split    <- as.numeric(which(dv[[variable]][[item]]==max(dv[[variable]][[item]]),arr.ind=TRUE)[,1])
    knoten   <- as.numeric(which(dv[[variable]][[item]]==max(dv[[variable]][[item]]),arr.ind=TRUE)[,2])
    if(length(split)>1){
      split  <- split[1]
      knoten <- knoten[1]
      warning(paste("Maximum in iteration ",count," not uniquely defined"))
    }
    ip_old   <- ip[[count]][[item]][knoten]
    level    <- length(strsplit(ip_old,":")[[1]])
    number   <- numbers[[count]][[item]][knoten]
    left     <- max(numbers[[count]][[item]])+1
    right    <- max(numbers[[count]][[item]])+2
    
    # compute permutation test 
    dev <- rep(NA,nperm)

    for(perm in 1:nperm){
      dv_perm <- rep(0,n_s[variable])
      obs_aktuell <- which_obs[[count]][[item]][knoten,]
      obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
      DM_kov_perm <- DM_kov
      DM_kov_perm[obs_aktuell,variable] <- sample(DM_kov_perm[obs_aktuell,variable],length(obs_aktuell))
      design_upper_perm      <- designlists(DM_kov_perm)[[1]]
      design_lower_perm      <- designlists(DM_kov_perm)[[2]]
      dv_perm <- allmodels(item,variable,knoten,design_lower_perm,design_upper_perm)
      dev[perm] <- max(dv_perm)
      if(trace){
        cat(".")
      }
    }
    
    # test decision 
    crit_val <- quantile(dev,1-(alpha/vars_evtl[[count]][[item]][knoten]))
    proof    <- max(dv[[variable]][[item]]) > crit_val
    devs[count]   <- max(dv[[variable]][[item]])
    crit[count]   <- crit_val
    pvalues[count] <- length(which(dev>max(dv[[variable]][[item]])))/nperm
    
    if(proof){
      
      # get new formula 
      help_kn2 <- ip[[count]][[item]][knoten]
      help5 <- paste0(unlist(ip[[count]])[-which(unlist(ip[[count]])==help_kn2)],collapse="+")
      help6 <- paste(ip[[count]][[item]][knoten],c(colnames(design_lower[[item]][[variable]])[split],colnames(design_upper[[item]][[variable]])[split]),sep=":")
      help7 <- paste(help6,collapse="+")
      help8 <- formula(paste("y~",help_p,"+",help5,"+",help7,"-1"))
      
      ######################
      if(level>1){
        help_kn4 <- lu(c(),1,level-1,c())
        help_kn5 <- unlist(strsplit(help_kn2,""))
        help_kn6 <- paste0(help_kn5[which(help_kn5=="_")+1],collapse="")
        knoten2  <- which(help_kn4==help_kn6)
      } else{
        knoten2 <- knoten
      }
      ######################
      
      splits <- rbind(splits,c(variable,item,split,level,knoten2,number,left,right))   

      
      # fit new model 
      dat <- dat0 <- data.frame(dat0,design_lower[[item]][[variable]][,split,drop=FALSE],design_upper[[item]][[variable]][,split,drop=FALSE])
      suppressWarnings(
        mod0  <- mod_potential[[count+1]] <- glm(help8,family=binomial(link="logit"),data=dat,etastart=start)
      )
      start <- predict(mod0)
    
      # generiere neue itemparameter 
      ip[[count+1]]                             <- ip[[count]]
      ip[[count+1]][[item]]                     <- rep("",length(ip[[count]][[item]])+1)
      ip[[count+1]][[item]][c(knoten,knoten+1)] <- help6
      ip[[count+1]][[item]][-c(knoten,knoten+1)]<- ip[[count]][[item]][-knoten]
      
      # passe splits_evtl an
      n_knots                                                       <- length(ip[[count+1]][[item]])
      splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
      for(var in 1:nvar){
        splits_evtl[[count+1]][[item]][[var]]                       <- matrix(0,nrow=n_knots,ncol=n_s[var])
        splits_evtl[[count+1]][[item]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[item]][[var]][knoten,],2),nrow=2,byrow=T)
        splits_evtl[[count+1]][[item]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[item]][[var]][-knoten,]
      }
      splits_evtl[[count+1]][[item]][[variable]][knoten,splits_evtl[[count+1]][[item]][[variable]][knoten,]>=split] <- NA 
      splits_evtl[[count+1]][[item]][[variable]][(knoten+1),splits_evtl[[count+1]][[item]][[variable]][(knoten+1),]<=split] <- NA
      
      # any split? 
      anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))
      
      # passe vars_evtl an 
      vars_evtl[[count+1]]                             <- vars_evtl[[count]]
      vars_evtl[[count+1]][[item]]                     <- rep(0,n_knots)
      vars_evtl[[count+1]][[item]][c(knoten,knoten+1)] <- rep(vars_evtl[[count]][[item]][knoten],2)
      vars_evtl[[count+1]][[item]][-c(knoten,knoten+1)]<- vars_evtl[[count]][[item]][-knoten]
      
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[variable]][knoten,])))==0){ 
        vars_evtl[[count+1]][[item]][knoten] <- vars_evtl[[count+1]][[item]][knoten]-1 
      }
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[variable]][knoten+1,])))==0){ 
        vars_evtl[[count+1]][[item]][knoten+1] <- vars_evtl[[count+1]][[item]][knoten+1]-1 
      }
      
      # passe which_obs an 
      which_obs[[count+1]]                               <- which_obs[[count]]
      which_obs[[count+1]][[item]]                       <- matrix(0,nrow=n_knots,ncol=npersons)
      which_obs[[count+1]][[item]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][[item]][knoten,],2),nrow=2,byrow=T)
      which_obs[[count+1]][[item]][-c(knoten,knoten+1),] <- which_obs[[count]][[item]][-knoten,]
      thresh <- ordered_values[[variable]][1:n_s[variable]][split]
      which_obs[[count+1]][[item]][knoten,DM_kov[,variable]>thresh] <- NA
      which_obs[[count+1]][[item]][(knoten+1),DM_kov[,variable]<=thresh] <- NA
      
      # passe numbers an 
      numbers[[count+1]]                              <- numbers[[count]]
      numbers[[count+1]][[item]]                      <- numeric(length=n_knots)
      numbers[[count+1]][[item]][c(knoten,knoten+1)]  <- c(left,right)
      numbers[[count+1]][[item]][-c(knoten,knoten+1)] <- numbers[[count]][[item]][-knoten] 
      
      # trace
      if(trace){
        cat(paste0("\n Split"," ",count,";"," ","Item"," ",item,"\n"))
      }
      
      # erhoehe counter
      count <- count+1 

    } else{
      sig <- FALSE
      if(penalize){
        if(count>1){
          #library("penalized")
          help9 <- formula(paste("~",0,"+",help_p))
          help10 <- formula(paste("~",help5,"+",help7))
          mod_potential[[count]] <- penalized(y,penalized=help10,unpenalized=help9,lambda2=1e-3,data=dat0,trace=FALSE)
        } else{
          help9  <- formula(paste("~",0,"+",help_p))
          help10 <- formula(paste("~",paste0(unlist(ip[[1]]),collapse="+")))
          mod_potential[[count]] <- penalized(y,penalized=help10,unpenalized=help9,lambda2=1e-3,data=dat0,trace=FALSE)
        }
      }
    }
  }
    
  ################################################################################### 

  mod_opt     <- mod_potential[[count]]
  ip_opt      <- ip[[count]]
  theta_hat   <- c(coef(mod_opt)[1:(npersons-1)],0)
  beta_hat    <- coef(mod_opt)[npersons:length(coef(mod_opt))]
  
  if(count>1){

    dif_items   <- unique(splits[,2])
    nodif_items <- c(1:nitems)[-dif_items]
  
    beta_hat_nodif <- sapply(nodif_items,function(j) beta_hat[ip_opt[[j]]])
    beta_hat_dif   <- lapply(dif_items, function(j)  beta_hat[ip_opt[[j]]])
    names(beta_hat_dif) <- dif_items
  
    help9 <- cumsum(c(0,(n_levels-1)))
    colnames(splits) <- c("var","item","split","level","node","number","left","right")
    splits <- data.frame(cbind(splits[,1:5,drop=FALSE],"variable"=rep(NA,nrow(splits)),"threshold"=rep(NA,nrow(splits)),splits[,6:8,drop=FALSE]))
    for(i in 1:nrow(splits)){
      splits[i,6] <- colnames(DM_kov)[splits[i,1]]
      v2 <- lapply(1:nvar,function(j) ordered_values[[j]][-length(ordered_values[[j]])])
      splits[i,7] <- v2[[splits[i,1]]][splits[i,3]]
    }
    splits <- splits[,-1]
  } else{
    
    beta_hat_nodif <- beta_hat
    beta_hat_dif   <- c()
  
  }
  
  to_return <- list("splits"=splits,
                    "theta_hat"=theta_hat,
                    "beta_hat_nodif"=beta_hat_nodif,
                    "beta_hat_dif"=beta_hat_dif,
                    "pvalues"=pvalues,
                    "Y"=Y,
                    "X"=DM_kov,
                    "persons"=npersons,
                    "items"=nitems,
                    "call"=match.call())
  class(to_return) <- "DIFtree"
  return(to_return)
}

#' @rdname DIFtree 
#' @method print DIFtree
#' @export
#' 
print.DIFtree <-
  function(x,...){
    npersons  <- x$persons
    nitems    <- x$items
    dif_items <- unique(x$splits[,"item"])
    splits    <- x$splits[,c("item","variable","threshold")]
    cat("\n")
    cat("Number of persons:", npersons, "\n")
    cat("Number of items:", nitems, "\n")
    if(!is.null(dif_items)){
      cat("DIF items:", dif_items, "\n")
    } else{
      cat("DIF items: no DIF item\n")
    }
    cat("\n")
    cat("Overview of executed splits:\n")
    if(!is.null(dif_items)){
      cat("\n")
      print(splits)
    } else{
      cat("no split performed")
    }
    invisible(x)
  } 

#' @rdname DIFtree 
#' @method summary DIFtree
#' @export 
#' 
summary.DIFtree <-
  function(object,...){
    
    nitems       <- object$items
    dif_items    <- unique(object$splits[,"item"])
    if(!is.null(dif_items)){
      no_dif_items <- seq(1,nitems)[-dif_items]
    } else{
      no_dif_items <- seq(1,nitems)
    }
    beta_hat_dif <- object$beta_hat_dif
    
    # compute number of nodes 
    if(!is.null(dif_items)){
      endnodes <- list()
      for(i in 1:length(dif_items)){
        item               <- dif_items[[i]]
        info               <- object$splits[which(object$splits[,"item"]==item),]
        endnodes[[i]]      <- list()
        endnodes[[i]][[1]] <- 1
        for(j in 1:nrow(info)){
          endnodes[[i]][[j+1]] <- numeric(length=(j+1))
          what <- c(info[j,"left"],info[j,"right"])
          delete <- info[j,"number"]
          where  <- which(endnodes[[i]][[j]]==delete)
          endnodes[[i]][[j+1]][c(where,where+1)] <- what
          endnodes[[i]][[j+1]][-c(where,where+1)] <- endnodes[[i]][[j]][-where]
        }
      }
      names(beta_hat_dif) <- paste("Item", dif_items)
      for(i in 1:length(dif_items)){
        beta_hat_dif[[i]] <- round(beta_hat_dif[[i]],3)
        names(beta_hat_dif[[i]]) <- paste("Node",endnodes[[i]][[length(endnodes[[i]])]])
      }
    }
    cat("\n")
    cat("Call:\n",paste(deparse(object$call)),"\n")
    cat("\n")
    cat("----------\n")
    cat("\n")
    if(!is.null(dif_items)){
      cat("Items with DIF:\n", dif_items, "\n")
    } else{
      cat("Items with DIF:\n no DIF item\n")
    }
    if(!is.null(dif_items)){
      cat("Estimated item difficulties:\n")
    } else{
      cat("Estimated item difficulties:\n no DIF item\n")
      cat("\n")
    }
    if(!is.null(dif_items)){
      for(i in 1:length(dif_items)){
        cat(names(beta_hat_dif)[[i]],"\n")
        print(beta_hat_dif[[i]])
        cat("\n")
      }
    }
    cat("----------\n")
    cat("\n")
    cat("Items without DIF:\n", no_dif_items, "\n")
    cat("Estimated item difficulties:\n")
    print(round(object$beta_hat_nodif,3))
    
    invisible(object)
  }
