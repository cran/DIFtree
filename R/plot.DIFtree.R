#' Plotting of Item focussed Trees
#' 
#' @description
#' Visualization of trees for items with DIF identified by item focussed recursive partitioning 
#' in dichotomous or polytomous items.
#' 
#' @param x Object of class \code{\link[DIFtree]{DIFtree}}
#' @param item Number of the item, for which the tree shall be plotted 
#' @param component Component of the model for which the tree shall be plotted; 
#' can be \code{"intercept"} or \code{"slope"}. For \code{"Rasch"} and \code{"PCM"} only one tree of item parameters
#' is available for each DIF item and therefore \code{component} will be ignored. 
#' @param cex.lines Width of branches of the tree
#' @param cex.branches Size of the labels of branches of the tree 
#' @param cex.coefs Size of coefficients in the terminal nodes of the tree
#' @param cex.main Size of the title of the tree
#' @param title Optional title, which is added to the tree;
#' if \code{title=NULL} the title is the number of the plotted item.
#' @param ... Further arguments passed to or from other methods
#' 
#' @author Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2016): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Journal of Educational and Behavioral Statistics 41(6), 559-592.
#' 
#' Bollmann, Stella, Berger, Moritz & Tutz, Gerhard (2018): Item-Focussed Trees for the Detection 
#' of Differential Item Functioning in Partial Credit Models, Educational and Psychological Measurement 78(5), 781-804.
#' 
#' Tutz, Gerhard and Berger, Moritz (2016): Item focussed Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika 81(3), 727-750. 
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
#' 
#' @examples 
#' data(data_sim_Rasch)
#'  
#' Y <- data_sim_Rasch[,1]
#' X <- data_sim_Rasch[,-1]
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#'  
#' plot(mod,item=1)
#' }
#
#' @method plot DIFtree
#' @export
#' @import grid
#' @import gridBase
#' @importFrom plotrix draw.ellipse
#' @importFrom grDevices grey
#' @importFrom graphics lines plot.new plot.window points rect text 

plot.DIFtree <-
function(x, # object of class DIFtree
                         item,
                         component="intercept",
                         cex.lines=2,cex.branches=1,cex.coefs=1,cex.main=1,title=NULL,...){
  
  if(missing(item)){
    stop("argument \"item\" is missing, with no default")
  }
  
  if(is.null(x$splits)){
    cat("There is no plot available in the case of no DIF item")
  } else{
    
    X <- x$X
    
    model <- which(c("Rasch","Logistic","PCM")==x$model)
    if(model==1){
      info     <- x$splits[which(x$splits[,"item"]==item),]
      if(nrow(info)==0){
        beta_item <- x$coefficients$betas_nodif[paste0("beta",item)]
        cat("Item", item, "is no DIF item. There is no tree to plot.\n")
        cat("Estimated item difficulty:", beta_item)
      } else{
        betas_item <- x$coefficients$betas_dif[[which(names(x$coefficients$betas_dif)==item)]]
        ptree(info,item,betas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
      }
    }
    if(model==2){
      type <- which(c("udif","dif","nudif")==x$type)
      if(type==1){
        info     <- x$splits[which(x$splits[,"item"]==item),]
        if(nrow(info)==0){
          gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
          cat("Item", item, "is no DIF item. There is no tree to plot.\n")
          cat("Estimated intercept:", gamma_item)
        } else{
          gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
          ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
        }
      }
      if(type==2){
      
        dif_items <- unique(c(x$splits[[1]][,"item"],x$splits[[2]][,"item"]))
        if(is.null(dif_items)){
          cat("There is no plot available in the case of no DIF item")
        } else{
          if(!item %in% dif_items){
            gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
            alpha_item <- x$coefficients$alphas_nodif[paste0("alpha",item)]
            cat("Item", item, "is no DIF item. There are no trees to plot.\n")
            cat("Estimated intercept:", gamma_item,"\n")
            cat("Estimated slope:", alpha_item)
          } else{
            
            if(!(component %in% c("intercept","slope"))){
              stop(paste("Component",component,"undefined. Must be 'intercept' or 'slope'."))
            }
            
            if(component=="intercept"){
              info     <- x$splits[[1]][which(x$splits[[1]][,"item"]==item),]
              if(is.null(info) || nrow(info)==0){
                gamma_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
                cat("Item", item, "has no DIF in intercept. There is no tree to plot.\n")
                cat("Estimated intercept:", gamma_item)
              } else{
                gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
                ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
              }
            } else{
              info     <- x$splits[[2]][which(x$splits[[2]][,"item"]==item),]
              if(is.null(info) || nrow(info)==0 ){
                alpha_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
                cat("Item", item, "has no DIF in slope. There is no tree to plot.\n")
                cat("Estimated slope:", alpha_item)
              } else{
                alphas_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
                ptree(info,item,alphas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
              }
            }
          }
        }
      }
      if(type==3){
        
        info     <- x$splits[which(x$splits[,"item"]==item),]
        if(nrow(info)==0){
          gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
          alpha_item <- x$coefficients$alphas_nodif[paste0("alpha",item)]
          cat("Item", item, "is no DIF item. There are no trees to plot.\n")
          cat("Estimated intercept:", gamma_item,"\n")
          cat("Estimated slope:", alpha_item)
        } else{
          
          if(!(component %in% c("intercept","slope"))){
            stop(paste("Component",component,"undefined. Must be 'intercept' or 'slope'."))
          }
          
          if(component=="intercept"){
            gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
            ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
          } else{
            alphas_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
            ptree(info,item,alphas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
          }
        }
      }
    }
    if(model==3){
      info     <- x$splits[which(x$splits[,"item"]==item),]
      if(nrow(info)==0){
        deltas_item <- x$coefficients$deltas_nodif[,paste0("delta",item)]
        deltas_item <- paste(round(deltas_item,3), collapse=", ")
        cat("Item", item, "is no DIF item. There is no tree to plot.\n")
        cat("Estimated item parameters:", deltas_item)
      } else{
        deltas_item <- x$coefficients$deltas_dif[[which(names(x$coefficients$deltas_dif)==item)]]
        ptree_PCM(info,item,deltas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
      }
    }
  }
  invisible(x)
}
