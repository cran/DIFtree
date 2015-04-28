#' Plotting Item focussed Trees
#' 
#' @description
#' Visualization of trees for items with DIF identified by item focussed recursive partitioning.  
#' 
#' @param x object of class \code{\link[DIFtree]{DIFtree}}
#' @param item number of the item, for which the tree shall be plotted 
#' @param cex.lines line width of lines of the tree
#' @param cex.main size of the title of the tree
#' @param cex.branches size of labeling of branches of the tree 
#' @param cex.coefs size of coefficients in the terminal nodes of the tree
#' @param title optional title, which is added to the tree, 
#' if \code{title=NULL} the title is the number of the plotted item.
#' @param ... further arguments passed to or from other methods
#' 
#' @author Moritz Berger <moritz.berger@@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Department of Statistics, LMU Munich
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}
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
#' plot(mod,item=1)
#' }
#
#' @method plot DIFtree
#' @export
#' @importFrom plotrix draw.ellipse

plot.DIFtree <-
function(x,item,cex.lines=2,cex.main=1,cex.branches=1,cex.coefs=1,title=NULL,...){
  
  #library("plotrix")
  
  if(is.null(x$splits)){
    cat("There is no plot available in the case of no DIF item")
  } else{
  
    info     <- x$splits[which(x$splits[,"item"]==item),]
  
    if(nrow(info)==0){
    
      beta_item <- x$beta_hat_nodif[paste0("beta",item)]
      cat("Item", item, "is no DIF item. There is no tree to plot.\n")
      cat("Estimated item difficulty:", beta_item)
    
    } else{
    
      # compute number of nodes 
      endnodes      <- list()
      endnodes[[1]] <- 1
      for(j in 1:nrow(info)){
        level <- info[j,"level"]
        node  <- info[j,"node"]
        endnodes[[j+1]] <- numeric(length=(j+1))
        what <- max(endnodes[[j]])+c(1,2)
        delete <- endnodes[[level]][node]
        where  <- which(endnodes[[j]]==delete)
        endnodes[[j+1]][c(where,where+1)] <- what
        endnodes[[j+1]][-c(where,where+1)] <- endnodes[[j]][-which(endnodes[[j]]==delete)]
      }
      endnodes <- endnodes[[length(endnodes)]]
  
      n_levels <- length(unique(info[,"level"]))
      
      hilfspunkte <- list()
      hilfspunkte[[1]] <- matrix(NA,nrow=2^n_levels,ncol=2)
      hilfspunkte[[1]][,1] <- 2^n_levels
      hilfspunkte[[1]][,2] <- rep(n_levels+1,2^n_levels)
      
      steps <- 2^((n_levels:1-1))
      
      for(i in 1:n_levels){
        
        hilfspunkte[[i+1]] <- hilfspunkte[[i]]
        hilfspunkte[[i+1]][,2] <- rep(n_levels+1-i,2^n_levels)
        
        help  <- c(-steps[i],steps[i])
        help1 <- rep(help,each=steps[i])
        help2 <- rep(help1,length=2^n_levels) 
        hilfspunkte[[i+1]][,1] <- hilfspunkte[[i]][,1]+help2 
        
        which_knots <- info[info[,"level"]==i,"node"]
        help3 <- seq(1,2^n_levels)
        help4 <- split(help3,rep(1:2^(i-1),each=2^n_levels/2^(i-1)))
        help5 <- unlist(lapply(which_knots, function(j) help4[[j]]))
        hilfspunkte[[i+1]][-help5,] <- hilfspunkte[[i]][-help5,]
    
      }
      
      
      plot.new()
      plot.window(ylim=c(0.5,n_levels+1),xlim=c(0,2^(n_levels+1)))
      rect(0,0.5,2^(n_levels+1),n_levels+1, border = grey(0.9),col = grey(0.9))
      
      
      for(j in length(hilfspunkte):2){
        for(i in 1:(2^n_levels)){
          lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
                lwd=cex.lines)
        }
      }
      if(is.null(title)){
        title <- paste("Item",item)
      }
      title(title,cex.main=cex.main)
    
      # Fuege Schaetzer in den Knoten hinzu
      beta_hat <- round(x$beta_hat_dif[[which(names(x$beta_hat_dif)==item)]],3)
      points_beta <- unique(hilfspunkte[[n_levels+1]])
      points_beta[,2] <- points_beta[,2]-0.2
      for(i in 1:length(beta_hat)){
        draw.ellipse(x=points_beta[i,1],y=points_beta[i,2],a=0.8,b=0.2,lwd=cex.lines,col=grey(0.8))
      }
      text(points_beta[,1],points_beta[,2],beta_hat,cex=cex.coefs)
      
      # Fuege Knotennummern hinzu 
      points_beta[,2] <- points_beta[,2]+0.2
      for(i in 1:length(beta_hat)){
      rect(points_beta[i,1]-0.2,points_beta[i,2]-0.1,points_beta[i,1]+0.2,points_beta[i,2]+0.1,col=grey(0.9),lwd=cex.lines)
      text(points_beta[i,1],points_beta[i,2],endnodes[i],cex=cex.branches)
      }
      
      # Fuege weitere Beschriftung hinzu 
      for(i in 1:nrow(info)){
        help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
        point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
        points(point_var[1],point_var[2],cex=cex.lines-1,pch=19)
        point_left  <- c(point_var[1]-steps[info[i,"level"]],point_var[2]-0.5)
        point_right <- c(point_var[1]+steps[info[i,"level"]],point_var[2]-0.5)
        var   <- info[i,"variable"]
        thres <- info[i,"threshold"]
        sort_values <- unique(sort(x$X[,var]))
        if(thres==min(sort_values)){
          text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches)
        } else{
          text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches)
        }
        if(thres==max(sort_values[-length(sort_values)])){
          text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches)
        } else{
          text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches)
        }
      }
    }
    invisible(x)
  }

}
