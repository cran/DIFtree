ptree_PCM <- 
function(info,
                      item,
                      params,
                      X,
                      cex.lines,
                      cex.main,
                      cex.branches,
                      cex.coefs,
                      title){
  
    n_levels <- length(unique(info[,"level"]))
    n_splits <- nrow(info)
    
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
  
    # add infos 
    for(i in 1:nrow(info)){
      help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
      point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
      points(point_var[1],point_var[2],cex=cex.lines-1,pch=19)
      point_left  <- c(point_var[1]-steps[info[i,"level"]],point_var[2]-0.5)
      point_right <- c(point_var[1]+steps[info[i,"level"]],point_var[2]-0.5)
      var   <- info[i,"variable"]
      thres <- info[i,"threshold"]
      sort_values <- unique(sort(X[,var]))
      if(thres==min(sort_values)){
        text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches,adj=c(1,0))
      } else{
        text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches,adj=c(1,0))
      }
      if(thres==max(sort_values[-length(sort_values)])){
        text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches,adj=c(0,0))
      } else{
        text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches,adj=c(0,0))
      }
    }
  
    params <- matrix(params,byrow=TRUE,nrow=n_splits+1)
    
    yscale <- c(min(params)-1,max(params)+1)
    if(min(params) < -4){
      yscale[1] <- -4 
    }
    if(max(params) > 4){
      yscale[2] <- 4 
    }
    
    vps <- baseViewports()
    pushViewport(vps$plot)
    
  
    points_params <- unique(hilfspunkte[[n_levels+1]])
    for(i in 1:nrow(points_params)){
      window <- viewport(x=points_params[i,1],y=points_params[i,2],
                         width=1,height=0.5,
                         just=c("center","top"),
                         xscale=c(0,ncol(params)+1),yscale=yscale,
                         default.units="native")
      pushViewport(window)
      grid.rect()
      
      if(any(params[i,] < -4)){
        params[i,which(params[i,]< -4)] <- -4 
      }
      if(any(params[i,] > 4)){
        params[i,which(params[i,]> 4)] <- 4 
      }
      
      for(j in 1:ncol(params)){
        if(-4 < params[i,j] & params[i,j] < 4){
          grid.points(j,params[i,j],pch=19,gp=gpar(cex=cex.coefs))
        }
        if(j < ncol(params)){
          grid.lines(c(j,j+1),params[i,(j:(j+1))],default.units="native",gp=gpar(lwd=cex.coefs*3))
        }
      }
      grid.yaxis(gp=gpar(cex=cex.branches))
      upViewport()
    }
}





