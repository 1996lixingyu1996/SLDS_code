#' @param x A Vector
#' @return numbers of obs in each factor levels.
#' @export
n_levels <- function(x){
  out <- NA
  if( is.factor(x) ){
    out <- sapply(X=levels(x), FUN=function(X,s) length(which(s==X)), s=x)
  } else if( is.character(x) ){
    x <- factor(x)
    out <- sapply(X=levels(x), FUN=function(X,s) length(which(s==X)), s=x)
  }
  return( out )
}

#' @param x A Vector

to_numeric <- function(x){ return(as.numeric(as.character(x))) }

#' describe splits
#'
#' @param p A DDD
#' @param k A DDD
#' @param v A DDD
#' @importFrom stringr str_sub str_length
#' @return DDDD
#' @export

parse_split <- function(p,k,v){

  s <- c(">=","< ","> ","<=","=")

  out <- temp <- NULL
  try( temp <- v[ which( sapply( v, FUN=function(x,d){ return( grepl( pattern=x, x=d, fixed=TRUE ) ) }, d=p[k] ) ) ], silent=FALSE )

  if( length(temp)>1 ){
    temp <- temp[which.max(stringr::str_length(temp))]
  }

  if( length(temp)>0 ){
    out <- temp
    temp <- NULL
    temp <- s[ which( sapply(s, FUN=function(x,d){ return( grepl( pattern=x, x=d, fixed=TRUE ) ) },
                             d=stringr::str_sub( p[k], stringr::str_length(out[1])+1, stringr::str_length(p[k]) ) ) ) ][1]

    if( length(temp)>0 ){

      out <- c(out, temp)

      temp <- NULL
      temp <- stringr::str_sub( p[k], stringr::str_length(out[1])+stringr::str_length(out[2])+1, stringr::str_length(p[k]) )

      if( length(temp)>0 ){
        out <- c(out,temp)
      }

    }

  }

  if( length(out)==3 ){

    if(out[2]=="="){ out[2] <- "==" }

  } else{

    out <- NULL

  }

  return(out)
}


#' @param tree A DDD
#' @param newdata A DDD
#' @importFrom  treeClust rpart.predict.leaves
#' @export
tree_where <- function(tree, newdata){

  out <- NULL

  if( inherits(tree,"rpart") ){

    suppressWarnings( try( out <- treeClust::rpart.predict.leaves( tree, newdata=newdata, type="where" ), silent=TRUE ) )

  } else{  ### Random assignment until Tree Build is complete ####

    #out <- factor(round(runif(nrow(newdata),0.51,4.49)))
    #names(out) <- rownames(newdata)

  }

  return(out)
}


#' @param tree A DDD
#' @param nodes A DDD
#' @param pretty A DDD
#' @param print.it A DDD
#' @importFrom stringr str_sub str_locate str_length
#' @import rpart

path_rpart <- function (tree, nodes, pretty = 0, print.it = FALSE){

  comp_digits <- function(tree){
    out <- 5
    oo <- capture.output( tree, type="output" )
    u <- grep(">",oo)
    if( length(u)>0 ){
      out1 <- max( sapply(u, function(x,p)
        stringr::str_length( stringr::str_sub( stringr::str_sub(p[x], stringr::str_locate(p[x],"\\>")[1]+2, stringr::str_length(p[x])), 1,
                                               stringr::str_locate( stringr::str_sub(p[x], stringr::str_locate(p[x],"\\>")[1]+2, stringr::str_length(p[x])), " ")[1]-1 ) ), p=oo) )+1
      if(out1>out){ out <- out1 }
    }
    return(out)
  }
  ##########################################

  if (!inherits(tree, "rpart"))
    stop("Not a legitimate \"rpart\" object")

  splits <- labels(tree, pretty = pretty, digits=comp_digits(tree))
  frame <- tree$frame
  n <- row.names(frame)
  node <- as.numeric(n)
  which <- rpart:::descendants(node)
  path <- list()

  if (length(nodes <- rpart:::node.match(nodes, node)) == 0L)
    return(invisible())
  for (i in nodes) {
    path[[n[i]]] <- path.i <- splits[which[, i]]
    if (print.it) {
      cat("\n", "node number:", n[i], "\n")
      cat(paste("  ", path.i), sep = "\n")
    }
  }
  invisible(path)
}

#' predict subtype and leaf assignments
#'
#' @param tree A decision tree from rpart object
#' @param newdata data.frame
#' @importFrom stats predict
#' @export

predict_path <- function(tree, newdata, newvarname="subtype"){

  temp <- NULL
  try( temp <- stats::predict(tree, newdata=newdata, type="class"), silent=TRUE )

  if( length(temp)>0 ){
    newdata[,paste0(newvarname)] <- as.factor(temp)
  }

  temp <- NULL
  try( temp <- tree_where(tree, newdata=newdata), silent=TRUE )

  if( length(temp)>0 ){
    newdata[,paste0("leaf")] <- as.factor(temp)
  }

  return(newdata)
}

#' Describe leaf splits
#'
#' @param tree A DDD
#' @param DAT data.frame
#' @importFrom stringr str_sub str_locate str_length
#' @export

get_valid_levels <- function(info){
  flag = 0
  info_summary = summary(info)
  for (i in 1:length(info_summary)){
    if (info_summary[i] > 0){
      flag = flag + 1
    }
  }
  return (flag)
}



#' @export
profile_leafs <- function(tree){

  oo <- utils::capture.output( lpath <- path_rpart(tree=tree, nodes=sort(as.numeric(rownames(tree$frame[as.numeric(names(n_levels(attributes(tree)[["traindata"]][, "leaf"]))),])))), type="output" )

  return(lpath)
}


#' Plot decsion tree (rpart decision tree object)
#'
#' @importFrom rpart.plot rpart.plot
#' @export

plot_profile <- function(tree){
  suppressWarnings( rpart.plot::rpart.plot(tree,
                                           type=5, clip.right.labs=FALSE, branch=.5,extra = 100,
                                           main = paste("Subtype Profile") ) )
}




#' Fit trees and sort by CV error
#'
#' @param Y vector or data.frame dependent vars
#' @param X data.frame independent vars
#' @param seed numeric seed
#' @param cv_select logical use CV to select best model
#' @param split split rules to evaluate
#' @param shrink hyperparameter
#' @param minsplit hyperparameter
#' @param cp hyperparameter
#' @param maxcompete hyperparameter
#' @param maxsurrogate hyperparameter
#' @param usesurrogate hyperparameter
#' @param maxdepth hyperparameter
#' @param tList hyperparameter
#' @param min.p.thresh hyperparameter
#' @importFrom survival Surv
#' @import rpart
#' @return list of trees
#' @export

tree_fit <- function(Y, X, seed=1234, cv_select=TRUE,
                     split = c("gini","information"),
                     shrink = 1, minsplit = c(0.10,0.15,0.2),
                     cp = c( 0.000001, seq(0.000005, 0.00005, 0.000005), 0.0001, 0.0005, 0.001, 0.005, 0.01 ),
                     maxcompete = 0,
                     maxsurrogate = 0,
                     min_leaf = 1,
                     usesurrogate = 1,
                     maxdepth = 30, tList=NULL, min.p.thresh=NULL ){

  if( is.null(tList) ){
    tList = list( split = split,
                  shrink = shrink,
                  minsplit = minsplit,
                  minbucket = min_leaf,
                  cp = cp,
                  maxcompete = maxcompete,
                  maxsurrogate = maxsurrogate,
                  usesurrogate = usesurrogate,
                  maxdepth = maxdepth ) }

  for(m in 1:length(tList$minsplit) ){
    if( tList$minsplit[m] < 1 ){
      tList$minsplit[m] <- round(length(Y)*tList$minsplit[m])
    }
  }

  tree_fit_1 <- function(arg, Y, Obs, cv_select){

    NAM <- colnames(X)

    xval <- 10
    if(!cv_select){ xval <- 0 }

    out <- fit <- NULL
    if( is.null(ncol(Y)) ){
      formu <- paste0("Y ~ ", NAM[1])
      if(length(NAM) > 1) {
        for(g in 2:length(NAM)) {
          formu <- paste(formu, " + ", NAM[g], sep = "")
        }
      }
      DAT <- as.data.frame( X, stringsAsFactors = TRUE )
      DAT$Y <- Y

      if( is.factor(DAT$Y) ){
        oo <- capture.output( try( fit <- rpart::rpart( formu, data=DAT, method="class",
                                                        parms=list( split=arg$split ),
                                                        control=rpart::rpart.control(minsplit=arg$minsplit,
                                                                                     minbucket=arg$minbucket,
                                                                                     cp=arg$cp,
                                                                                     maxcompete=arg$maxcompete,
                                                                                     maxsurrogate=arg$maxsurrogate,
                                                                                     usesurrogate=arg$usesurrogate,
                                                                                     maxdepth=arg$maxdepth, xval=xval ) ), silent=TRUE ), type="output" )
      } else{
        oo <- capture.output( try( fit <- rpart::rpart( formu, data=DAT, method="anova",
                                                        control=rpart::rpart.control(minsplit=arg$minsplit,
                                                                                     minbucket=arg$minbucket,
                                                                                     cp=arg$cp,
                                                                                     maxcompete=arg$maxcompete,
                                                                                     maxsurrogate=arg$maxsurrogate,
                                                                                     usesurrogate=arg$usesurrogate,
                                                                                     maxdepth=arg$maxdepth, xval=xval ) ), silent=TRUE ), type="output" )
      }
    } else{
      colnames(Y) <- c("Y","E")
      formu <- paste0("survival::Surv(Y,E) ~ ", NAM[1])
      if(length(NAM) > 1) {
        for(g in 2:length(NAM)) {
          formu <- paste(formu, " + ", NAM[g], sep = "")
        }
      }
      DAT <- as.data.frame( X, stringsAsFactors = TRUE )
      DAT$Y <- Y[,1]
      DAT$E <- Y[,2]

      oo <- capture.output( try( fit <- rpart::rpart( formu, data=DAT, method="exp",
                                                      parms=list( shrink=arg$shrink ),
                                                      control=rpart::rpart.control(minsplit=arg$minsplit,
                                                                                   minbucket=arg$minbucket,
                                                                                   cp=arg$cp,
                                                                                   maxcompete=arg$maxcompete,
                                                                                   maxsurrogate=arg$maxsurrogate,
                                                                                   usesurrogate=arg$usesurrogate,
                                                                                   maxdepth=arg$maxdepth, xval=xval ) ), silent=TRUE ), type="output" )
    }

    if( !is.null(fit) ){
      if( cv_select ){ pfit <- NULL
      oo <- capture.output( try( pfit <- rpart::prune(fit, cp = fit$cptable[ which.min( fit$cptable[,"xerror"] ), 1 ] ), silent=TRUE  ), type="output" )
      if( !is.null(pfit) ){
        temp <- NULL
        try( temp <- predict_path(tree=pfit, newdata=DAT), silent=TRUE )
        if( length(temp)>0 ){
          attributes(pfit)[["traindata"]] <- temp
          if( "leaf" %in% names(attributes(pfit)[["traindata"]]) ){
            if( length(levels(attributes(pfit)[["traindata"]]$leaf))>1 ){
              out <- list( key=paste0("Tree_",paste(unlist(arg),collapse="-")), fit=pfit, xerror=min(fit$cptable[,"xerror"]),
                           min.p=min( n_levels(factor(pfit$where))/length(pfit$where) ) )
            }
          }
        }
      }
      } else{
        temp <- NULL
        try( temp <- predict_path(tree=fit, newdata=DAT), silent=TRUE )
        if( length(temp)>0 ){
          attributes(fit)[["traindata"]] <- temp
          if( "leaf" %in% names(attributes(pfit)[["traindata"]]) ){
            if( length(levels(attributes(fit)[["traindata"]]$leaf))>1 ){
              out <- list( key=paste0("Tree_",paste(unlist(arg),collapse="-")), fit=fit, xerror=fit$cptable[nrow(fit$cptable),"rel error"],
                           min.p=min( n_levels(factor(fit$where))/length(fit$where) ) )
            }
          }
        }
      }
    }

    return( out )
  }
  ###################

  set.seed(seed)
  out <- NULL

  tList. <- list( minsplit = tList$minsplit,
                  minbucket = tList$minbucket,
                  cp = tList$cp,
                  maxcompete = tList$maxcompete,
                  maxsurrogate = tList$maxsurrogate,
                  usesurrogate = tList$usesurrogate,
                  maxdepth = tList$maxdepth )

  if( !is.null(ncol(Y)) ){ tList.$shrink <- tList$shrink
  } else if( is.factor(Y) ){  tList.$split <- tList$split }

  runList <- expand.grid( tList. )

  min.p <- xerror <- splitvar <- tnode <- n.tnode <- NULL; fit <- list()
  for(k in 1:nrow(runList)){
    temp <- NULL
    oo <- capture.output( try( temp <- tree_fit_1( arg=as.list(runList[k,]), Y=Y, Obs=X, cv_select=cv_select ), silent=TRUE ), type="output" )
    if(!is.null(temp)){
      tnode. <- paste(temp$fit$where, collapse=",")
      if( sum( tnode==tnode. )==0 ){
        tnode <- c(tnode, tnode. )
        n.tnode <- c(n.tnode, length(levels(factor(temp$fit$where))) )
        xerror <- c(xerror, temp$xerror)
        min.p <- c(min.p, temp$min.p)
        splitvar <- c( splitvar, paste(sort(unique(rownames(temp$fit$splits))), collapse=", ") )
        #splitvar <- c( splitvar, paste(sort(unique(rownames(summary(temp$fit)$splits))), collapse=", ") )
        fit[[temp$key]] <- temp$fit
      }
    }
  }
  if( length(n.tnode)<1 ){ return(NULL) } #stop("tree_fit failed to produce trees") }
  stats <- rep(NA,4); s.fit <- list()
  o <- order(n.tnode, decreasing=FALSE)
  for(k in 1:length(o)){
    stats <- rbind(stats, c(xerror[o[k]], n.tnode[o[k]], min.p[o[k]], splitvar[o[k]]) )
    s.fit[[k]] <- fit[[o[k]]]
    rownames(stats)[nrow(stats)] <- names(fit)[o[k]]
    names(s.fit)[k] <- names(fit)[o[k]]
  }
  stats <- as.data.frame(stats, stringsAsFactors = TRUE)[-1,]
  stats[,1] <- to_numeric(stats[,1])
  stats[,2] <- to_numeric(stats[,2])
  stats[,3] <- to_numeric(stats[,3])
  names(stats) <- c("xerror","t.node.size","min.p","split.vars")

  if( length(min.p.thresh)>0 ){
    stats <- base::subset(stats, stats[,3] > min.p.thresh)
  }

  v <- as.numeric(as.character(levels(factor(stats[,2]))))
  if( nrow(stats)>1 ){

    keep <- rep(0, nrow(stats))
    for(k in v){
      stats. <- stats[ which(stats[,2]==k), ]
      if( nrow(stats.)==1 ){ keep[ which( rownames(stats)==rownames(stats.) ) ] <- 1
      } else{
        keep[ which( rownames(stats)==rownames(stats.)[which.min(stats.[,1])] ) ] <- 1
      }
    }
    if( sum(keep)>0 ){
      trim.stats <- stats[ which(keep==1), ]

      rm <- rep(0,nrow(trim.stats)); c.min <- Inf
      for(k in 1:length(rm)){
        if( c.min > trim.stats[k,1] ){ c.min <- trim.stats[k,1]
        } else{ rm[k] <- 1
        }
      }

      final.stats <- NULL
      if( sum(1-rm)==length(rm) ){
        final.stats <- trim.stats
      } else{
        final.stats <- trim.stats[-which(rm==1),]
      }
      if( !is.null(final.stats) ){
        use.i <- which( rownames(stats) %in% rownames(final.stats) )
        trees <- list()
        for(k in 1:length(use.i)){
          trees[[k]] <- s.fit[[ use.i[k] ]]
          names(trees)[k] <- names(s.fit)[use.i[k]]
        }

        out <- list( trees=trees, stats=final.stats,
                     best.tree=trees[[ length(trees) ]],
                     best.tree.vars=sort(unique(rownames(trees[[ length(trees) ]]$splits))) )
      }
    }
  } else{

    out <- list( trees=s.fit, stats=stats,
                 best.tree=s.fit[[1]],
                 best.tree.vars=sort(unique(rownames(s.fit[[1]]$splits))) )
  }

  return( out )
}



#' @export
count_var_in_single_path <- function(single_path, var_names){
  n_leaf = length(single_path)  # length of this path(including root)
  n_var_in_profile = length(var_names)
  count_var = rep(0,n_var_in_profile)
  for (j in 1:n_var_in_profile)
  {
    # count variable emerging times in the string
    count_var[j] = sum(stringr::str_detect(single_path, var_names[j]))
  }
  count_var
}


#' @export
maximum_var_split <- function(profile){
  var_names = names(profile$variable.importance)
  leafs = profile_leafs(profile)
  n_cluster = length(leafs)
  maximum_var_used = 0
  for (i in 1:n_cluster){
    single_path_result = count_var_in_single_path(single_path = leafs[[i]], var_names = var_names)
    maximum_var_used = max(maximum_var_used, max(single_path_result))
    #cat(maximum_var_used)
  }

  maximum_var_used
}

has_repeated_vars <- function(fit) {
  frame <- fit$frame
  nodes <- as.numeric(row.names(frame))
  vars <- frame$var
  leaf_nodes <- nodes[frame$var == "<leaf>"]

  for (leaf in leaf_nodes) {
    path_vars <- c()
    node <- leaf

    while (node >= 1) {
      if (node %in% nodes) {
        v <- vars[match(node, nodes)]
        if (v != "<leaf>") path_vars <- c(path_vars, v)
      }
      node <- floor(node / 2)
    }

    #
    if (any(duplicated(path_vars))) return(TRUE)
  }

  return(FALSE)
}
