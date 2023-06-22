tp <- function(comp,gt){return (length(intersect(which(comp != 0),which(gt != 0))))}
fp <- function(comp,gt){return (length(setdiff(which(comp != 0),which(gt != 0))))}
tn <- function(comp,gt){return (length(intersect(which(comp == 0),which(gt == 0))))}
fn <- function(comp,gt){return (length(setdiff(which(comp == 0),which(gt == 0))))}
acc <- function(comp,gt){
  accuracy <- (tp(comp,gt) + tn(comp,gt)) /(tp(comp,gt) + tn(comp,gt) 
                                            + fp(comp,gt) + fn(comp,gt))
  return(accuracy)
  }
tpr <- function(comp,gt){
  return(tp(comp, gt) / (tp(comp, gt) + fn(comp, gt)))
}
fpr <- function(comp,gt){
  return(fp(comp, gt) / (fp(comp, gt) + tn(comp, gt)))
}
f1score <- function(comp,gt){
  return(tp(comp, gt) / (tp(comp, gt) + 0.5 * (fp(comp, gt) + fn(comp, gt))))
}

sigmoid <- function(x){
  return(exp(abs(x)) / (1 + exp(abs(x))))
}
precision <- function(comp,gt){
  return(tp(comp,gt) / (tp(comp,gt) + fp(comp,gt)))
}


AUCPR <- function(comp,gt){

  comp <- sigmoid(abs(comp))

  x = list()
  comp_sort <- c(0,sort(comp))
  for (i in 1:length(comp_sort)) {
    comp_i <- (comp > comp_sort[i])
    x$precision[i] <- precision(comp_i, gt)
    if (is.nan(x$precision[i])) x$precision[i] <- 1
    x$recall[i] <- tpr(comp_i,gt)
    if (is.nan(x$recall[i])) x$recall[i] <- 1
  }

  # Trapezoidal quadrature rule
  auc <- sum(diff(x$recall) * (x$precision[-1] + x$precision[-length(x$precision)]) / 2)

  return(auc)
}

AUCROC <- function(comp, gt) {

  comp <- sigmoid(abs(comp))
  x <- list()
  comp_sort <- c(0, sort(comp))
  for (i in 1:length(comp_sort)) {
    comp_i <- (comp > comp_sort[i])
    x$tpr[i] <- tpr(comp_i, gt)
    if (is.nan(x$tpr[i])) x$tpr[i] <- 1
    x$fpr[i] <- fpr(comp_i,gt)
    if (is.nan(x$fpr[i])) x$fpr[i] <- 1
  }

  # Trapezoidal quadrature rule
  auc <- sum(diff(x$fpr) * (x$tpr[-1] + x$tpr[-length(x$tpr)]) / 2)

  return(auc)
}

AUCROC_path <- function(tpr, fpr){

  sorted_indices <- order(fpr)
  sorted_tpr <- tpr[sorted_indices]
  sorted_fpr <- fpr[sorted_indices]
  if (max(sorted_fpr) != 1) {
    sorted_tpr <- c(sorted_tpr, max(tpr))
    sorted_fpr <- c(sorted_fpr, 1)
  }
  if (min(sorted_fpr) != 0) {
    sorted_tpr <- c(max(fpr), sorted_tpr)
    sorted_fpr <- c(0, sorted_fpr)
  }

  # Trapezoidal quadrature rule
  auc <- sum(diff(sorted_fpr) * (sorted_tpr[-1] + sorted_tpr[-length(sorted_tpr)]) / 2)

  return(auc)
}

AUCPR_path <- function(precision, recall){

sorted_indices <- order(recall)
sorted_precision <- precision[sorted_indices]
sorted_recall <- recall[sorted_indices]

if (max(sorted_precision) != 1) {
    sorted_precision <- c(sorted_precision, max(precision))
    sorted_recall <- c(sorted_recall, 1)
}
if (min(sorted_recall) != 0) {
    sorted_precision <- c(max(precision), sorted_precision)
    sorted_recall <- c(0, sorted_recall)
}
auc <- sum(diff(sorted_recall) * (sorted_precision[-1] + sorted_precision[-length(sorted_precision)]) / 2)

return(auc)
}
