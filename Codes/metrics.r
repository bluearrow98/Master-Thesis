# True positive
tp <- function(comp,gt){return (length(intersect(which(comp != 0),which(gt != 0))))}
# False positive
fp <- function(comp,gt){return (length(setdiff(which(comp != 0),which(gt != 0))))}
# True negatives
tn <- function(comp,gt){return (length(intersect(which(comp == 0),which(gt == 0))))}
# False negatives
fn <- function(comp,gt){return (length(setdiff(which(comp == 0),which(gt == 0))))}
# Accuracy
acc <- function(comp,gt){
  accuracy <- (tp(comp,gt) + tn(comp,gt)) /(tp(comp,gt) + tn(comp,gt) 
                                            + fp(comp,gt) + fn(comp,gt))
  return(accuracy)
  }
# True positive rate or recall  
tpr <- function(comp,gt){
  return(tp(comp, gt) / (tp(comp, gt) + fn(comp, gt)))
}
# False positive rate
fpr <- function(comp,gt){
  return(fp(comp, gt) / (fp(comp, gt) + tn(comp, gt)))
}
# Precision
precision <- function(comp,gt){
  return(tp(comp,gt) / (tp(comp,gt) + fp(comp,gt)))
}
# F1 score
f1score <- function(comp,gt){
  return(tp(comp, gt) / (tp(comp, gt) + 0.5 * (fp(comp, gt) + fn(comp, gt))))
}
# Sigmoid function
sigmoid <- function(x){
  return(exp(abs(x)) / (1 + exp(abs(x))))
}


# Compute area under precision recall curve
AUCPR <- function(comp, gt){

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

  auc <- AUCPR_path(x$precision, x$recall)

  return(auc)
}

# Compute area under receiver operating characteristic curve
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

  auc <- AUCROC_path(x$tpr, x$fpr)

  return(auc)
}

AUCROC_path <- function(tpr, fpr, num_points = 10){

  sorted_indices <- order(fpr)
  sorted_tpr <- tpr[sorted_indices]
  sorted_fpr <- fpr[sorted_indices]

  # Extrapolation by max value
  if (max(sorted_fpr) != 1) {
    sorted_tpr <- c(sorted_tpr, max(tpr))
    sorted_fpr <- c(sorted_fpr, 1)
  }
  if (min(sorted_fpr) != 0) {
    sorted_tpr <- c(max(tpr), sorted_tpr)
    sorted_fpr <- c(0, sorted_fpr)
  }

  if (all(sorted_fpr[order(sorted_tpr)] == sorted_fpr)) {
  sorted_tpr <- sorted_tpr[order(sorted_tpr)]
  }

  # Interpolation by nearest neighbour 
  fpr_vals <- unique(sorted_fpr)
  if (length(fpr_vals) < num_points) {
    set.seed(length(fpr_vals))

    new_fpr <- runif(num_points, min = min(fpr_vals), max = max(fpr_vals))
    interp_fpr <- c(sorted_fpr, new_fpr)
    interp_tpr <- c(sorted_tpr, sapply(new_fpr, function(x) {
      max(sorted_tpr[which(abs(x - sorted_fpr) ==
     min(abs(x - sorted_fpr)))])
    }))
    newSortInd <- order(interp_fpr)
    sorted_tpr <- interp_tpr[newSortInd]
    sorted_fpr <- interp_fpr[newSortInd]
  }

  # Trapezoidal quadrature rule
  auc <- sum(diff(sorted_fpr) * (sorted_tpr[-1] + sorted_tpr[-length(sorted_tpr)]) / 2)

  return(auc)
}

AUCPR_path <- function(precision, recall, num_points = 10){

sorted_indices <- order(recall)
sorted_precision <- precision[sorted_indices]
sorted_recall <- recall[sorted_indices]

# Extrapolation by max value
if (max(sorted_precision) != 1) {
    sorted_precision <- c(sorted_precision, max(precision))
    sorted_recall <- c(sorted_recall, 1)
}
if (min(sorted_recall) != 0) {
    sorted_precision <- c(max(precision), sorted_precision)
    sorted_recall <- c(0, sorted_recall)
}

if (all(sorted_recall[order(-sorted_precision)] == sorted_recall)) {
  sorted_precision <- sorted_precision[order(-sorted_precision)]
}

# Interpolation by nearest neighbour 
recall_vals <- unique(sorted_recall)
if (length(recall_vals) < num_points) {
  set.seed(length(recall_vals))

  new_recall <- runif(num_points, min = min(recall_vals), max = max(recall_vals))
  interp_recall <- c(sorted_recall, new_recall)
  interp_precision <- c(sorted_precision, sapply(new_recall, function(x) {
    max(sorted_precision[which(abs(x - sorted_recall) ==
     min(abs(x - sorted_recall)))])
  }))
  newSortInd <- order(interp_recall)
  sorted_precision <- interp_precision[newSortInd]
  sorted_recall <- interp_recall[newSortInd]
}

auc <- sum(diff(sorted_recall) * (sorted_precision[-1] + sorted_precision[-length(sorted_precision)]) / 2)

return(auc)
}
