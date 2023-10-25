obtain.single.Effects <- function(model_result, component){
  beta <- summary(model_result)$coefficients[component,"Estimate"]
  beta_se <- summary(model_result)$coefficients[component,"Std. Error"]
  p <- summary(model_result)$coefficients[component,"Pr(>|t|)"]
  return(c(beta,beta_se,p))
}

# calculation of IR
Calc_IR <- function(model_result, component){
  beta <- summary(model_result)$coefficients[component,"Estimate"]
  beta_se <- summary(model_result)$coefficients[component,"Std. Error"]
  p <- summary(model_result)$coefficients[component,"Pr(>|t|)"]
  RR <- exp(beta)
  RR_low <- exp(beta - 1.96*beta_se)
  RR_up <- exp(beta + 1.96*beta_se)
  ir <- (RR-1)*100
  ir_low <- (RR_low-1)*100
  ir_up <- (RR_up-1)*100
  return(c(ir,ir_low,ir_up,p))
}
Calc_IR1 <- function(model_result, component){
  beta <- summary(model_result)$coefficients[component,"Estimate"]
  beta_se <- summary(model_result)$coefficients[component,"Std. Error"]
  p <- summary(model_result)$coefficients[component,"Pr(>|z|)"]
  RR <- exp(beta)
  RR_low <- exp(beta - 1.96*beta_se)
  RR_up <- exp(beta + 1.96*beta_se)
  ir <- (RR-1)*100
  ir_low <- (RR_low-1)*100
  ir_up <- (RR_up-1)*100
  return(c(ir,ir_low,ir_up,p))
}
Calc_IR2 <- function(beta, beta_low, beta_up){
  RR <- exp(beta)
  RR_low <- exp(beta_low)
  RR_up <- exp(beta_up)
  ir <- (RR-1)*100
  ir_low <- (RR_low-1)*100
  ir_up <- (RR_up-1)*100
  return(c(ir,ir_low,ir_up))
}

getIR_inter <- function(model, var){
  modfit = model$msmfit
  X <- model.matrix(modfit)
  dof <- nrow(X) - ncol(X)
  group2 <- paste0("psi1:",var)
  coefs <- modfit$coefficients[c("psi1",group2)]
  coefs_var <- vcov(modfit)
  inter_1 <- qt(0.975, dof) * sqrt(coefs_var['psi1','psi1'] + coefs_var[group2,group2] +
                                     2*coefs_var['psi1',group2])
  beta_G2 <- c(coefs[1]+coefs[2],coefs[1]+coefs[2]-inter_1,coefs[1]+coefs[2]+inter_1)
  ir = exp(coefs[1])
  ir_low = exp(model$ci.coef["psi1",1])
  ir_up = exp(model$ci.coef["psi1",2])
  
  r1 <- c(ir,ir_low,ir_up)
  r2 <- Calc_IR2(beta = beta_G2[1], beta_low = beta_G2[2], beta_up = beta_G2[3])
  # Cochran's Q test
  beta1 = coefs[1]
  beta2 = coefs[1]+coefs[2] 
  var1 = model$var.coef[2]
  var2 = (inter_1/qt(0.975, dof))^2
  pooled_beta = ((beta1/var1)+(beta2/var2))/(1/var1 +1/var2)
  CochranQ = ((beta1-pooled_beta)^2)/var1 + ((beta2-pooled_beta)^2)/var2
  Pvalue = pchisq(CochranQ,df = 1, lower.tail = F)
  
  Result <- rbind(r1,r2)
  colnames(Result) <- c("Est","Low","Up")
  l = list(Result = Result, CochranQ = CochranQ, Pvalue = Pvalue)
  return(l)
}



# get quantiles for variable
quantile.fn <- function(data, n.quantiles){
  data <- as.matrix(data)
  q <- matrix(0, dim(data)[1], dim(data)[2])
  I <- dim(data)[2]
  for(i in 1:I){
    a <- rank(data[,i], ties.method = "first")
    q[,i] <- cut(a, stats::quantile(a, probs = c(0:n.quantiles/n.quantiles)), include.lowest=TRUE)}
  q <- q-1
  colnames(q) <- colnames(data)
  return(q)
}

#### qgcomp.boot ####
qgcomp.boot.quasi <- function (f, data, expnms = NULL, q = 4, breaks = NULL, id = NULL, 
                               weights, alpha = 0.05, B = 200, rr = TRUE, degree = 1, seed = NULL, 
                               bayes = FALSE, MCsize = nrow(data), parallel = FALSE, parplan = FALSE, 
                               ...) 
{
  oldq = NULL
  if (is.null(seed)) 
    seed = round(runif(1, min = 0, max = 1e+08))
  newform <- terms(f, data = data)
  hasintercept = as.logical(attr(newform, "intercept"))
  class(newform) <- "formula"
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 
             0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, 
                      FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if (hasweights) {
    data$weights <- as.vector(model.weights(thecalle))
  }
  else data$weights = rep(1, nobs)
  if (is.null(expnms)) {
    expnms <- attr(newform, "term.labels")
    message("Including all model terms as exposures of interest\n")
  }
  lin = qgcomp:::checknames(expnms)
  if (!lin) 
    stop("Model appears to be non-linear and I'm having trouble parsing it:\n                  please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)) {
    oldq = q
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)) {
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if (is.null(q)) {
      nvals <- length(br[[1]]) - 1
    }
    else {
      nvals <- q
    }
    intvals <- (seq_len(nvals)) - 1
  }
  else {
    qdata <- data
    nvals = length(table(unlist(data[, expnms])))
    if (nvals < 10) {
      message("\nNote: using all possible values of exposure as the\n              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[, expnms]))))
      br <- lapply(seq_len(p), function(x) c(-1e+16, intvals[2:nvals] - 
                                               1e-16, 1e+16))
    }
    else {
      message("\nNote: using quantiles of all exposures combined in order to set\n          proposed intervention values for overall effect (25th, 50th, 75th %ile)\n        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
      intvals = as.numeric(quantile(unlist(data[, expnms]), 
                                    c(0.25, 0.5, 0.75)))
      br <- NULL
    }
  }
  if (is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  msmfit <- msm.fit.quasi(newform, qdata, intvals, expnms, rr, main = TRUE, 
                          degree = degree, id = id, weights, bayes, MCsize = MCsize, 
                          hasintercept = hasintercept, ...)
  estb <- as.numeric(msmfit$msmfit$coefficients)
  nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[, id, drop = TRUE]))
  starttime = Sys.time()
  psi.only <- function(i = 1, f = f, qdata = qdata, intvals = intvals, 
                       expnms = expnms, rr = rr, degree = degree, nids = nids, 
                       id = id, weights, MCsize = MCsize, hasintercept = hasintercept, 
                       ...) {
    if (i == 2 & !parallel) {
      timeiter = as.numeric(Sys.time() - starttime)
      if ((timeiter * B/60) > 0.5) 
        message(paste0("Expected time to finish: ", round(B * 
                                                            timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp = sort(sample(unique(qdata[, 
                                                          id, drop = TRUE]), nids, replace = TRUE)))
    names(bootids) <- id
    qdata_ <- merge(qdata, bootids, by = id, all.x = FALSE, 
                    all.y = TRUE)
    ft = msm.fit.quasi(f, qdata_, intvals, expnms, rr, main = FALSE, 
                       degree, id, weights = weights, bayes, MCsize = MCsize, 
                       hasintercept = hasintercept, ...)
    yhatty = data.frame(yhat = predict(ft$msmfit), psi = ft$msmfit$data[, 
                                                                        "psi"])
    as.numeric(c(ft$msmfit$coefficients, with(yhatty, tapply(yhat, 
                                                             psi, mean))))
  }
  set.seed(seed)
  if (parallel) {
    if (parplan) {
      oplan <- future::plan(strategy = future::multisession, workers = 4)
      on.exit(future::plan(oplan), add = TRUE)
    }
    bootsamps <- future.apply::future_lapply(X = seq_len(B), 
                                             FUN = psi.only, f = f, qdata = qdata, intvals = intvals, 
                                             expnms = expnms, rr = rr, degree = degree, nids = nids, 
                                             id = id, weights = qdata$weights, MCsize = MCsize, 
                                             hasintercept = hasintercept, future.seed = TRUE, 
                                             ...)
  }
  else {
    bootsamps <- lapply(X = seq_len(B), FUN = psi.only, f = f, 
                        qdata = qdata, intvals = intvals, expnms = expnms, 
                        rr = rr, degree = degree, nids = nids, id = id, weights = weights, 
                        MCsize = MCsize, hasintercept = hasintercept, ...)
  }
  bootsamps = do.call("cbind", bootsamps)
  hats = t(bootsamps[-c(seq_len(degree + ifelse(hasintercept, 
                                                1, 0))), , drop = FALSE])
  cov.yhat = cov(hats)
  bootsamps = bootsamps[seq_len(degree + ifelse(hasintercept, 
                                                1, 0)), , drop = FALSE]
  seb <- apply(bootsamps, 1, sd)
  covmat <- cov(t(bootsamps))
  cnms = c(paste0("psi", seq_len(nrow(bootsamps) - 1)))
  if (hasintercept) 
    cnms = c("(intercept)", cnms)
  colnames(covmat) <- rownames(covmat) <- names(estb) <- cnms
  tstat <- estb/seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 
    1 - degree
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - 
                                                                alpha/2))
  if (!is.null(oldq)) {
    q = oldq
  }
  psiidx = seq_len(degree) + ifelse(hasintercept, 1, 0)
  qx <- qdata[, expnms]
  res <- qgcomp:::.qgcomp_object(qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
                                 psi = estb[psiidx], var.psi = seb[psiidx]^2, covmat.psi = covmat[psiidx, 
                                                                                                  psiidx, drop = FALSE], ci = ci[psiidx, ], coef = estb, 
                                 var.coef = seb^2, covmat.coef = covmat, ci.coef = ci, 
                                 expnms = expnms, q = q, breaks = br, degree = degree, 
                                 bootstrap = TRUE, y.expected = msmfit$Ya, y.expectedmsm = msmfit$Yamsm, 
                                 index = msmfit$A, bootsamps = bootsamps, cov.yhat = cov.yhat, 
                                 alpha = alpha, call = origcall, hasintercept = hasintercept)
  if (msmfit$fit$family$family == "gaussian") {
    res$tstat <- tstat
    res$df <- df
    res$pval <- pval
  }
  if (msmfit$fit$family$family %in% c("binomial", "poisson","quasipoisson")) {
    res$zstat <- tstat
    res$pval <- pvalz
  }
  res
}

#### msm.fit ####
msm.fit.quasi <- function (f, qdata, intvals, expnms, rr = TRUE, main = TRUE, 
                           degree = 1, id = NULL, weights, bayes = FALSE, MCsize = nrow(qdata), 
                           hasintercept = TRUE, ...) 
{
  newform <- terms(f, data = qdata)
  nobs = nrow(qdata)
  thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("qdata", "data", names(thecall))
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 
             0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, 
                      FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if (hasweights) {
    qdata$weights <- as.vector(model.weights(thecalle))
  }
  else qdata$weights = rep(1, nobs)
  if (is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  nidx = which(!(names(qdata) %in% id))
  if (!bayes) 
    fit <- glm(newform, data = qdata, weights = weights, 
               ...)
  if (bayes) {
    requireNamespace("arm")
    fit <- bayesglm(f, data = qdata[, nidx, drop = FALSE], 
                    weights = weights, ...)
  }
  if (fit$family$family %in% c("gaussian", "poisson","quasipoisson")) 
    rr = FALSE
  if (is.null(intvals)) {
    intvals <- (seq_len(length(table(qdata[expnms[1]])))) - 
      1
  }
  predit <- function(idx, newdata) {
    newdata[, expnms] <- idx
    suppressWarnings(predict(fit, newdata = newdata, type = "response"))
  }
  if (MCsize == nrow(qdata)) {
    newdata <- qdata
  }
  else {
    newids <- data.frame(temp = sort(sample(unique(qdata[, 
                                                         id, drop = TRUE]), MCsize, replace = TRUE)))
    names(newids) <- id
    newdata <- merge(qdata, newids, by = id, all.x = FALSE, 
                     all.y = TRUE)[seq_len(MCsize), ]
  }
  predmat = lapply(intvals, predit, newdata = newdata)
  msmdat <- data.frame(cbind(Ya = unlist(predmat), psi = rep(intvals, 
                                                             each = MCsize), weights = rep(newdata$weights, times = length(intvals))))
  msmforms = paste0("Ya ~ ", ifelse(hasintercept, "1 +", "-1 +"), 
                    "poly(psi, degree=", degree, ", raw=TRUE)")
  msmform = as.formula(msmforms)
  if (bayes) {
    if (!rr) 
      suppressWarnings(msmfit <- bayesglm(msmform, data = msmdat, 
                                          weights = weights, x = TRUE, ...))
    if (rr) 
      suppressWarnings(msmfit <- bayesglm(msmform, data = msmdat, 
                                          family = binomial(link = "log"), start = rep(-1e-04, 
                                                                                       degree + 1), weights = weights, x = TRUE))
  }
  if (!bayes) {
    if (!rr) 
      suppressWarnings(msmfit <- glm(msmform, data = msmdat, 
                                     weights = weights, x = TRUE, ...))
    if (rr) 
      suppressWarnings(msmfit <- glm(msmform, data = msmdat, 
                                     family = binomial(link = "log"), start = rep(-1e-04, 
                                                                                  degree + 1), weights = weights, x = TRUE))
  }
  res <- list(fit = fit, msmfit = msmfit)
  if (main) {
    res$Ya <- msmdat$Ya
    res$Yamsm <- predict(msmfit, type = "response")
    res$Yamsml <- predict(msmfit, type = "link")
    res$A <- msmdat$psi
  }
  res
}
#### .intmaker ####
copy.intmaker <- function (f, expnms, emmvars) 
{
  rightside = as.character(f)[3]
  trms = strsplit(gsub(" ", "", rightside), "+", fixed = TRUE)[[1]]
  expidx <- which(trms %in% expnms)
  newtrmsl = lapply(emmvars, function(x) paste(paste0(trms[expidx], 
                                                      ":", x), collapse = "+"))
  newtrms = paste0(newtrmsl, collapse = "+")
  newrightside = paste(rightside, "+", newtrms)
  newf <- as.formula(paste0(as.character(f)[2], as.character(f)[1], 
                            newrightside))
  newf
}
#### vc_multiscomb ####
copy.vc_multiscomb <- function (inames = c("(Intercept)"), emmvars, expnms, addedintsl, 
          covmat, grad = NULL) 
{
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  expidx = ifelse(is.null(inames), 1, 2)
  dimnew <- expidx + length(emmvars) + length(addedintsl) -1
  dimold <- dim(covmat)[1]
  if (!is.null(grad[1])) 
    grad = NULL
  if (is.null(grad[1])) 
    grad <- 1
  nms = list(expnms)
  if (!is.null(inames)) 
    nms = c(inames, nms)
  for (i in seq_len(length(emmvars))) {
    nms = c(nms, addedintsl[i])
  }
  weightvec <- list()
  for (j in seq_len(dimnew)) {
    weightvec[[j]] = rep(0, dimold)
    vars = nms[[j]]
    if (j == expidx) {
      weightvec[[j]][which(colnames(covmat) %in% vars)] <- grad
    }
    else {
      weightvec[[j]][which(colnames(covmat) %in% vars)] <- 1
    }
  }
  outcov = matrix(NA, nrow = dimnew, ncol = dimnew)
  for (jj in seq_len(dimnew)) {
    for (ii in jj:dimnew) {
      outcov[jj, ii] <- outcov[ii, jj] <- weightvec[[jj]] %*% 
        covmat %*% weightvec[[ii]]
    }
  }
  outcov
}


#### qgcomp.emm.noboot ####
qgcomp.modi.noboot <- function (f, data, expnms = NULL, emmvar = NULL, q = 4, breaks = NULL, 
          id = NULL, weights, alpha = 0.05, bayes = FALSE, errcheck = TRUE, 
          ...) 
{
  if (errcheck) {
    if (is.null(expnms)) {
      stop("'expnms' must be specified explicitly\n")
    }
    if (is.null(emmvar)) {
      stop("'emmvar' must be specified explicitly\n")
    }
  }
  allemmvals <- unique(data[, emmvar])
  emmlev <- length(allemmvals)
  zdata = qgcompint:::zproc(data[, emmvar], znm = emmvar)
  emmvars = names(zdata)
  data = cbind(data, zdata)
  if (errcheck) {
  }
  originalform <- terms(f, data = data)
  hasintercept = as.logical(attr(originalform, "intercept"))
  (f <- copy.intmaker(f, expnms, emmvars))
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(originalform, 
                                                           "term.labels"))
  addedmain <- setdiff(addedterms, grep(":", addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, 
                                                 value = TRUE))
  addedintsord = addedints
  nobs = nrow(data)
  data$weights = rep(1, nobs)
  lin = qgcompint:::.intchecknames(expnms)
  if (!lin) 
    stop("Model appears to be non-linear: this is not yet implemented")
  if (!is.null(q) | !is.null(breaks)) {
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  } else {
    qdata <- data
    br <- breaks
  }
  if (is.null(id)) {
    id = "id__"
    qdata$id__ = seq_len(dim(qdata)[1])
  }
  if (!bayes) 
    fit <- glm(newform, data = qdata, weights = weights, 
               ...)
  if (bayes) {
    requireNamespace("arm")
    fit <- arm::bayesglm(newform, data = qdata, weights = weights, 
                         ...)
  }
  mod <- summary(fit)
  if (length(setdiff(expnms, rownames(mod$coefficients))) > 0) {
    stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures.\n           Try one of the following:\n             1) set 'bayes' to TRUE in the qgcomp function (recommended)\n             2) set 'q' to a higher value in the qgcomp function (recommended)\n             3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)\n           ")
  }
  estb <- sum(mod$coefficients[expnms, 1, drop = TRUE])
  seb <- qgcompint:::se_comb2(expnms, covmat = mod$cov.scaled)
  if (hasintercept) {
    estb <- c(fit$coefficients[1], estb)
    seb <- c(sqrt(mod$cov.scaled[1, 1]), seb)
  }
  tstat <- estb/seb
  df <- mod$df.null - length(expnms) - length(addedints) - 1
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - alpha/2))
  estb.prod <- do.call(c, lapply(1:length(emmvars), function(x) sum(mod$coefficients[addedintsl[[x]], 1, drop = TRUE])))
  names(estb.prod) <- do.call(c, lapply(1:length(emmvars), 
                                        function(x) paste0(emmvars[x], ":mixture")))
  seb.prod <- do.call(c, lapply(1:length(emmvars), function(x) qgcompint:::se_comb2(addedintsl[[x]], covmat = mod$cov.scaled)))
  tstat.prod <- estb.prod/seb.prod
  pval.prod <- 2 - 2 * pt(abs(tstat.prod), df = df)
  pvalz.prod <- 2 - 2 * pnorm(abs(tstat.prod))
  ci.prod <- cbind(estb.prod + seb.prod * qnorm(alpha/2), estb.prod + 
                     seb.prod * qnorm(1 - alpha/2))
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  pos.coef <- which(wcoef > 0)
  neg.coef <- which(wcoef <= 0)
  pos.weights <- abs(wcoef[pos.coef])/sum(abs(wcoef[pos.coef]))
  neg.weights <- abs(wcoef[neg.coef])/sum(abs(wcoef[neg.coef]))
  pos.psi <- sum(wcoef[pos.coef])
  neg.psi <- sum(wcoef[neg.coef])
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  psiidx = 1 + hasintercept
  cnms = "psi1"
  inames = NULL
  if (hasintercept) {
    inames = "(Intercept)"
    cnms = c(inames, cnms)
  }
  covmat.coef = copy.vc_multiscomb(inames = inames, emmvars = emmvars, 
                              expnms = expnms, addedintsl = addedintsl, covmat = mod$cov.scaled, 
                              grad = NULL)
  names(estb) <- cnms
  colnames(covmat.coef) <- rownames(covmat.coef) <- names(c(estb, 
                                                            estb.prod))
  res <- qgcompint:::.qgcompemm_object(qx = qx, fit = fit, psi = estb[psiidx], 
                           psiint = estb.prod[2 * (1:length(emmvars))], var.psi = seb[psiidx]^2, 
                           var.psiint = seb.prod[2 * (1:length(emmvars))]^2, covmat.psi = covmat.coef["psi1", 
                                                                                                      "psi1"], covmat.psiint = covmat.coef[grep("mixture", 
                                                                                                                                                colnames(covmat.coef)), grep("mixture", colnames(covmat.coef))], 
                           ci = ci[psiidx, ], ciint = ci.prod[(1:length(emmvars)),], coef = c(estb, estb.prod), var.coef = c(seb^2, seb.prod^2), covmat.coef = covmat.coef, 
                           ci.coef = rbind(ci, ci.prod), expnms = expnms, intterms = addedintsord, 
                           q = q, breaks = br, degree = 1, pos.psi = pos.psi, neg.psi = neg.psi, 
                           pos.weights = sort(pos.weights, decreasing = TRUE), neg.weights = sort(neg.weights, 
                                                                                                  decreasing = TRUE), pos.size = sum(abs(wcoef[pos.coef])), 
                           neg.size = sum(abs(wcoef[neg.coef])), bootstrap = FALSE, 
                           cov.yhat = NULL, alpha = alpha, emmlev = emmlev, mod.with.se = mod)
  if (emmlev == 2) {
    ww = getstratweights(res, emmval = 1)
    ff = getstrateffects(res, emmval = 1)
    cl = class(res)
    res = c(res, list(pos.weights1 = ww$pos.weights, neg.weights1 = ww$neg.weights, 
                      pos.psi1 = ww$pos.psi, neg.psi1 = ww$neg.psi, pos.size1 = ww$pos.size, 
                      neg.size1 = ww$neg.size), list(effect = ff$eff, vareffect = ff$se^2, 
                                                     cieffect = ff$ci))
    class(res) = cl
  }
  if (fit$family$family == "gaussian") {
    res$tstat <- c(tstat, tstat.prod)
    res$df <- df
    res$pval <- c(pval, pval.prod)
  }
  if (fit$family$family %in% c("binomial", "poisson")) {
    res$zstat <- c(tstat, tstat.prod)
    res$pval <- c(pvalz, pvalz.prod)
  }
  res
}


#### qgcomp.modi.boot ####
qgcomp.modi.boot <- function (f, data, expnms = NULL, emmvar = "", q = 4, breaks = NULL, 
                              id = NULL, weights, alpha = 0.05, B = 200, rr = TRUE, degree = 1, 
                              seed = NULL, bayes = FALSE, MCsize = nrow(data), parallel = FALSE, 
                              parplan = FALSE, errcheck = FALSE, ...) 
{
  oldq = NULL
  if (is.null(seed)) 
    seed = round(runif(1, min = 0, max = 1e+08))
  if (errcheck) {
    if (is.null(expnms)) {
      stop("'expnms' must be specified explicitly\n")
    }
    if (is.null(emmvar)) {
      stop("'emmvar' must be specified explicitly\n")
    }
  }
  allemmvals <- unique(data[, emmvar])
  emmlev <- length(allemmvals)
  zdata = qgcompint:::zproc(data[, emmvar], znm = emmvar)
  emmvars = names(zdata)
  data = cbind(data, zdata)
  data = data[, unique(names(data)), drop = FALSE]
  if (errcheck) {
  }
  originalform <- terms(f, data = data)
  hasintercept = as.logical(attr(originalform, "intercept"))
  (f <- copy.intmaker(f, expnms, emmvars))
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(originalform, 
                                                           "term.labels"))
  addedmain <- setdiff(addedterms, grep(":", addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, 
                                                 value = TRUE))
  addedintsord = addedints
  class(newform) <- "formula"
  nobs = nrow(data)
  data$weights = rep(1, nobs)
  if (is.null(expnms)) {
    expnms <- attr(newform, "term.labels")
    message("Including all model terms as exposures of interest\n")
  }
  lin = qgcompint:::.intchecknames(expnms)
  if (!lin) 
    stop("Model appears to be non-linear and I'm having trouble parsing it:\n                  please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)) {
    oldq = q
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)) {
    ql <- qgcomp::quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if (is.null(q)) {
      nvals <- length(br[[1]]) - 1
    } else {
      nvals <- q
    }
    intvals <- (seq_len(nvals)) - 1
  } else {
    qdata <- data[unique(names(data)), ]
    nvals = length(table(unlist(data[, expnms])))
    if (nvals < 10) {
      message("\nNote: using all possible values of exposure as the\n              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[, expnms]))))
      br <- lapply(seq_len(p), function(x) c(-1e+16, intvals[2:nvals] - 
                                               1e-16, 1e+16))
    } else {
      message("\nNote: using quantiles of all exposures combined in order to set\n          proposed intervention values for overall effect (25th, 50th, 75th %ile)\n        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
      intvals = as.numeric(quantile(unlist(data[, expnms]), 
                                    c(0.25, 0.5, 0.75)))
      br <- NULL
    }
  }
  if (is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  msmfit <- copy.msm.emm.fit(newform, qdata, intvals, emmvar = emmvar, weights = weights,
                        emmvars = emmvars, expnms = expnms, rr, main = T, 
                        degree = degree, id = id, bayes, MCsize = MCsize, 
                        ...)
  msmcoefnames <- attr(msmfit, "term.labels")
  estb <- as.numeric(msmfit$msmfit$coefficients)
  nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[, id, drop = TRUE]))
  starttime = Sys.time()
  psi.emm.only <- function(i = 1, f = f, qdata = qdata, intvals = intvals, 
                           emmvar = emmvar, emmvars = emmvars, expnms = expnms, 
                           rr = rr, degree = degree, nids = nids, id = id, weights, 
                           MCsize = MCsize, ...) {
    if (i == 2 & !parallel) {
      timeiter = as.numeric(Sys.time() - starttime)
      if ((timeiter * B/60) > 0.5) 
        message(paste0("Expected time to finish: ", round(B * 
                                                            timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp = sort(sample(unique(qdata[, 
                                                          id, drop = TRUE]), nids, replace = TRUE)))
    names(bootids) <- id
    qdata_ <- merge(qdata, bootids, by = id, all.x = FALSE, 
                    all.y = TRUE)
    ft = copy.msm.emm.fit(f, qdata_, intvals = intvals, expnms = expnms, weights = weights,
                     emmvar = emmvar, emmvars = emmvars, rr, main = FALSE, 
                     degree, id, bayes, MCsize = MCsize, 
                     ...)
    yhatty = data.frame(yhat = predict(ft$msmfit), psi = ft$msmfit$data[,"psi"])
    as.numeric(c(with(yhatty, tapply(yhat, psi, mean)), ft$msmfit$coefficients))
  }
  set.seed(seed)
  if (parallel) {
    if (parplan) {
      oplan <- future::plan(strategy = future::multisession,workers = 10)
      on.exit(future::plan(oplan), add = TRUE)
    }
    bootsamps <- future.apply::future_lapply(X = seq_len(B), 
                                             FUN = psi.emm.only, f = f, qdata = qdata, intvals = intvals, 
                                             emmvar = emmvar, emmvars = emmvars, expnms = expnms, 
                                             rr = rr, degree = degree, nids = nids, id = id, weights = qdata$weights, 
                                             MCsize = MCsize, future.seed = TRUE, ...)
  }
  else {
    bootsamps <- lapply(X = seq_len(B), FUN = psi.emm.only, 
                        f = f, qdata = qdata, intvals = intvals, emmvar = emmvar, 
                        emmvars = emmvars, expnms = expnms, rr = rr, degree = degree, 
                        nids = nids, id = id, weights = weights, MCsize = MCsize, 
                        ...)
  }
  bootsamps = do.call("cbind", bootsamps)
  hatidx = seq_len(length(intvals))
  hats = t(bootsamps[hatidx, ])
  cov.yhat = cov(hats)
  bootsamps = bootsamps[-hatidx, ]
  seb <- apply(bootsamps, 1, sd)
  covmat.coef <- cov(t(bootsamps))
  colnames(covmat.coef) <- rownames(covmat.coef) <- names(estb) <- c("(Intercept)", 
                                                                     msmcoefnames)
  tstat <- estb/seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 
    1 - degree
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha/2), estb + seb * qnorm(1 - 
                                                                alpha/2))
  if (!is.null(oldq)) {
    q = oldq
  }
  psidx = 1:(hasintercept + 1)
  qx <- qdata[, expnms]
  res <- qgcompint:::.qgcompemm_object(qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit, 
                           psi = estb[-1], var.psi = seb[-1]^2, covmat.psi = covmat.coef["psi1", 
                                                                                         "psi1"], covmat.psiint = covmat.coef[grep("mixture", 
                                                                                                                                   colnames(covmat.coef)), grep("mixture", colnames(covmat.coef))], 
                           ci = ci[-1, ], coef = estb, var.coef = seb^2, covmat.coef = covmat.coef, 
                           ci.coef = ci, expnms = expnms, intterms = addedintsord, 
                           q = q, breaks = br, degree = degree, pos.psi = NULL, 
                           neg.psi = NULL, pos.weights = NULL, neg.weights = NULL, 
                           pos.size = NULL, neg.size = NULL, bootstrap = TRUE, y.expected = msmfit$Ya, 
                           y.expectedmsm = msmfit$Yamsm, index = msmfit$A, emmvar.msm = msmfit[[emmvar]], 
                           bootsamps = bootsamps, cov.yhat = cov.yhat, alpha = alpha, emmlev = emmlev)
  if (msmfit$fit$family$family == "gaussian") {
    res$tstat <- tstat
    res$df <- df
    res$pval <- pval
  }
  if (msmfit$fit$family$family %in% c("binomial", "poisson")) {
    res$zstat <- tstat
    res$pval <- pvalz
  }
  res
}


#### msm.emm.fit ####
copy.msm.emm.fit <- function (f, qdata, intvals, expnms, emmvar, emmvars, rr = TRUE, 
                              main = T, degree = 1, id = NULL,  bayes = FALSE, weights = weights,
                              MCsize = nrow(qdata), hasintercept = TRUE, ...) 
{
  newform <- terms(f, data = qdata)
  nobs = nrow(qdata)
  qdata$weights = rep(1, nobs)
  if (is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  nidx = which(!(names(qdata) %in% id))
  if (!bayes) 
    fit <- glm(newform, data = qdata, weights = weights,
               ...)
  if (bayes) {
    requireNamespace("arm")
    fit <- bayesglm(f, data = qdata[, nidx, drop = FALSE], 
                    weights = weights, ...)
  }
  if (fit$family$family %in% c("gaussian", "poisson")) 
    rr = FALSE
  if (is.null(intvals)) {
    intvals <- (seq_len(length(table(qdata[expnms[1]])))) - 
      1
  }
  predit <- function(idx, newdata) {
    newdata[, expnms] <- idx
    suppressWarnings(predict(fit, newdata = newdata, type = "response"))
  }
  if (MCsize == nrow(qdata)) {
    newdata <- qdata
  } else {
    newids <- data.frame(temp = sort(sample(unique(qdata[, id, drop = TRUE]), MCsize, replace = TRUE)))
    names(newids) <- id
    newdata <- merge(qdata, newids, by = id, all.x = FALSE, all.y = TRUE)[seq_len(MCsize), ]
  }
  predmat = lapply(intvals, predit, newdata = newdata)
  msmdat <- data.frame(cbind(Ya = do.call(c, predmat), psi = rep(intvals, 
                                                                 each = MCsize), weights = rep(newdata$weights, times = length(intvals))))
  msmdat[, emmvars] <- newdata[, emmvars]
  polydat <- as.data.frame(poly(msmdat$psi, degree = degree, raw = TRUE))
  newexpnms <- paste0("psi", 1:degree)
  names(polydat) <- newexpnms
  msmdat <- cbind(msmdat, polydat)
  msmf <- paste0("Ya ~ ", ifelse(hasintercept, "1 +", "-1 +"), 
                 paste0(newexpnms, collapse = "+"))
  msmform <- copy.intmaker(as.formula(msmf), expnms = newexpnms, emmvars = emmvars)
  class(msmform) <- "formula"
  newterms <- terms(msmform)
  nterms = length(attr(newterms, "term.labels"))
  nterms = nterms + attr(newterms, "intercept")
  if (bayes) {
    if (!rr) 
      suppressWarnings(msmfit <- bayesglm(msmform, data = msmdat, 
                                          weights = weights, x = TRUE, ...))
    if (rr) 
      suppressWarnings(msmfit <- bayesglm(msmform, data = msmdat, 
                                          family = binomial(link = "log"), start = rep(-1e-04, 
                                                                                       nterms), weights = weights, x = TRUE))
  }
  if (!bayes) {
    if (!rr) 
      suppressWarnings(msmfit <- glm(msmform, data = msmdat, 
                                     weights = weights, x = TRUE, ...))
    if (rr) 
      suppressWarnings(msmfit <- glm(msmform, data = msmdat, 
                                     family = binomial(link = "log"), start = rep(-1e-04, 
                                                                                  nterms), weights = weights, x = TRUE))
  }
  res <- list(fit = fit, msmfit = msmfit)
  if (main) {
    res$Ya <- msmdat$Ya
    res$Yamsm <- as.numeric(predict(msmfit, type = "response"))
    res$Yamsml <- as.numeric(predict(msmfit, type = "link"))
    res$A <- msmdat$psi
    res[[emmvar]] <- do.call(c, lapply(intvals, function(x) newdata[, 
                                                                    emmvar]))
  }
  newtermlabels <- attr(newterms, "term.labels")
  for (emmv in emmvars) {
    newtermlabels <- gsub(paste0("psi([0-9]):", emmv), paste0(emmv, 
                                                              ":", "mixture", "^\\1"), newtermlabels)
  }
  newtermlabels <- gsub("\\^1", "", newtermlabels)
  attr(res, "term.labels") <- newtermlabels
  res
}

#### gWQS with repeated holdhout validation ####
copy.gwqs = function (formula, data, na.action, weights, mix_name, stratified, 
                      valid_var, b = 100, b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE, 
                      q = 4, validation = 0.6, family = gaussian, signal = c("t2", 
                                                                             "one", "abst", "expt"), rs = FALSE, n_vars = NULL, zilink = c("logit", 
                                                                                                                                           "probit", "cloglog", "cauchit", "log"), seed = NULL, 
                      plan_strategy = "sequential", optim.method = c("BFGS", "Nelder-Mead", 
                                                                     "CG", "SANN"), control = list(trace = FALSE, maxit = 2000, 
                                                                                                   reltol = 1e-09), ...) 
{
  wqsGaussBin <- function(initp, kw, bXm, bY, boffset, bQ, 
                          kx, Xnames, n_levels, level_names, wqsvars, family, zilink, 
                          zero_infl, formula, ff, bwghts, stratified, b1_pos, b1_constr) {
    if (b1_constr) 
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = sum(family$dev.resids(y = bY, mu = family$linkinv(term), 
                              wt = bwghts))
    return(f)
  }
  wqsPoisson <- function(initp, kw, bXm, bY, boffset, bQ, kx, 
                         Xnames, n_levels, level_names, wqsvars, family, zilink, 
                         zero_infl, formula, ff, bwghts, stratified, b1_pos, b1_constr) {
    if (b1_constr) 
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum(dpois(bY, lambda = exp(term), log = TRUE))
    return(f)
  }
  wqsNegBin <- function(initp, kw, bXm, bY, boffset, bQ, kx, 
                        Xnames, n_levels, level_names, wqsvars, family, zilink, 
                        zero_infl, formula, ff, bwghts, stratified, b1_pos, b1_constr) {
    if (b1_constr) 
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum((suppressWarnings(dnbinom(bY, size = theta, 
                                       mu = exp(term), log = TRUE))) * bwghts)
    return(f)
  }
  wqsMultinom <- function(initp, kw, bXm, bY, boffset, bQ, 
                          kx, Xnames, n_levels, level_names, wqsvars, family, zilink, 
                          zero_infl, formula, ff, bwghts, stratified, b1_pos, b1_constr) {
    if (b1_constr) {
      par_pos <- which(grepl("wqs", names(initp)))
      initp[par_pos] <- sapply(1:length(b1_pos), function(i) ifelse(b1_pos[i], 
                                                                    abs(initp[par_pos[i]]), -abs(initp[par_pos[i]])))
    }
    w <- matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels - 
                  1)
    w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
    bXm[, wqsvars] <- bQ %*% w
    b_covs = matrix(0, kx, n_levels - 1)
    i = 1:(kx/(n_levels - 1))
    for (j in 0:(n_levels - 2)) {
      b_covs[(kx/(n_levels - 1)) * j + i, j + 1] = initp[(kx/(n_levels - 
                                                                1)) * j + i]
    }
    term = bXm %*% b_covs + boffset
    f = -sum((diag(bY %*% t(term)) - log(1 + rowSums(exp(term)))) * 
               bwghts)
    return(f)
  }
  one <- function(x) {
    rep(1, length(x))
  }
  t2 <- function(x) {
    x^2
  }
  if (is.character(family)) {
    if (family %in% c("multinomial", "negbin")) 
      family <- list(family = family)
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
  }
  objfn <- switch(family$family, gaussian = wqsGaussBin, binomial = wqsGaussBin, 
                  poisson = wqsPoisson, quasipoisson = wqsPoisson, negbin = wqsNegBin, 
                  multinomial = wqsMultinom)
  optim.method <- match.arg(optim.method)
  signal <- match.arg(signal)
  signal <- switch(signal, one = one, abst = abs, t2 = t2, 
                   expt = exp)
  if (!any(grepl("wqs", rownames(attr(terms(formula), "factors"))))) 
    stop("'formula' must contain 'wqs' term: e.g. y ~ wqs + ...\n")
  if (zero_infl) {
    zilink <- make.link(match.arg(zilink))
    ff = formula
    if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
                                              as.name("|"))) 
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", 
               "weights", "stratified", "valid_var"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("selectdatavars")
  dtf <- eval(dtf, parent.frame())
  l <- list(...)
  solve_dir_issue <- ifelse(is.null(l$solve_dir_issue), FALSE, 
                            l$solve_dir_issue)
  if (is.null(q)) {
    Q = as.matrix(dtf[, mix_name])
    qi = NULL
  }
  else {
    q_f = gWQS:::gwqs_rank(dtf, mix_name, q)
    Q = q_f$Q
    qi = q_f$qi
  }
  m <- match(c("stratified"), names(mc), 0)
  if (m) {
    strtfd_out = gWQS:::stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name = strtfd_out$mix_name
  }
  else stratified <- NULL
  if (!is.null(l$stratified_rh)) 
    stratified <- l$stratified_rh
  if (!is.null(seed) & !is.numeric(seed)) 
    stop("seed must be numeric or NULL\n")
  if (!is.null(seed)) {
    set.seed(seed)
    fseed <- TRUE
  }
  else fseed <- FALSE
  N <- nrow(dtf)
  m <- match(c("valid_var"), names(mc), 0)
  rindex = gWQS:::create_rindex(dtf, N, validation, valid_var, m, 
                         family)
  if (is.null(n_vars)) 
    n_vars = round(sqrt(length(mix_name)))
  if (family$family %in% c("gaussian", "quasipoisson")) 
    ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", 
                                "negbin")) 
    ts = "z"
  if (!is.numeric(b)) 
    stop("'b' must be a number\n")
  m <- match(c("weights"), names(mc), 0)
  if (m[1]) 
    dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else dtf$wghts <- wghts <- rep(1, N)
  Xnames <- gWQS:::parnames(dtf, formula, NULL)
  kx <- length(Xnames)
  if (family$family == "multinomial") {
    n_levels <- nlevels(eval(formula[[2]], envir = dtf))
    if (n_levels == 0) 
      stop("y must be of class factor when 'family = \"multinomial\"'\n")
    level_names <- levels(eval(formula[[2]], envir = dtf))
    Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:kx], 
                                                           "_", i, "_vs_", level_names[1])))
    kx <- kx * (n_levels - 1)
    wqsvars = Xnames[grepl("^wqs_", Xnames)]
    dtf[, wqsvars] <- 0
  }
  else {
    n_levels <- ifelse(is.null(stratified), 2, nlevels(unlist(dtf[, 
                                                                  stratified])))
    level_names <- wqsvars <- NULL
  }
  kw <- ncol(Q)
  initp <- gWQS:::values.0(kw, Xnames, kx, n_levels, formula, ff, 
                    wghts, dtf, stratified, b1_pos, family, zilink, zero_infl)
  mf <- model.frame(formula, dtf)
  Y <- model.response(mf, "any")
  if (family$family == "binomial" & any(class(Y) %in% c("factor", 
                                                        "character"))) {
    if (class(Y) == "character") 
      Y = factor(Y)
    Y <- as.numeric(Y != levels(Y)[1])
  }
  if (family$family == "multinomial") 
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == 
                                                       level_names[i], 1, 0)))
  offset <- model.offset(mf)
  if (is.null(offset)) 
    offset <- rep(0, nrow(dtf))
  if (family$family == "multinomial") {
    Xl = lapply(wqsvars, function(i) {
      fmi <- as.formula(gsub("wqs", i, format(formula)))
      model.matrix(fmi, data = dtf)
    })
    Xm = do.call("cbind", Xl)
  }
  else Xm <- model.matrix(formula, dtf)
  future:::plan(plan_strategy, workers = 4)
  if (control$trace) 
    cat("start opt\n")
  param <- future.apply:::future_lapply(X = 1:b, FUN = gWQS:::optim.f, objfn = objfn, 
                         Y = if (family$family == "multinomial") 
                           Y[rindex$it, ]
                         else Y[rindex$it], Xm = Xm[rindex$it, ], Q = Q[rindex$it, 
                         ], offset = offset[rindex$it], wghts = wghts[rindex$it], 
                         initp = initp, n_levels = n_levels, level_names = level_names, 
                         wqsvars = wqsvars, b1_pos = b1_pos, b1_constr = b1_constr, 
                         n_vars = n_vars, family = family, rs = rs, zilink = zilink, 
                         zero_infl = zero_infl, formula = formula, ff = ff, kx = kx, 
                         kw = kw, Xnames = Xnames, stratified = stratified, b = b, 
                         optim.method = optim.method, control = control, future.seed = fseed)
  conv <- c(sapply(param, function(i) i$conv))
  counts <- c(sapply(param, function(i) i$counts))
  val <- c(sapply(param, function(i) i$val))
  mex <- lapply(param, function(i) i$mex)
  bindex <- lapply(param, function(i) i$bindex)
  slctd_vars <- lapply(param, function(i) i$slctd_vars)
  if (rs) {
    future:::plan(plan_strategy, workers = 4)
    param <- future.apply:::future_lapply(X = 1:b, FUN = set_par_names, 
                           slctd_vars, param, q_name = colnames(Q), family = family, 
                           future.seed = FALSE)
  }
  if (family$family == "multinomial") {
    n_levels <- dim(param[[1]]$par_opt)[2] + 1
    wqs_site <- which(grepl("^wqs_", rownames(param[[1]]$mfit$m_f$coefficients)))
    wght_matrix <- lapply(1:(n_levels - 1), function(j) do.call("rbind", 
                                                                lapply(param, function(i) i$par_opt[, j])))
    b1 <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$nlm_out$estimate[j]))
    se <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$Standard_Error[j]))
    stat <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$stat[j]))
    p_val <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$p_value[j]))
  }
  else {
    wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))
    b1 <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", 
                                                                     "Estimate"])
    se <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", 
                                                                     "Std. Error"])
    stat <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", 
                                                                       paste0(ts, " value")])
    p_val <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", 
                                                                        gsub("x", ts, "Pr(>|x|)")])
    n_levels <- 1
  }
  n_non_conv = sum(conv == 1)
  if (n_non_conv == 0 & control$trace) 
    cat(paste0("The optimization function always converged\n"))
  else if (n_non_conv == b) 
    stop("The optimization function never converged\n")
  else if (control$trace) 
    cat(paste0("The optimization function did not converge ", 
               n_non_conv, " time/times\n"))
  if (control$trace) 
    cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 0, 
                                                na.rm = T), sum(b1 <= 0, na.rm = T)), ifelse(b1_pos, 
                                                                                             " positive", " negative"), " bootstrapped b1 out of ", 
               b, "\n"))
  if (family$family == "multinomial") {
    bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
    bres <- lapply(bres, as.data.frame)
    bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", 
                                     "stat", "p_val"))
    strata_names <- gsub("wqs_", "", rownames(param[[1]]$mfit$m_f$coefficients)[wqs_site])
    names(bres) <- strata_names
  }
  else {
    bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, 
                                p_val))
    names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", 
                     "p_val")
    strata_names <- NULL
  }
  if (family$family == "multinomial") {
    mean_weight <- lapply(1:(n_levels - 1), function(i) {
      if (rs) 
        bres[[i]][mix_name][is.na(bres[[i]][mix_name])] <- 0
      if (b1_pos[i]) 
        w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv == 
                                0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 > 
                                                                                   0 & conv == 0, "stat"]))
      else if (!b1_pos[i]) 
        w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv == 
                                0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 < 
                                                                                   0 & conv == 0, "stat"]))
      if (all(is.nan(w_t))) {
        if (solve_dir_issue == "inverse") 
          w_t <- 1/apply(bres[[i]][, mix_name], 2, weighted.mean, 
                         signal(bres[[i]][, "stat"]))/sum(1/apply(bres[[i]][, 
                                                                            mix_name], 2, weighted.mean, signal(bres[[i]][, 
                                                                                                                          "stat"])))
        else if (solve_dir_issue == "average") 
          w_t <- rep(1/length(mix_name), length(mix_name))
        if (!(solve_dir_issue %in% c("average", "inverse"))) 
          stop(paste0("There are no ", ifelse(b1_pos[i], 
                                              "positive", "negative"), " b1 in the bootstrapped models for ", 
                      strata_names[i], "\n"))
      }
      return(w_t)
    })
    mean_weight <- list.cbind(mean_weight)
  }
  else {
    if (rs) 
      bres[mix_name][is.na(bres[mix_name])] <- 0
    if (b1_pos) 
      mean_weight = apply(bres[bres$b1 > 0 & conv == 0, 
                               mix_name], 2, weighted.mean, signal(bres[bres$b1 > 
                                                                          0 & conv == 0, "stat"]))
    else mean_weight = apply(bres[bres$b1 < 0 & conv == 0, 
                                  mix_name], 2, weighted.mean, signal(bres[bres$b1 < 
                                                                             0 & conv == 0, "stat"]))
    if (all(is.nan(mean_weight))) {
      if (solve_dir_issue == "inverse") 
        mean_weight <- 1/apply(bres[, mix_name], 2, weighted.mean, 
                               signal(bres[, "stat"]))/sum(1/apply(bres[, 
                                                                        mix_name], 2, weighted.mean, signal(bres[, 
                                                                                                                 "stat"])))
      else if (solve_dir_issue == "average") 
        mean_weight <- rep(1/length(mix_name), length(mix_name))
      if (!(solve_dir_issue %in% c("average", "inverse"))) 
        stop("There are no ", ifelse(b1_pos, "positive", 
                                     "negative"), " b1 in the bootstrapped models\n")
    }
  }
  wqs_model = model.fit(mean_weight, dtf[rindex$iv, ], Q[rindex$iv, 
  ], if (family$family == "multinomial") 
    Y[rindex$iv, ]
  else Y[rindex$iv], family, zilink, formula, ff, wghts[rindex$iv], 
  offset[rindex$iv], initp, Xnames, n_levels, level_names, 
  wqsvars, stratified, b1_pos, zero_infl, kx, kw)
  if (all(grepl("wqs", attr(terms(formula), "term.labels")))) 
    y_plot <- model.response(model.frame(formula, dtf[rindex$iv, 
    ]), "any")
  else {
    formula_wo_wqs = remove_terms(formula, "wqs")
    if (family$family != "multinomial") {
      if (zero_infl) {
        if (length(ff[[3]]) > 1 && identical(ff[[3]][[1]], 
                                             as.name("|"))) {
          if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                            " ~ ", deparse(ff[[3]][[2]])))), "term.labels")))) 
            f1 <- as.formula(paste0(ff[[2]], " ~ ", 1))
          else f1 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                    " ~ ", deparse(ff[[3]][[2]]))), "wqs")
          if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                            " ~ ", deparse(ff[[3]][[3]])))), "term.labels")))) 
            f2 <- as.formula(paste0(ff[[2]], " ~ ", 1))
          else f2 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                    " ~ ", deparse(ff[[3]][[3]]))), "wqs")
          formula_wo_wqs <- as.formula(paste0(deparse(f1), 
                                              " | ", f2[[3]]))
        }
        fit <- zeroinfl(formula_wo_wqs, dtf[rindex$iv, 
        ], dist = family$family, link = zilink$name)
      }
      else {
        if (family$family == "negbin") 
          fit = glm.nb(formula_wo_wqs, dtf[rindex$iv, 
          ])
        else fit = glm(formula_wo_wqs, dtf[rindex$iv, 
        ], family = family)
      }
      if (family$family == "binomial") 
        y_plot = fit$y
      else y_plot = mean(fit$y) + resid(fit, type = "pearson")
    }
  }
  data_plot <- data.frame(mix_name, mean_weight, stringsAsFactors = TRUE)
  if (family$family == "multinomial") {
    Y <- model.response(model.frame(formula, dtf[rindex$iv, 
    ]), "any")
    level_names <- levels(Y)
    names(data_plot)[2:n_levels] = strata_names
    data_plot <- data_plot[order(data_plot[, strata_names[1]], 
                                 decreasing = TRUE), ]
    wqs_index <- wqs_model$wqs
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == 
                                                       level_names[i], 1, 0)))
    colnames(Y) <- strata_names
    y_adj_wqs_df <- do.call("rbind", lapply(strata_names, 
                                            function(i) {
                                              data.frame(level = i, y = Y[rowSums(Y[, -which(colnames(Y) == 
                                                                                               i), drop = F]) == 0, i], wqs = wqs_model$wqs[rowSums(Y[, 
                                                                                                                                                      -which(colnames(wqs_model$wqs) == i), drop = F]) == 
                                                                                                                                              0, i], stringsAsFactors = TRUE)
                                            }))
    dtf <- cbind(dtf, Q %*% mean_weight)
    names(dtf)[(ncol(dtf) - ncol(mean_weight) + 1):ncol(dtf)] <- colnames(wqs_index)
    dtf$wqs <- NULL
  }
  else {
    y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs_model$wqs))
    names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", 
                                    "y", "y_adj"), "wqs")
    data_plot <- data_plot[order(data_plot$mean_weight, decreasing = TRUE), 
    ]
    wqs_index <- as.numeric(unlist(wqs_model$wqs))
    dtf$wqs <- as.numeric(Q %*% mean_weight)
  }
  results = list(fit = wqs_model$m_f, final_weights = data_plot, 
                 conv = conv, bres = bres, wqs = wqs_index, qi = qi, bindex = bindex, 
                 slctd_vars = slctd_vars, tindex = rindex$it, vindex = rindex$iv, 
                 y_wqs_df = y_adj_wqs_df, family = family, call = cl, 
                 formula = formula, mix_name = mix_name, stratified = stratified, 
                 q = q, n_levels = n_levels, zero_infl = zero_infl, zilink = zilink, 
                 levelnames = strata_names, data = dtf, objfn_values = val, 
                 optim_messages = mex)
  if (zero_infl) 
    results$formula <- ff
  class(results) <- "gwqs"
  return(results)
}

