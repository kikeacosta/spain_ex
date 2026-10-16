rm(list=ls()); gc()
library(Epi)
library(MASS)
library(splines)
library(tidyverse)
library(MortalitySmooth)
library(viridis)

select <- dplyr::select

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")

# amin=20; amax=90; pmin=1970; pmax=2019; gr=1; sx = "female"; cd = "ESP"

### Most of this function is an adaptation of the apc.fit function from the Epi package (Carstensen 2019, https://CRAN.R-project.org/package=Epi)
curvature_dAPC <- function(amin=0, amax=60, pmin=1990, pmax=2016, gr=1, sx = "female", cd = "ESP"){
  
  amax <- amin - 1 + (floor((amax-amin+1)/gr))*gr
  pmax <- min(pmax, 2019)
  pmax <- pmin - 1 + (floor((pmax-pmin+1)/gr))*gr
  
  data0 <- 
    hmd %>% 
    filter(year %in% pmin:pmax, age %in% amin:amax, sex == sx, code == cd) %>% 
    select(age, year, e = exposure, d = dts)
  
  data1 <- 
    data0 %>% 
    mutate(A = as.integer(amin+floor((age-amin)/gr)*gr), 
           P = as.integer(pmin+floor((year-pmin)/gr)*gr))
  
  data2 <- data1 %>% 
    reframe(D = sum(d), 
            Y = sum(e),
            .by = c(A, P)) 
  
  data <- 
    data2 %>% 
    mutate(C = as.integer(P - A),
           D = round(D),
           D = ifelse(D < 1, 1, D),
           D = as.integer(D),
           Y = as.integer(Y))
  
  # rm(list=setdiff(ls(), "data"))
  # select <- dplyr::select
  
  # ref.c <- NA
  # ref.p <- NA
  dist <- "poisson" 
  model<- "factor" 
  dr.extr <- "weighted"
  parm = "APC"
  npar=c(15,10,15)
  scale = 1
  alpha = 0.05 
  print.AOV = TRUE
  
  
  # dist <- match.arg(dist)
  # model <- match.arg(model)
  drtyp <- deparse(substitute(dr.extr))
  # parm <- toupper(match.arg(parm))
  if (!missing(data)) {
    if (length(match(c("A", "P", "D", "Y"), names(data))) != 
        4) 
      stop("Data frame ", deparse(substitute(data)), " has columns:\n", 
           names(data), "\nmust have variables:\n", "A (age), P (period), D (cases) and Y (person-time)")
    data <- data[, c("A", "P", "D", "Y")]
    data <- data[complete.cases(data), ]
    A <- data$A
    P <- data$P
    D <- data$D
    Y <- data$Y
  } else {
    nm <- logical(4)
    nm[1] <- missing(A)
    nm[2] <- missing(P)
    nm[3] <- missing(D)
    nm[4] <- missing(Y)
    if (any(nm)) 
      stop("Variable", if (sum(nm) > 1) 
        "s", paste(c(" A", " P", " D", " Y")[nm], collapse = ","), 
        " missing from input")
    if (diff(range(lv <- sapply(list(A = A, P = P, D = D, 
                                     Y = Y), length))) != 0) 
      stop("\nLengths of variables (", paste(paste(names(lv), 
                                                   lv, sep = ":"), collapse = ", "), ") are not the same.")
  }
  med <- function(x, y) {
    o <- order(x)
    a <- y[o]
    names(a) <- x[o]
    return(as.numeric(names(a[cumsum(a)/sum(a) > 0.5][1])))
  }
  p0 <- med(P, D)
  c0 <- med(P - A, D)
  ref.p <- F
  ref.c <- F
  if (is.list(npar) & length(npar) < 3) 
    stop("npar as a list should have length 3! \n")
  if (!is.list(npar) & length(npar) != 3) {
    npar <- rep(npar, 3)[1:3]
    names(npar) = c("A", "P", "C")
    cat("NOTE: npar is specified as:")
    print(npar)
  }
  if (is.null(names(npar))) 
    names(npar) <- c("A", "P", "C")
  lu <- paste(formatC(c(alpha/2, 1 - alpha/2) * 100, format = "f", 
                      digits = 1), "%", sep = "")
  proj.ip <- function(X, M, orth = FALSE, weight = rep(1, nrow(X))) {
    if (nrow(X) != length(weight)) 
      stop("Dimension of space and length of i.p. weights differ!")
    if (nrow(X) != nrow(M)) 
      stop("Dimension of space and rownumber of model matrix differ!")
    Pp <- solve(crossprod(X * sqrt(weight)), t(X * weight)) %*% 
      M
    PM <- X %*% Pp
    if (orth) 
      PM <- M - PM
    else PM
  }
  Thin.col <- function(X, tol = 0.000001) {
    QR <- qr(X, tol = tol, LAPACK = FALSE)
    X[, QR$pivot[seq(length = QR$rank)], drop = FALSE]
  }
  detrend <- function(M, t, weight = rep(1, nrow(M))) {
    Thin.col(proj.ip(cbind(1, t), M, orth = TRUE, weight = weight))
  }
  if (is.list(model)) {
    if (!all(sapply(model, is.function))) 
      stop("'model' is a list, but not all elements are functions as they should be.")
    if ((lmod <- length(model)) < 3) 
      stop("'model' is a list, with", lmod, "elements, it should have three.")
    if (is.null(names(model))) 
      names(model) <- c("A", "P", "C")
    MA <- model[["A"]](A)
    MP <- model[["P"]](P)
    MC <- model[["C"]](P - A)
    Rp <- model[["P"]](p0)
    Rc <- model[["C"]](c0)
  } else {
    if (model == "factor") {
      MA <- model.matrix(~factor(A) - 1)
      MP <- model.matrix(~factor(P) - 1)
      MC <- model.matrix(~factor(P - A) - 1)
      Rp <- MP[abs(P - p0) == min(abs(P - p0)), , drop = FALSE][1, 
      ]
      Rc <- MC[abs(P - A - c0) == min(abs(P - A - c0)), 
               , drop = FALSE][1, ]
    }
    if (model == "ns") {
      knl <- is.list(npar)
      if (!knl & length(npar) == 1) 
        npar <- rep(npar, 3)
      if (is.null(names(npar))) 
        names(npar) <- c("A", "P", "C")
      names(npar) <- toupper(substr(names(npar), 1, 1))
      MA <- if (knl) 
        Ns(A, knots = npar[["A"]])
      else Ns(A, knots = quantile(rep(A, D), probs = (1:npar["A"] - 
                                                        0.5)/npar["A"]))
      MP <- if (knl) 
        Ns(P, knots = npar[["P"]])
      else Ns(P, knots = quantile(rep(P, D), probs = (1:npar["P"] - 
                                                        0.5)/npar["P"]))
      MC <- if (knl) 
        Ns(P - A, knots = npar[["C"]])
      else Ns(P - A, knots = quantile(rep(P - A, D), probs = (1:npar["C"] - 
                                                                0.5)/npar["C"]))
      Rp <- ns(p0, knots = attr(MP, "knots"), Boundary.knots = attr(MP, 
                                                                    "Boundary.knots"))
      Rc <- ns(c0, knots = attr(MC, "knots"), Boundary.knots = attr(MC, 
                                                                    "Boundary.knots"))
      Knots <- list(Age = sort(c(attr(MA, "knots"), attr(MA, 
                                                         "Boundary.knots"))), Per = sort(c(attr(MP, "knots"), 
                                                                                           attr(MP, "Boundary.knots"))), Coh = sort(c(attr(MC, 
                                                                                                                                           "knots"), attr(MC, "Boundary.knots"))))
    }
    if (model %in% c("bs", "ls")) {
      deg <- switch(model, ls = 1, bs = 3)
      knl <- is.list(npar)
      if (knl) 
        nk <- sapply(npar, length)
      MA <- if (knl) 
        bs(A, knots = npar[["A"]][-c(1, nk[1])], Boundary.knots = npar[["A"]][c(1, 
                                                                                nk[1])], degree = deg)
      else bs(A, df = npar[["A"]], degree = deg)
      MP <- if (knl) 
        bs(P, knots = npar[["P"]][-c(1, nk[2])], Boundary.knots = npar[["P"]][c(1, 
                                                                                nk[2])], degree = deg)
      else bs(P, df = npar[["P"]], degree = deg)
      MC <- if (knl) 
        bs(P - A, knots = npar[["C"]][-c(1, nk[3])], 
           Boundary.knots = npar[["C"]][c(1, nk[3])], 
           degree = deg)
      else bs(P - A, df = npar[["C"]], degree = deg)
      Rp <- bs(p0, knots = attr(MP, "knots"), Boundary.knots = attr(MP, 
                                                                    "Boundary.knots"), degree = attr(MP, "degree"))
      Rc <- bs(c0, knots = attr(MC, "knots"), Boundary.knots = attr(MC, 
                                                                    "Boundary.knots"), degree = attr(MC, "degree"))
      Knots <- list(Age = sort(c(attr(MA, "knots"), attr(MA, 
                                                         "Boundary.knots"))), Per = sort(c(attr(MP, "knots"), 
                                                                                           attr(MP, "Boundary.knots"))), Coh = sort(c(attr(MC, 
                                                                                                                                           "knots"), attr(MC, "Boundary.knots"))))
    }
  }
  if (tolower(substr(dist, 1, 2)) == "po") {
    m.APC <- glm(D ~ MA + I(P - p0) + MP + MC, offset = log(Y), 
                 family = poisson)
    Dist <- "Poisson with log(Y) offset"
  }
  if (tolower(substr(dist, 1, 3)) %in% c("bin")) {
    m.APC <- glm(cbind(D, Y - D) ~ MA + I(P - p0) + MP + 
                   MC, family = binomial)
    Dist <- "Binomial regression (logistic) of D/Y"
  }
  m.AP <- update(m.APC, . ~ . - MC)
  m.AC <- update(m.APC, . ~ . - MP)
  m.Ad <- update(m.AP, . ~ . - MP)
  m.A <- update(m.Ad, . ~ . - I(P - p0))
  m.0 <- update(m.A, . ~ . - MA)
  AOV <- anova(m.A, m.Ad, m.AC, m.APC, m.AP, m.Ad, test = "Chisq")
  attr(AOV, "heading") <- "\nAnalysis of deviance for Age-Period-Cohort model\n"
  attr(AOV, "row.names") <- c("Age", "Age-drift", "Age-Cohort", 
                              "Age-Period-Cohort", "Age-Period", "Age-drift")
  A.pt <- unique(A)
  A.pos <- match(A.pt, A)
  P.pt <- unique(P)
  P.pos <- match(P.pt, P)
  C.pt <- unique(P - A)
  C.pos <- match(C.pt, P - A)
  MA <- cbind(1, MA)
  if (!mode(dr.extr) %in% c("character", "numeric")) 
    stop("\"dr.extr\" must be of mode \"character\" or \"numeric\".\n")
  if (is.character(dr.extr)) {
    wt <- rep(1, length(D))
    drtyp <- "1-weights"
    if (toupper(substr(dr.extr, 1, 1)) %in% c("T", "D")) {
      wt <- D
      drtyp <- "D-weights"
    }
    else if (toupper(substr(dr.extr, 1, 1)) %in% c("L", "R")) {
      wt <- (Y^2)/D
      drtyp <- "Y^2/D-weights"
    }
    else if (toupper(substr(dr.extr, 1, 1)) %in% c("Y")) {
      wt <- Y
      drtyp <- "Y-weights"
    }
    else if (dr.extr %in% names(data)) {
      wt <- data[, dr.extr]
      drtyp <- paste(dr.extr, "weights")
    }
  }
  if (is.numeric(dr.extr)) 
    wt <- dr.extr
  Rp <- matrix(Rp, nrow = 1)
  Rc <- matrix(Rc, nrow = 1)
  xP <- detrend(rbind(Rp, MP), c(p0, P), weight = c(0, wt))
  xC <- detrend(rbind(Rc, MC), c(c0, P - A), weight = c(0, 
                                                        wt))
  MPr <- xP[-1, , drop = FALSE] - ref.p * xP[rep(1, nrow(MP)), 
                                             , drop = FALSE]
  MCr <- xC[-1, , drop = FALSE] - ref.c * xC[rep(1, nrow(MC)), 
                                             , drop = FALSE]
  if (length(grep("-", parm)) == 0) {
    if (parm %in% c("ADPC", "ADCP", "APC", "ACP")) 
      m.APC <- update(m.0, . ~ . - 1 + MA + I(P - p0) + 
                        MPr + MCr)
    drift <- rbind(ci.exp(m.APC, subset = "I\\(", alpha = alpha), 
                   ci.exp(m.Ad, subset = "I\\(", alpha = alpha))
    rownames(drift) <- c(paste("APC (", drtyp, ")", sep = ""), 
                         "A-d")
    if (parm == "ADCP") 
      m.APC <- update(m.0, . ~ . - 1 + MA + I(P - A - c0) + 
                        MPr + MCr)
    if (parm == "APC") {
      MPr <- cbind(P - p0, MPr)
      m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    if (parm == "ACP") {
      MCr <- cbind(P - A - c0, MCr)
      m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    Age <- cbind(Age = A.pt, ci.exp(m.APC, subset = "MA", 
                                    ctr.mat = MA[A.pos, , drop = FALSE], alpha = alpha))[order(A.pt), 
                                    ]
    Per <- cbind(Per = P.pt, ci.exp(m.APC, subset = "MPr", 
                                    ctr.mat = MPr[P.pos, , drop = FALSE], alpha = alpha))[order(P.pt), 
                                    ]
    Coh <- cbind(Coh = C.pt, ci.exp(m.APC, subset = "MCr", 
                                    ctr.mat = MCr[C.pos, , drop = FALSE], alpha = alpha))[order(C.pt), 
                                    ]
    colnames(Age)[-1] <- c("Rate", lu)
    colnames(Per)[-1] <- c("P-RR", lu)
    colnames(Coh)[-1] <- c("C-RR", lu)
    Type <- paste("ML of APC-model", Dist, ": (", parm, "):\n")
    Model <- m.APC
  } else {
    adc <- update(m.0, . ~ . - 1 + MA + I(P - A - c0))
    adp <- update(m.0, . ~ . - 1 + MA + I(P - p0))
    drift <- ci.exp(adc, subset = "I\\(")
    rownames(drift) <- "A-d"
    xP <- cbind(1, P - p0, MPr)
    xC <- cbind(1, P - A - c0, MCr)
    lP <- cbind(P - p0, MPr)
    lC <- cbind(P - A - c0, MCr)
    if (parm == "AD-C-P") {
      rc <- update(m.0, . ~ . - 1 + xC, offset = predict(adc, 
                                                         type = "link"))
      rp <- update(m.0, . ~ . - 1 + xP, offset = predict(adc, 
                                                         type = "link"))
      A.eff <- ci.exp(adc, subset = "MA", ctr.mat = MA[A.pos, 
      ], alpha = alpha)
      C.eff <- ci.exp(rc, subset = "xC", ctr.mat = xC[C.pos, 
      ], alpha = alpha)
      P.eff <- ci.exp(rp, subset = "xP", ctr.mat = xP[P.pos, 
      ], alpha = alpha)
      Model <- list(adc, rc, rp)
    } else if (parm == "AD-P-C") {
      rp <- update(m.0, . ~ . - 1 + xP, offset = predict(adp, 
                                                         type = "link"))
      rc <- update(m.0, . ~ . - 1 + xC, offset = predict(rp, 
                                                         type = "link"))
      A.eff <- ci.exp(adp, subset = "MA", ctr.mat = MA[A.pos, 
      ], alpha = alpha)
      P.eff <- ci.exp(rp, subset = "xP", ctr.mat = xP[P.pos, 
      ], alpha = alpha)
      C.eff <- ci.exp(rc, subset = "xC", ctr.mat = xC[C.pos, 
      ], alpha = alpha)
      Model <- list(adp, rp, rc)
    } else if (parm == "AC-P") {
      ac <- update(m.0, . ~ . - 1 + MA + lC)
      rp <- update(m.0, . ~ . - 1 + xP, offset = predict(ac, 
                                                         type = "link"))
      A.eff <- ci.exp(ac, subset = "MA", ctr.mat = MA[A.pos, 
      ], alpha = alpha)
      C.eff <- ci.exp(ac, subset = "lC", ctr.mat = lC[C.pos, 
      ], alpha = alpha)
      P.eff <- ci.exp(rp, subset = "xP", ctr.mat = xP[P.pos, 
      ], alpha = alpha)
      Model <- list(ac, rp)
    } else if (parm == "AP-C") {
      ap <- update(m.0, . ~ . - 1 + MA + lP)
      rc <- update(m.0, . ~ . - 1 + xC, offset = predict(ap, 
                                                         type = "link"))
      A.eff <- ci.exp(ap, subset = "MA", ctr.mat = MA[A.pos, 
      ], alpha = alpha)
      P.eff <- ci.exp(ap, subset = "lP", ctr.mat = lP[P.pos, 
      ], alpha = alpha)
      C.eff <- ci.exp(rc, subset = "xC", ctr.mat = xC[C.pos, 
      ], alpha = alpha)
      Model <- list(ap, rc)
    }
    Age <- cbind(Age = A.pt, A.eff)[order(A.pt), ]
    Per <- cbind(Per = P.pt, P.eff)[order(P.pt), ]
    Coh <- cbind(Cph = C.pt, C.eff)[order(C.pt), ]
    colnames(Age)[-1] <- c("A.eff", lu)
    colnames(Per)[-1] <- c("P.eff", lu)
    colnames(Coh)[-1] <- c("C.eff", lu)
    Type <- paste("Sequential modelling", Dist, ": (", parm, 
                  "):\n")
  }
  res <- list(Type = Type, Model = Model, Age = Age, Per = Per, 
              Coh = Coh, Drift = drift, Ref = c(Per = if (parm %in% 
                                                          c("APC", "ADPC", "Ad-P-C", "AP-C")) p0 else NA, Coh = if (parm %in% 
                                                                                                                    c("ACP", "ADCP", "Ad-C-P", "AC-P")) c0 else NA), 
              Anova = AOV)
  if (model %in% c("ns", "bs")) 
    res <- c(res, list(Knots = Knots))
  res$Age[, -1] <- res$Age[, -1] * scale
  if (print.AOV) {
    print(res$Type)
    print(res$Anova)
  }
  if (!ref.p & parm %in% c("APC", "ADPC")) 
    cat("No reference period given:\n", "Reference period for age-effects is chosen as\n", 
        "the median date of event: ", p0, ".\n")
  if (!ref.c & parm %in% c("ACP", "ADCP")) 
    cat("No reference period given:\n", "Reference period for age-effects is chosen as\n", 
        "the median date of birth for persons  with event: ", 
        c0, ".\n")
  class(res) <- "apc"
  invisible(res)
  
  # MCr backup
  MCr2 <- MCr
  
  data2 <- 
    data %>%
    mutate(d_spl = m.APC$fitted.values %>% as.double())
  
  # Fiting with cohort effects set to zero
  MCr <- MCr * 0
  
  m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
  
  data3 <- 
    data2 %>% 
    mutate(d_ap_spl = m.APC$fitted.values %>% as.double())
    
  ####
  #### (smoothed observed mortality) - (AP fitting from a full APC model)
  ####
  
  mtx <- data %>% 
    arrange(P, A)
  ylist <- unique(mtx$P) %>% sort()
  alist <- unique(mtx$A) %>% sort()
  deaths <- matrix(mtx$D, nrow=length(alist), ncol=length(ylist), byrow=F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  # population at risk (population/exposure)
  exposure <- matrix(mtx$Y, nrow=length(alist), ncol=length(ylist), byrow=F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  # smoothing mortality
  
  fit <- Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                      overdispersion = TRUE, method = 1)
  
  # fit <- Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
  #                     overdispersion = TRUE, method = 3, lambda = 0)
  
  summary(fit)
  mx.smooth <- (exp(fit[26]$logmortality))*100000
  d_smth <- (exp(fit[26]$logmortality))*exposure
  
  
  smt <- d_smth %>%
    as_tibble() %>% 
    mutate(A = alist) %>% 
    gather(-A, key = P, value = d_sth) %>% 
    mutate(A = as.integer(A),
           P = as.integer(P))
  
  data4 <- 
    data3 %>% 
    left_join(smt) %>% 
    rename(age = A,
           year = P,
           d = D,
           e = Y) %>% 
    mutate(sex = sx,
           code = cd)
  
  return(data4)
  
}

library(Epi)
library(MASS)
select <- dplyr::select

amin <- 0
amax <- 99
pmin <- 1940
pmax <- 2019

curv_f <- curvature_dAPC(amin=amin, amax=amax, pmin=pmin, pmax=pmax, gr=1, sx = "female", cd = "ESP")
curv_m <- curvature_dAPC(amin=amin, amax=amax, pmin=pmin, pmax=pmax, gr=1, sx = "male", cd = "ESP")

detach(package:MASS)
detach(package:Epi)

curvs <- 
  bind_rows(curv_m, curv_f) %>% 
  mutate(mx = 100000 * d / e,
         mx_sth = 100000 * d_sth / e,
         mx_spl = 100000 * d_spl / e,
         mx_ap_spl = 100000 * d_ap_spl / e,
         cohort = year - age,
         mx_exc = mx_sth - mx_ap_spl,
         mx_exc2 = mx - mx_ap_spl,
         rr = mx_sth / mx_ap_spl,
         rr2 = mx / mx_ap_spl,
         rr3 = ((mx / mx_ap_spl) - 1) * 100) 

sx <- "male"
nplot <- "4c"

{
qt <- 25

db_per <- 
  curvs %>% 
  filter(sex == sx)

# one for decrease of mortality (from green to blue), 
# another for mortality increases (from yellow to red)
col_scale <- 
  c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
    "springgreen",
    colorRampPalette(c("yellow", "red"), space = "Lab")(qt))

# Definition of brackets for the scale of change
val <- unique(db_per$rr2)
val <- log(unique(db_per$rr2))
# separate negative and positive values
pval <- val[val>0] # all positive values of change (mortality deterioration)
nval <- val[val<0] # all negative values of change (mortality improvement)

# identification of brackets for positive values (minimum + 23 quantiles + maximum)
# pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
pcop <- quantile(pval, prob = seq(0, 1, 1/qt))
pcop[length(pcop)] <- pcop[length(pcop)]*1.01
# the same as above but for negative values (minimum + 23 quantiles + maximum)
# ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
ncop <- quantile(nval, prob = seq(0, 1, 1/qt))
ncop[1] <- ncop[1]*1.01
# chain of brackets 25 ranges for negative changes, central value of no change (0), 
# and 25 ranges for positive changes
breaks_mc <- c(ncop, pcop) 
as.numeric(breaks_mc)
# transforming them to percentage values
(exp(breaks_mc)-1)*100

# adding to each value of change (continuous) the corresponding 
# bracket (a discrete interval)
db_per2 <- 
  db_per %>% 
  mutate(ch_cut = cut(rr3, breaks = (exp(breaks_mc)-1)*100)) 

cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
col_brks <- cuts[c(1, 10, 20, 25, 30, 40, 50)]

# assigning a color from our scale to each bracket
col_values <- setNames(col_scale, cuts)


library(cowplot)
lg <- 
  get_legend(db_per2 %>%
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  scale_fill_manual(name = "Excess\nMortality %",
                    breaks = col_brks,
                    values = col_values))

db_per2 %>%
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  scale_fill_manual(name = "Mortality deviance (%)",
                    breaks = col_brks,
                    values = col_values)+
  geom_vline(xintercept = seq(1950, 2020, 10), 
             linewidth = 0.2, linetype = "dashed", 
             alpha = 0.8, color = "grey30")+
  geom_hline(yintercept = seq(0, 100, 10), 
             linewidth = 0.2, linetype = "dashed", 
             alpha = 0.8, color = "grey30")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1950, 2020, 10)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0, 100, 5)) +
  
  geom_rect(xmin = 2019.65, xmax = 2032, ymin = -1, ymax = 112, 
            fill="white", col = "white")+
  geom_rect(xmin = 1949, xmax = 2032, ymin = 100.6, ymax = 110.5, 
            fill="white", col = "white")+
  
  geom_abline(intercept = seq(-1990, -1900, by = 10), slope = 1, 
              linetype="dashed", color="grey30", size=.15, alpha = 0.8)+
  
  annotate("text", x = 2022, y = seq(34, 94, by = 10),
           label = as.character(rev(seq(1930, 1990, by=10))),
           size = 2.2, angle = 45, color="grey20") +
  
  annotate("text", x = seq(1990, 2020, 10), y = 103,
           label = as.character(seq(1890, 1920, by=10)),
           size = 2.2, angle = 45, color="grey20") +
  annotate("text", x = 2023.5, y = 28.5, label = "Cohort", size = 2.6, 
           angle = 45, color="grey20") +
  
  annotate("text", x = 1945, y = 105, label = sx %>% str_to_title(), size = 3.2, 
           color="black", hjust = 0) +
  
  
  coord_equal(expand = 0,
              xlim = c(1939.5, 2028),
              ylim = c(0, 108))+
  labs(
    # subtitle = sx %>% str_to_title(),
       y = "Age", x = "Period")+
  theme_minimal()+
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.key.size = unit(0.5,"line"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 12),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    plot.background = element_rect(fill = "white", colour = "transparent")
  )
ggsave(paste0("figures/manuscript/fig0", nplot, "_cohort_patterns_over_time_", sx, "_", amin, "_", amax, "_", 
              pmin, "_", pmax, "_tst2.png"),
       w = 4, h = 6)


}



# db_per_m <- db_per2
# db_per_f <- db_per2











db_per2 %>%
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  scale_fill_manual(name = "Excess\nMortality %",
                    breaks = col_brks,
                    values = col_values)+
  geom_vline(xintercept = seq(1950, 2020, 10), 
             linewidth = 0.2, linetype = "dashed", 
             alpha = 0.8, color = "grey30")+
  geom_hline(yintercept = seq(0, 100, 10), 
             linewidth = 0.2, linetype = "dashed", 
             alpha = 0.8, color = "grey30")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1950, 2020, 10)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0, 100, 5)) +
  
  geom_rect(xmin = 2019.65, xmax = 2032, ymin = -1, ymax = 112, 
            fill="white", col = "white")+
  geom_rect(xmin = 1949, xmax = 2032, ymin = 100.6, ymax = 110.5, 
            fill="white", col = "white")+
  
  geom_abline(intercept = seq(-1990, -1900, by = 10), slope = 1, 
              linetype="dashed", color="grey30", size=.15, alpha = 0.8)+
  
  annotate("text", x = 2022, y = seq(34, 94, by = 10),
           label = as.character(rev(seq(1930, 1990, by=10))),
           size = 2, angle = 45, color="grey20") +
  
  annotate("text", x = seq(1990, 2020, 10), y = 103,
           label = as.character(seq(1890, 1920, by=10)),
           size = 2, angle = 45, color="grey20") +
  annotate("text", x = 2023.5, y = 28.5, label = "Cohort", size = 2.6, 
           angle = 45, color="grey20") +
  coord_equal(expand = 0,
              xlim = c(1939.5, 2028),
              ylim = c(0, 110))+
  labs(title = sx,
       y = "Age", x = "Period")+
  theme_minimal()+
  theme(
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 12),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    plot.background = element_rect(fill = "white", colour = "transparent")
  )


# sx <- "male"

curvs %>%
  filter(sex == sx,
         age %in% seq(0, 50, 10)) %>%
  select(year, age, starts_with("mx"), -mx_spl) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~age, scales = "free")+
  scale_y_log10()+
  geom_line(aes(year, mx, col = type))+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == "female",
         # year >= 1970,
         age %in% seq(0, 50, 10)) %>%
  select(year, age, cohort, starts_with("mx"), -mx_spl) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~age, scales = "free")+
  scale_y_log10()+
  geom_line(aes(cohort, mx, col = type))+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == "female",
         year %in% seq(1940, 1980, 10)) %>%
  select(year, age, cohort, starts_with("mx")) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~year, scales = "free")+
  geom_line(aes(age, mx, col = type))+
  scale_y_log10()+
  theme_bw()

# curvs %>%
#   mutate(cohort = year - age) %>%
#   filter(sex == sx,
#          year >= 1970,
#          age %in% seq(20, 80, 10)) %>%
#   select(year, age, cohort, starts_with("mx")) %>%
#   gather(starts_with("mx"), key = type, value = mx) %>%
#   filter(type %in% c("mx_exc2", "mx_exc")) %>%
#   ggplot()+
#   facet_wrap(~age, scales = "free")+
#   geom_line(aes(cohort, mx, col = type))+
#   geom_hline(yintercept = 0, lty = "dashed")+
#   theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == sx,
         year >= 1970,
         age %in% seq(20, 80, 10)) %>%
  # select(year, age, cohort, starts_with("mx")) %>%
  # gather(starts_with("mx"), key = type, value = mx) %>%
  # mutate(log_mx = log(mx))+
  # filter(type %in% c("mx_exc2", "mx_exc")) %>%
  ggplot()+
  facet_wrap(~age, scales = "free_x")+
  geom_line(aes(cohort, rr))+
  geom_hline(yintercept = 1, lty = "dashed")+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == sx,
         year >= 1960) %>%
  ggplot()+
  geom_tile(aes(year, age, fill = rr2))+
  scale_fill_gradient2(midpoint = 1, low = "blue", mid = "white", high = "red")+
  scale_y_continuous(breaks=seq(20, 90, 10)) +
  # geom_vline(xintercept = seq(1960, 2010, by=10), linetype="dashed", color="grey30", size=.15) +
  # geom_hline(yintercept = seq(20, 90, by=10), linetype="dashed", color="grey30", size=.15) +
  # geom_abline(intercept = seq(-2010, -1860, by=10),
  #             slope = 1, linetype="dashed", color="grey30", size=.15)+
  coord_fixed(expand = 0)+
  theme_bw()


