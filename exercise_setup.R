# Setup only: synthetic data and function definitions.
# Self-contained teaching example. Not an empirical ECFIN replication.
# Requires base R plus the recommended splines package.
simulate_panel <- function(N=100, TT=240, burn=100,
                           state=c("iid","endogenous"), seed=20261016) {
  state <- match.arg(state)
  stopifnot(N>=2, TT>=10, burn>=0)
  set.seed(seed)
  x <- rnorm(TT+burn); common <- rnorm(TT+burn,sd=.3)
  y <- rep(0,N); z <- runif(N,-1,1); pieces <- vector("list",TT)
  for(t in seq_len(TT+burn)) {
    oldy <- y; oldz <- z
    y <- .6*oldy+(1-oldz^2)*x[t]+common[t]+rnorm(N,sd=.5)
    if(t>burn) pieces[[t-burn]] <- data.frame(id=seq_len(N),t=t-burn,
       y=y,x=x[t],z=oldz,ylag=oldy)
    # z saved above is PRE-shock. This update defines the next period's state.
    z <- if(state=="iid") runif(N,-1,1) else
      tanh(.75*oldz+.25*x[t]+.3*rnorm(N))
  }
  panel <- do.call(rbind,pieces)
  panel <- panel[order(panel$id,panel$t),]
  rownames(panel) <- NULL
  panel
}
d <- simulate_panel()
out <- "output_R_paper"
dir.create(out,showWarnings=FALSE)

# 1. Fit a horizon LP and form the covariance of the full coefficient vector.
# This full-regressor sandwich gives the same sieve block as partialling out W
# in equations (13)-(16), with all sample-size normalizations cancelled.
ols_date_hac <- function(D,y,time,L) {
  stopifnot(length(L)==1, is.finite(L), L>=0, L==floor(L),
            all(is.finite(D)), all(is.finite(y)), qr(D)$rank==ncol(D))
  b <- drop(qr.solve(D,y)); u <- drop(y-D%*%b)
  S <- rowsum(D*u,time,reorder=TRUE)
  TT <- nrow(S); stopifnot(TT>L,all(diff(sort(unique(time)))==1))
  meat <- crossprod(S)
  if(L>0) for(k in seq_len(L)) {
    C <- crossprod(S[(k+1):TT,,drop=FALSE],S[1:(TT-k),,drop=FALSE])
    meat <- meat+(1-k/(L+1))*(C+t(C))
  }
  bread <- solve(crossprod(D)); V <- bread%*%meat%*%bread
  list(b=b,V=(V+t(V))/2,u=u)
}

# 2. Construct J cubic B-spline columns, including a constant in their span.
# J=4 is a global cubic; J>4 adds J-4 empirical-quantile interior knots.
make_basis <- function(z,J) {
  stopifnot(length(J)==1, J>=4, J==floor(J), all(is.finite(z)), diff(range(z))>0)
  knots <- if(J>4) as.numeric(quantile(z,seq_len(J-4)/(J-3),type=7)) else numeric(0)
  if(anyDuplicated(knots)) stop("Repeated knots: reduce J or reconsider the state.")
  splines::bs(z,knots=knots,degree=3,intercept=TRUE,Boundary.knots=range(z))
}

# 3. Horizon-specific estimation and AIC choice, for the teaching example.
estimate_horizon <- function(h,Jset=4:20,bandwidth=NULL,data=d) {
  d <- data
  stopifnot(length(h)==1, h>=0, h==floor(h), length(Jset)>0,
            all(Jset>=4), all(Jset==floor(Jset)),
            !anyDuplicated(d[c("id","t")]))
  key <- paste(d$id,d$t,sep=":")
  yh <- d$y[match(paste(d$id,d$t+h,sep=":"),key)]
  ok <- complete.cases(yh,d$x,d$z,d$ylag)
  z <- d$z[ok]; x <- d$x[ok]; y <- yh[ok]; time <- d$t[ok]
  W <- cbind(ylag=d$ylag[ok],constant=1)
  TT <- length(unique(time))
  # Teaching bandwidth rule. Examine sensitivity to longer lags.
  L <- if(is.null(bandwidth)) floor(4*(TT/100)^(2/9)) else bandwidth
  models <- lapply(Jset,function(J) {
    B <- make_basis(z,J); D <- cbind(B*x,W)
    b <- drop(qr.solve(D,y)); u <- drop(y-D%*%b)
    list(J=J,B=B,D=D,AIC=length(y)*log(mean(u^2))+2*ncol(D))
  })
  aic <- vapply(models,function(m)m$AIC,numeric(1))
  best <- models[[which.min(aic)]]
  fit <- ols_date_hac(best$D,y,time,L)
  J <- best$J
  list(h=h,J=J,B=best$B,b=fit$b[seq_len(J)],V=fit$V[seq_len(J),seq_len(J),drop=FALSE],
       AIC=data.frame(J=Jset,AIC=aic),n=length(y),dates=TT,L=L,
       support=range(z),central=quantile(z,c(.05,.95)))
}

# 4. Pointwise intervals and a studentized Gaussian band across a stated grid.
# A grid approximates a continuum. These bands are for ONE horizon, not all h.
# Sieve approximation bias must be small for coverage of the true function.
# AIC selection alone does not establish that condition in an empirical dataset.
infer_curve <- function(fit,grid,Bdraws=1999,alpha=.05,seed=20261014) {
  stopifnot(length(grid)>1, all(is.finite(grid)), Bdraws>=99,
            alpha>0, alpha<1, min(grid)>=fit$support[1],max(grid)<=fit$support[2])
  P <- predict(fit$B,newx=grid)
  g <- drop(P%*%fit$b)
  se <- sqrt(pmax(rowSums((P%*%fit$V)*P),0))
  stopifnot(all(se>0))
  ev <- eigen(fit$V,symmetric=TRUE)
  if(min(ev$values)< -1e-8*max(abs(ev$values))) stop("Covariance not positive semidefinite.")
  Vroot <- sweep(ev$vectors,2,sqrt(pmax(ev$values,0)),"*")
  set.seed(seed+fit$h)
  G <- P%*%Vroot%*%matrix(rnorm(fit$J*Bdraws),nrow=fit$J)
  sup_t <- apply(abs(G)/se,2,max)
  crit <- unname(quantile(sup_t,1-alpha,type=7))
  zcrit <- qnorm(1-alpha/2)
  data.frame(h=fit$h,z=grid,estimate=g,se,
    point_lower=g-zcrit*se,point_upper=g+zcrit*se,
    band_lower=g-crit*se,band_upper=g+crit*se,
    band_critical=crit,J=fit$J,dates=fit$dates,bandwidth=fit$L)
}

# 5. A pre-specified contrast: response at z=0 minus response at z=.75.
# Both states must be inside the data support; use the SAME fitted function.
infer_contrast <- function(fit,zA=0,zB=.75) {
  stopifnot(min(zA,zB)>=fit$support[1],max(zA,zB)<=fit$support[2])
  P <- predict(fit$B,newx=c(zA,zB)); a <- P[1,]-P[2,]
  est <- sum(a*fit$b); se <- sqrt(max(drop(a%*%fit$V%*%a),0))
  data.frame(h=fit$h,zA,zB,estimate=est,se,lower=est-qnorm(.975)*se,upper=est+qnorm(.975)*se)
}



# Exact population calculations, with a Monte Carlo illustration, for a
# finite SHIFT of the realized shock. Distinct from an intervention from 0.
nonlinear_shift_demo <- function(z=0,delta=1,cubic=.25,
                                 draws=200000,seed=20261108) {
  set.seed(seed)
  x <- rnorm(draws); g <- 1-z^2
  f <- function(v) g*v+cubic*v^3
  data.frame(z=z,delta=delta,cubic=cubic,
    LP_times_shift=(g+3*cubic)*delta,
    finite_shift=g*delta+cubic*(3*delta+delta^3),
    discrepancy=cubic*delta^3,
    MC_shift=mean(f(x+delta)-f(x)))
}
