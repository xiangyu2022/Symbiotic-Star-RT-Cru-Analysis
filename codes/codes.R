# Load the data and necessary library
library(LPBkg)
library(LPsmooth)
library(Hmisc)
library(VGAM)
options(digits=15)
bs1<-read.table("/Users/zxysmacbook/Downloads/16688_src.txt")
bs2<-read.table("/Users/zxysmacbook/Downloads/18710_src.txt")
bg1<-read.table("/Users/zxysmacbook/Downloads/16688_bkg.txt")
bg2<-read.table("/Users/zxysmacbook/Downloads/18710_bkg.txt")
bs_full<-c(bs1[,1],bs2[,1])
bs_abs<-abs(bs_full)
bg_full<-c(bg1[,1],bg2[,1])
bg_abs<-abs(bg_full)
length(bg_abs)
length(bs_abs)
B=10000;M=10

# MERGE LAST THREE REGIONS. TRY 16.8 - 17.2 
# n_1 n_2 SHOULD BE ADDED TO THE TABLE, 
# want to make the notations different I_r, Band of interest
# W/B background notation not f_br rather than bg bkg, 
# \newcomand 


# Set the Wavelength range for regions Ir, with r=1,...,9.
L1=1.65;U1=2.05;L2=5.5;U2=10.2;L3=11.5;U3=13;
L4=14.6;U4=15.15;L5=16.2;U5=17.4;L6=18.5;U6=19.5
L7=21.3;U7=21.75;L8=21.65;U8=22;L9=21.9;U9=22.4

# Set the Wavelength range for combined regions Cw, with w=1,...,5.
L11=1.65; U11=2.05
L21=5.5; U21=10.2
L31=11.5; U31=13
L41=14.6; U41=17.4
L51=18.5; U51=22.4

# Set the data for each region Ir, Cw; bs corresponds to the physics samples while bg corresponds 
# to the background only data. 
bs1<-bs_abs[bs_abs>=L1&bs_abs<=U1]
bs2<-bs_abs[bs_abs>=L2&bs_abs<=U2]
bs3<-bs_abs[bs_abs>=L3&bs_abs<=U3]
bs4<-bs_abs[bs_abs>=L4&bs_abs<=U4]
bs5<-bs_abs[bs_abs>=L5&bs_abs<=U5]
bs6<-bs_abs[bs_abs>=L6&bs_abs<=U6]
bs7<-bs_abs[bs_abs>=L7&bs_abs<=U7]
bs8<-bs_abs[bs_abs>=L8&bs_abs<=U8]
bs9<-bs_abs[bs_abs>=L9&bs_abs<=U9]

bg1<-bg_abs[bg_abs>=L1&bg_abs<=U1]
bg2<-bg_abs[bg_abs>=L2&bg_abs<=U2]
bg3<-bg_abs[bg_abs>=L3&bg_abs<=U3]
bg4<-bg_abs[bg_abs>=L4&bg_abs<=U4]
bg5<-bg_abs[bg_abs>=L5&bg_abs<=U5]
bg6<-bg_abs[bg_abs>=L6&bg_abs<=U6]
bg7<-bg_abs[bg_abs>=L7&bg_abs<=U7]
bg8<-bg_abs[bg_abs>=L8&bg_abs<=U8]
bg9<-bg_abs[bg_abs>=L9&bg_abs<=U9]

bg11<-bg_abs[bg_abs>=L11&bg_abs<=U11]
bg21<-bg_abs[bg_abs>=L21&bg_abs<=U21]
bg31<-bg_abs[bg_abs>=L31&bg_abs<=U31]
bg41<-bg_abs[bg_abs>=L41&bg_abs<=U41]
bg51<-bg_abs[bg_abs>=L51&bg_abs<=U51]


# Number of samples for two datasets of each region Ir. This is used to generate the 
# n1 and n2 of Table 2.
(n1=length(bg1));(n2=length(bg2));(n3=length(bg3));(n4=length(bg4));(n5=length(bg5));
(n6=length(bg6));(n7=length(bg7));(n8=length(bg8));(n9=length(bg9))
(n1=length(bs1));(n2=length(bs2));(n3=length(bs3));(n4=length(bs4));(n5=length(bs5));
(n6=length(bs6));(n7=length(bs7));(n8=length(bs8));(n9=length(bs9))


# Number of samples for two datasets of each region Cw. This is used to generate the 
# n1 Table 3.
(n1=length(bg11));(n2=length(bg21));(n3=length(bg31));(n4=length(bg41));(n5=length(bg51));


# Postulated uniform background for region Cw, equivalant with the function fbr(x) in paper. 
unibkg11<-function(x)ifelse(x>=L11&x<=U11,1/(U11-L11),0)
unibkg21<-function(x)ifelse(x>=L21&x<=U21,1/(U21-L21),0)
unibkg31<-function(x)ifelse(x>=L31&x<=U31,1/(U31-L31),0)
unibkg41<-function(x)ifelse(x>=L41&x<=U41,1/(U41-L41),0)
unibkg51<-function(x)ifelse(x>=L51&x<=U51,1/(U51-L51),0)

# A function which achieves the process of methods in section 3 and returns the three
# adjusted p-values and the number of non-zero coefficients.  
dhatL2_new <- function (data, g, M = 9, Mmax = NULL, smooth = TRUE, criterion = "BIC", 
                        hist.u = TRUE, breaks = 20, ylim = c(0, 2.5), range = c(min(data), 
                                                                                max(data)), sigma = 2) 
{
  bluetrans <- rgb(0, 0, 250, 50, maxColorValue = 300)
  pinktrans <- rgb(30, 0, 10, 30, maxColorValue = 300)
  G <- function(y) integrate(g, lower = range[1], upper = y)$value
  G <- Vectorize(G)
  u <- G(data)
  xx <- seq(range[1], range[2], by = 0.01)
  uu <- G(xx)
  n <- length(u)
  S <- as.matrix(Legj(u = u, m = M))
  LP <- apply(S, FUN = "mean", 2)
  if (smooth == TRUE) 
    LP <- denoise(LP, n = n, criterion)
  IDS0 <- which(LP == 0)
  IDS1 <- which(LP != 0)
  z.L2 <- approxExtrap(u, 1 + S %*% LP,ties = "mean", xout = c(0, u, 1))
  dhat <- approxfun(z.L2$x, z.L2$y, rule = 2,ties = "mean")
  covs <- cov(S) * (n - 1)/n^2
  covs[IDS0, ] <- 0
  covs[, IDS0] <- 0
  vec1 <- rep(0, length(u))
  for (j in 1:M) {
    for (k in 1:M) {
      vec1 <- vec1 + S[, j] * S[, k] * covs[j, k]
    }
  }
  z.SE.L2 <- approxExtrap(u, sqrt(vec1), xout = c(0, u, 1),ties = "mean")
  SEdhat <- approxfun(z.SE.L2$x, z.SE.L2$y, rule = 2,ties = "mean")
  sigmas0 <- rep(1/n, M)
  sigmas0[IDS0] <- 0
  z.SE0.L2 <- approxExtrap(u, sqrt(apply(S^2 %*% sigmas0, 1, 
                                         sum)), xout = c(0, u, 1),ties = "mean")
  SEdhatH0 <- approxfun(z.SE0.L2$x, z.SE0.L2$y, rule = 2,ties = "mean")
  alpha = pnorm(sigma, 0, 1, lower.tail = FALSE)
  if (Mmax > 1 & smooth == FALSE) 
    alpha = alpha/Mmax
  if (Mmax > 1 & smooth == TRUE) 
    alpha = 2 * alpha/(Mmax + M * (M - 1))
  qa <- c_alpha2(M, IDS1, alpha = alpha, c_interval = c(1, 
                                                        20))
  LBf1 <- function(u) 1 - qa * SEdhatH0(u)
  UBf1 <- function(u) 1 + qa * SEdhatH0(u)
  SE <- function(u) SEdhat(u)
  deviance <- n * sum(LP^2)
  Kstat <- n*max((LP^2))
  K_adj_pval <- 1-pchisq(Kstat,df=1,ncp=0)^Mmax
  Naive_adj_pval <- pchisq(deviance,Mmax, lower.tail = FALSE)
  pval <- pchisq(deviance, sum(LP != 0), lower.tail = FALSE)
  b_adj_pval = min(pval * Mmax,1)
  
  if (hist.u == TRUE) {
    oldpar <- par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
    on.exit(par(oldpar))
    par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
    hist(u, breaks = breaks, prob = TRUE, main = " ", ylim = ylim, 
         ylab = "Comparison density", xlab = expression(G(x)), 
         cex.axis = 1, cex.lab = 1.4, xaxt = "n")
    lines(uu, dhat(uu), col = "dodgerblue1", lwd = 2)
    polygon(c(uu, rev(uu)), c(LBf1(uu), rev(UBf1(uu))), col = pinktrans, 
            border = FALSE)
    polygon(c(uu, rev(uu)), c(dhat(uu) - SE(uu), rev(dhat(uu) + 
                                                       SE(uu))), col = bluetrans, border = FALSE)
    abline(h = 1, lwd = 2, lty = 2, col = "tomato3")
    Gx <- quantile(data, probs = seq(0, 1, by = 0.1))
    Axis(side = 1, at = seq(0, 1, by = 0.1), labels = paste("G(", 
                                                            round(Gx, 4), ")", sep = ""), cex.axis = 1, col = "black", 
         col.lab = "black")
  }
  dhat_num <- function(x) dhat(G(x))
  dhat.x <- function(x) dhat_num(x)
  f <- function(x) g(x) * dhat(G(x))
  return(list(Bon_adj_pvalue = b_adj_pval, Naive_adj_pval=Naive_adj_pval,
              K_adj_pval=K_adj_pval,
              kstar = sum(LP != 0),LP = LP))
}

# This gives the background correction results as listed in Table 3.
bkg_results1=dhatL2_new(bg11, unibkg11, Mmax = 10)
bkg_results2=dhatL2_new(bg21, unibkg21, Mmax = 10)
bkg_results3=dhatL2_new(bg31, unibkg31, Mmax = 10)
bkg_results4=dhatL2_new(bg41, unibkg41, Mmax = 10)
bkg_results5=dhatL2_new(bg51, unibkg51, Mmax = 10)
bkg_results1;bkg_results2;bkg_results3;bkg_results4;bkg_results5
sidak_like = function(x,M=9){return(1-(1-x)^(M))}

sidak_like(3.76874938784572e-10,5)
sidak_like(3.67970409078824e-06,5)
sidak_like(3.76875197716231e-10,5)



# After getting the results above, we can correct the background as:
bkg1<-function(x)ifelse(x>=L1&x<=U1,1/(U1-L1),0)
bkg2<-function(x)ifelse(x>=L21&x<=U21,0.242353 - 0.00376911*x,0)
bkg3<-function(x)ifelse(x>=L3&x<=U3,1/(U3-L3),0)
bkg4<-function(x)ifelse(x>=L4&x<=U4,2.4749 - 0.0441494*x,0)
bkg5<-function(x)ifelse(x>=L5&x<=U5,1.18995 - 0.0212274*x,0) #16.2-17.4
bkg6<-function(x)ifelse(x>=L6&x<=U6,1/(U6-L6),0)
bkg7<-function(x)ifelse(x>=L7&x<=U7,1/(U7-L7),0)
bkg8<-function(x)ifelse(x>=L8&x<=U8,1/(U8-L8),0)
bkg9<-function(x)ifelse(x>=L9&x<=U9,1/(U9-L9),0)


# This gives the signal detection results as listed in Table 4.
dhat_result1_new=dhatL2_new(bs1,bkg1,Mmax=10)
dhat_result2_new=dhatL2_new(bs2,bkg2,Mmax=10)
dhat_result3_new=dhatL2_new(bs3,bkg3,Mmax=10) # New region
dhat_result4_new=dhatL2_new(bs4,bkg4,Mmax=10)
dhat_result5_new=dhatL2_new(bs5,bkg5,Mmax=10)
dhat_result6_new=dhatL2_new(bs6,bkg6,Mmax=10)
dhat_result7_new=dhatL2_new(bs7,bkg7,Mmax=10)
dhat_result8_new=dhatL2_new(bs8,bkg8,Mmax=10)
dhat_result9_new=dhatL2_new(bs9,bkg9,Mmax=10)

sidak_like(dhat_result1_new$Bon_adj_pvalue,9)
sidak_like(dhat_result1_new$Naive_adj_pval,9)
sidak_like(dhat_result1_new$K_adj_pval,9)

sidak_like(dhat_result2_new$Bon_adj_pvalue,9)
sidak_like(dhat_result2_new$Naive_adj_pval,9)
sidak_like(dhat_result2_new$K_adj_pval,9)

# This is how we generate the comparison density plot, as listed in Figure 1.
# The function here involves simulations, and will be very time-consuming. 
bkg2_n=function(x) bkg2(x)/integrate(bkg2,5.5,10.2)$value
sampler_xb2=function(n){
  xq=runif(100*n,L2,U2)
  V=runif(100*n,0,1)
  m=max(bkg2_n(seq(L2,U2,length=10000)))*(U2-L2)
  x=xq[V<= bkg2_n(xq)*(U2-L2)/m]
  return (x[1:n])
}

g=bkg2_n
L=L2;U=U2;data=bs2;range=c(L2,U2); lattice=NULL;B=10000;R=1000;criterion="BIC";selection=T;M=10
samplerG=sampler_xb2
n <- length(data)
uu <- seq(0, 1, length = R)
D <- c(); K<-c()
Mat_ddG <- Mat_ddF <- matrix(rep(0, R * B), nrow = B, ncol = R)
dd <- d_hat(data, m=M, g, range, lattice, selection, criterion)
ddx <- function(x) dd$dx(x)
ddx <- Vectorize(ddx)
ddu <- function(u) dd$du(u)
ddu <- Vectorize(ddu)
D_obs <- sum(dd$LPj^2)
K_obs <- max(dd$LPj^2)
G <- function(x) integrate(g, L, x)$value
G <- Vectorize(G)
fhat <- function(x) g(x) * ddx(x)
fhat <- Vectorize(fhat)
MM <- max(dd$du(uu))
for (b in 1:B) {
  message(paste("Performing smoothed bootstrap...iteration n.", 
                b))
  xG <- samplerG(round(n * 2 * MM))
  dataG <- sample(xG, n)
  ddG <- d_hat(dataG, m=10, g, range, lattice, selection, 
               criterion)
  D[b] <- sum(ddG$LPj^2)
  K[b] <- max(ddG$LPj^2)
  dduG <- function(u) ddG$du(u)
  dduG <- Vectorize(dduG)
  Mat_ddG[b, ] <- dduG(uu)
  v <- runif(round(n * 2 * MM))
  dataF <- sample(xG[v * MM <= dd$du(G(xG))], n)
  ddF <- d_hat(dataF, m=10, g, range, lattice, selection, 
               criterion)
  dduF <- function(u) ddF$du(u)
  dduF <- Vectorize(dduF)
  Mat_ddF[b, ] <- dduF(uu)
}

SD <- SDf <- SDftrue <- SDfemp <- c()
processG <- matrix(rep(0, R * B), nrow = B, ncol = R)
for (r in 1:R) {
  SD[r] <- sd(Mat_ddG[, r])
  SDf[r] <- sd(Mat_ddF[, r])
  processG[, r] <- (Mat_ddG[, r] - 1)/sd(Mat_ddG[, r])
}
MaxDistr <- apply(abs(processG), 1, max)
calpha <- quantile(MaxDistr, 0.99)
p.value <- mean(n * D >= n * D_obs)
p.value_K <- mean(n * K >= n * K_obs)
CD.plot = TRUE;ylim=c(0,2)
if (CD.plot == TRUE) {
  Transparency <- function(color, shade = 80) {
    color2 <- col2rgb(color)
    apply(color2, 2, function(color3) {
      rgb(red = color3[1], green = color3[2], blue = color3[3], 
          alpha = shade, maxColorValue = 250)
    })
  }
  transgreen <- Transparency("chartreuse2")
  transgrey <- Transparency("grey67")
  oldpar <- par(mfrow = c(1, 1), mar = c(5, 6, 1, 1))
  on.exit(par(oldpar))
  par(mfrow = c(1, 1), mar = c(5, 6, 1, 1))
  plot(uu, ddu(uu), type = "l", main = " ", ylim = ylim, 
       ylab = expression(hat(d)(u, G, F)), xlab = "u", cex.axis = 2, 
       cex.lab = 2, col = "darkgreen", lwd = 2, xaxt = "n")
  abline(h = 1, lwd = 2, lty = 2, col = "tomato3")
  if (is.null(lattice)) {
    polygon(c(uu, rev(uu)), c(rep(1, length(uu)) - calpha * 
                                SD, rev(rep(1, length(uu)) + calpha * SD)), col = transgrey, 
            border = FALSE)
    polygon(c(uu, rev(uu)), c(ddu(uu) - SDf, rev(ddu(uu) + 
                                                   SDf)), col = transgreen, border = FALSE)
    Gx <- quantile(data, probs = seq(0, 1, by = 0.1))
    Axis(side = 1, at = seq(0, 1, by = 0.1), labels = paste("G(", 
                                                            round(Gx, 4), ")", sep = ""), cex.axis = 1, col = "black", 
         col.lab = "black")
  }
  else {
    polygon.step <- function(x, y1, y2, border = FALSE, 
                             ...) {
      nx <- length(x)
      ny <- length(y1)
      if (length(y2) != ny) 
        stop("y1 and y2 must be the same length")
      if (nx != (ny + 1)) 
        stop("x must be one longer than y")
      xx <- c(x[1], rep(x[-c(1, nx)], rep(2, nx - 2)), 
              x[nx])
      xxx <- c(xx, rev(xx))
      yy1 <- rep(y1, rep(2, ny))
      yy2 <- rep(y2, rep(2, ny))
      yyy <- c(yy1, rev(yy2))
      polygon(xxx, yyy, border = border, ...)
    }
    polygon.step(c(0, uu), 1 - calpha * SD, 1 + calpha * 
                   SD, col = transgrey)
    polygon.step(c(0, uu), ddu(uu) - SDf, ddu(uu) + SDf, 
                 col = transgreen)
  }
  legend("topleft", legend = c(expression(hat(d)(u, G, 
                                                 F)), expression("99% confidence bands"), "Standard error"), 
         lwd = 1, col = c("darkgreen", transgrey, transgreen), 
         pch = c(NA, 22, 22), bty = "n", lty = c(1, NA, NA), 
         pt.bg = c(transgreen, transgrey), cex = 1)
  Axis(side = 1, at = seq(0, 1, by = 0.1), cex.axis = 2, 
       cex = 2, col = "black", col.lab = "black")
}
list(Deviance = D_obs, Kstatistic=K_obs,p_value_dev = p.value,p_value_K=p.value_K)
results2=list(Deviance = D_obs, Kstatistic=K_obs,p_value_dev = p.value,p_value_K=p.value_K)

# Examples of codes used for the simulation to generate upper limits and power curves. 
sampler_xb3=function(n){
  xq=runif(20*n,L3,U3)
  V=runif(20*n,0,1)
  m=max(bkg3(seq(L3,U3,length=10000)))*(U3-L3)
  x=xq[V<= bkg3(xq)*(U3-L3)/m]
  return (x[1:n])
}
sampler_xb4=function(n){
  xq=runif(20*n,L4,U4)
  V=runif(20*n,0,1)
  m=max(bkg4(seq(L4,U4,length=10000)))*(U4-L4)
  x=xq[V<= bkg4(xq)*(U4-L4)/m]
  return (x[1:n])
}
sampler_xb5=function(n){
  xq=runif(20*n,L5,U5)
  V=runif(20*n,0,1)
  m=max(bkg5(seq(L5,U5,length=10000)))*(U5-L5)
  x=xq[V<= bkg5(xq)*(U5-L5)/m]
  return (x[1:n])
}
sampler_xb6=function(n){
  xq=runif(20*n,L6,U6)
  V=runif(20*n,0,1)
  m=max(bkg6(seq(L6,U6,length=10000)))*(U6-L6)
  x=xq[V<= bkg6(xq)*(U6-L6)/m]
  return (x[1:n])
}
sampler_xb7=function(n){
  xq=runif(20*n,L7,U7)
  V=runif(20*n,0,1)
  m=max(bkg7(seq(L7,U7,length=10000)))*(U7-L7)
  x=xq[V<= bkg7(xq)*(U7-L7)/m]
  return (x[1:n])
}
sampler_xb8=function(n){
  xq=runif(20*n,L8,U8)
  V=runif(20*n,0,1)
  m=max(bkg8(seq(L8,U8,length=10000)))*(U8-L8)
  x=xq[V<= bkg8(xq)*(U8-L8)/m]
  return (x[1:n])
}
sampler_xb9=function(n){
  xq=runif(20*n,L9,U9)
  V=runif(20*n,0,1)
  m=max(bkg9(seq(L9,U9,length=10000)))*(U9-L9)
  x=xq[V<= bkg9(xq)*(U9-L9)/m]
  return (x[1:n])
}

f_signal3=function(x) 1/((1+((x-rnorm(1,12.131,sd=0.005))/0.05)^2)^2.5)
nor3=integrate(f_signal3,L3,U3)$value
f_signal_nor3=function(x)f_signal3(x)/nor3

f_signal4=function(x) 1/((1+((x-rnorm(1,15.014,sd=0.005))/0.05)^2)^2.5)
nor4=integrate(f_signal4,L4,U4)$value
f_signal_nor4=function(x)f_signal4(x)/nor4

f_signal5=function(x) 1/((1+((x-rnorm(1,17.051,sd=0.005))/0.05)^2)^2.5)
nor5=integrate(f_signal5,L5,U5)$value
f_signal_nor5=function(x)f_signal5(x)/nor5

# f_signal5=function(x) 1/((1+((x-rnorm(1,16.951,sd=0.005))/0.05)^2)^2.5)
# nor5=integrate(f_signal5,L5,U5)$value
# f_signal_nor5=function(x)f_signal5(x)/nor5

f_signal6=function(x) 1/((1+((x-rnorm(1,18.967,sd=0.005))/0.05)^2)^2.5)
nor6=integrate(f_signal6,L6,U6)$value
f_signal_nor6=function(x)f_signal6(x)/nor6

f_signal7=function(x) 1/((1+((x-rnorm(1,21.602,sd=0.005))/0.05)^2)^2.5)
nor7=integrate(f_signal7,L7,U7)$value
f_signal_nor7=function(x)f_signal7(x)/nor7

f_signal8=function(x) 1/((1+((x-rnorm(1,21.804,sd=0.005))/0.05)^2)^2.5)
nor8=integrate(f_signal8,L8,U8)$value
f_signal_nor8=function(x)f_signal8(x)/nor8

f_signal9=function(x) 1/((1+((x-rnorm(1,22.101,sd=0.005))/0.05)^2)^2.5)
nor9=integrate(f_signal9,L9,U9)$value
f_signal_nor9=function(x)f_signal9(x)/nor9


sampler_xs3=function(n){
  xq=runif(100*n,L3,U3)
  V=runif(100*n,0,1)
  m=max(f_signal_nor3(seq(L3,U3,length=1000)))*(U3-L3)
  x=xq[V<= f_signal_nor3(xq)*(U3-L3)/m]
  return (x[1:n])
}
sampler_xs4=function(n){
  xq=runif(50*n,L4,U4)
  V=runif(50*n,0,1)
  m=max(f_signal_nor4(seq(L4,U4,length=1000)))*(U4-L4)
  x=xq[V<= f_signal_nor4(xq)*(U4-L4)/m]
  return (x[1:n])
}
sampler_xs5=function(n){
  xq=runif(40*n,L5,U5)
  V=runif(40*n,0,1)
  m=max(f_signal_nor5(seq(L5,U5,length=2000)))*(U5-L5)
  x=xq[V<= f_signal_nor5(xq)*(U5-L5)/m]
  return (x[1:n])
}
sampler_xs6=function(n){
  xq=runif(40*n,L6,U6)
  V=runif(40*n,0,1)
  m=max(f_signal_nor6(seq(L6,U6,length=2000)))*(U6-L6)
  x=xq[V<= f_signal_nor6(xq)*(U6-L6)/m]
  return (x[1:n])
}
sampler_xs7=function(n){
  xq=runif(40*n,L7,U7)
  V=runif(40*n,0,1)
  m=max(f_signal_nor7(seq(L7,U7,length=1000)))*(U7-L7)
  x=xq[V<= f_signal_nor7(xq)*(U7-L7)/m]
  return (x[1:n])
}
sampler_xs8=function(n){
  xq=runif(40*n,L8,U8)
  V=runif(40*n,0,1)
  m=max(f_signal_nor8(seq(L8,U8,length=1000)))*(U8-L8)
  x=xq[V<= f_signal_nor8(xq)*(U8-L8)/m]
  return (x[1:n])
}
sampler_xs9=function(n){
  xq=runif(40*n,L9,U9)
  V=runif(40*n,0,1)
  m=max(f_signal_nor9(seq(L9,U9,length=1000)))*(U9-L9)
  x=xq[V<= f_signal_nor9(xq)*(U9-L9)/m]
  return (x[1:n])
}

G3=function(x) ifelse(x>=L3&x<=U3,integrate(bkg3,L3,upper=x)$value,0)
G3=Vectorize(G3)
G4=function(x) ifelse(x>=L4&x<=U4,integrate(bkg4,L4,upper=x)$value,0)
G4=Vectorize(G4)
G5=function(x) ifelse(x>=L5&x<=U5,integrate(bkg5,L5,upper=x)$value,0)
G5=Vectorize(G5)
G6=function(x) ifelse(x>=L6&x<=U6,integrate(bkg6,L6,upper=x)$value,0)
G6=Vectorize(G6)
G7=function(x) ifelse(x>=L7&x<=U7,integrate(bkg7,L7,upper=x)$value,0)
G7=Vectorize(G7)
G8=function(x) ifelse(x>=L8&x<=U8,integrate(bkg8,L8,upper=x)$value,0)
G8=Vectorize(G8)
G9=function(x) ifelse(x>=L9&x<=U9,integrate(bkg9,L9,upper=x)$value,0)
G9=Vectorize(G9)

#Likelihood Ratio Test
ns=c(n1,n2,n3,n4,n5,n6,n7,n8,n9)

lsf<-function(w0,w){1/(1+((w-w0)/0.05)^2)^2.5}
k_lsf<-function(w)integrate(lsf,lower=L,upper=U,w=w)$value
k_lsf<-Vectorize(k_lsf)
fe<-function(w0,w){1/(1+((w-w0)/0.05)^2)^2.5}
line_int<-function(w0,w,c,sigma){dnorm(w0,c,sigma)*fe(w0,w)}
fw<-function(w,c,sigma){integrate(line_int,lower=L,upper=U,w=w,c=c,sigma=sigma)$value/((pnorm(U,c,sigma)-pnorm(L,c,sigma))*k_lsf(w))}
fw<-Vectorize(fw)

bs_model_new<-function(x,pars){
  (1-sum(pars[1:3]))*bkg1(x)+pars[1]*fw(x,1.78499,pars[4])+pars[2]*fw(x,1.85247,pars[5])+ 
    pars[3]*fw(x,1.94365,pars[6])
}

fe<-function(x,w0){1/(1+((x-w0)/0.05)^2)^2.5}
int1 <- function(w0) integrate(fe,lower=L1,upper=U1,w0=w0)$value
fe1_norm <- function(x,w0) fe(x,w0)/int1(w0)
ll_bs_new1<-function(pars,x){
  -sum(log((1-sum(pars[1:3]))*bkg1(x)+pars[1]*fe1_norm(x,1.78499)+
             pars[2]*fe1_norm(x,1.85247)+pars[3]*fe1_norm(x,1.94365)))
}

int2 <- function(w0) integrate(fe,lower=L2,upper=U2,w0=w0)$value
fe2_norm <- function(x,w0) fe(x,w0)/int2(w0)
ll_bs_new2<-function(pars,x){
  -sum(log((1-sum(pars[1:3]))*bkg2(x)+pars[1]*fe2_norm(x,6.180)+
             pars[2]*fe2_norm(x,7.171)+pars[3]*fe2_norm(x,8.419)))
}

int3 <- function(w0) integrate(fe,lower=L3,upper=U3,w0=w0)$value
fe3_norm <- function(x,w0) fe(x,w0)/int3(w0)
ll_bs_new3<-function(pars,x){
  -sum(log((1-pars[1])*bkg3(x)+pars[1]*fe3_norm(x,rnorm(1,12.131,0.005))))
}

int4 <- function(w0) integrate(fe,lower=L4,upper=U4,w0=w0)$value
fe4_norm <- function(x,w0) fe(x,w0)/int4(w0)
ll_bs_new4<-function(pars,x){
  -sum(log((1-pars[1])*bkg4(x)+pars[1]*fe4_norm(x,rnorm(1,15.014,0.005))))
}

int5 <- function(w0) integrate(fe,lower=L5,upper=U5,w0=w0)$value
fe5_norm <- function(x,w0) fe(x,w0)/int5(w0)
ll_bs_new5<-function(pars,x){
  -sum(log((1-pars[1])*bkg5(x)+pars[1]*fe5_norm(x,rnorm(1,17.051,0.005))))
}


int6 <- function(w0) integrate(fe,lower=L6,upper=U6,w0=w0)$value
fe6_norm <- function(x,w0) fe(x,w0)/int6(w0)
ll_bs_new6<-function(pars,x){
  -sum(log((1-pars[1])*bkg6(x)+pars[1]*fe6_norm(x,rnorm(1,18.967,0.005))))
}

int7 <- function(w0) integrate(fe,lower=L7,upper=U7,w0=w0)$value
fe7_norm <- function(x,w0) fe(x,w0)/int7(w0)
ll_bs_new7<-function(pars,x){
  -sum(log((1-pars[1])*bkg7(x)+pars[1]*fe7_norm(x,rnorm(1,21.602,0.005))))
}

int8 <- function(w0) integrate(fe,lower=L8,upper=U8,w0=w0)$value
fe8_norm <- function(x,w0) fe(x,w0)/int8(w0)
ll_bs_new8<-function(pars,x){
  -sum(log((1-pars[1])*bkg8(x)+pars[1]*fe8_norm(x,rnorm(1,21.804,0.005))))
}

int9 <- function(w0) integrate(fe,lower=L9,upper=U9,w0=w0)$value
fe9_norm <- function(x,w0) fe(x,w0)/int9(w0)
ll_bs_new9<-function(pars,x){
  -sum(log((1-pars[1])*bkg9(x)+pars[1]*fe9_norm(x,rnorm(1,22.101,0.005))))
}



MLE3 = optimize(f=ll_bs_new3, interval=c(0,1), x=bs3)
MLE4 = optimize(f=ll_bs_new4, interval=c(0,1), x=bs4)
MLE5 = optimize(f=ll_bs_new5, interval=c(0,1), x=bs5)
MLE6 = optimize(f=ll_bs_new6, interval=c(0,1), x=bs6)
MLE7 = optimize(f=ll_bs_new7, interval=c(0,1), x=bs7)
MLE8 = optimize(f=ll_bs_new8, interval=c(0,1), x=bs8)
MLE9 = optimize(f=ll_bs_new9, interval=c(0,1), x=bs9)

(pval_like3=pchisq(-2*(sum(log(bkg3(bs3))) + MLE3$objective),df=1,lower.tail = FALSE)/2)
(pval_like4=pchisq(-2*(sum(log(bkg4(bs4)))+ MLE4$objective),df=1,lower.tail = FALSE)/2) #
(pval_like5=pchisq(-2*(sum(log(bkg5(bs5)))+ MLE5$objective),df=1,lower.tail = FALSE)/2) #
(pval_like6=pchisq(-2*(sum(log(bkg6(bs6))) + MLE6$objective),df=1,lower.tail = FALSE)/2)
(pval_like7=pchisq(-2*(sum(log(bkg7(bs7))) + MLE7$objective),df=1,lower.tail = FALSE)/2) #
(pval_like8=pchisq(-2*(sum(log(bkg8(bs8))) + MLE8$objective),df=1,lower.tail = FALSE)/2)
(pval_like9=pchisq(-2*(sum(log(bkg9(bs9)))+ MLE9$objective),df=1,lower.tail = FALSE)/2)





# The simulation will be very time-consuming, which we provide the simulation results below directly. 
# Here, instead of running the simulation yourself, we recommend you to load the results obtained 
# by us directly MNRAS224808.RData (provided in our github page)
etas=c(0,0.005, 0.025, 0.05,0.075,0.1,0.125,0.150,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4)
for (i in etas){
  for (k in 3:9){
    init1="D_eta."
    init2="Msel_eta."
    init3="T1_eta."
    init11=paste(init1,i,"_",k,"=c()",sep="")
    init21=paste(init2,i,"_",k,"=c()",sep="")
    init31=paste(init3,i,"_",k,"=c()",sep="")
    init41=paste("LRT_eta.",i,"_",k,"=c()",sep="")
    eval(parse(text=init11))
    eval(parse(text=init21))
    eval(parse(text=init31))
    eval(parse(text=init41))
  }
}

(n3=length(bs3));(n4=length(bs4));(n5=length(bs5));
(n6=length(bs6));(n7=length(bs7));(n8=length(bs8));(n9=length(bs9))
ns=c(n3,n4,n5,n6,n7,n8,n9)
B=10000
M=10

load("MNRAS224808.RData")

# for (i in 3:9){
#   for(eta in etas){
#     for (b in 1:B){
#       init1=paste("xb=sampler_xb",i,"(",ns[i-2],")",sep="")
#       eval(parse(text=init1))
#       init2=paste("xs=sampler_xs",i,"(",ns[i-2],")",sep="")
#       eval(parse(text=init2))
#       init3=paste("Mat<-matrix(rep(0,",ns[i-2],"*","2),ncol=2,nrow=",ns[i-2],")",sep="")
#       eval(parse(text=init3))
#       Mat[,1]=xb
#       Mat[,2]=xs
#       init4=paste("flags=t(rmultinom(",ns[i-2], ",size=rep(1,2), prob=c(1-",eta,",",eta,")))",sep="")
#       eval(parse(text=init4))
#       # flags <- t(rmultinom(10*n, size=rep(1,2), prob=c(1-eta,eta)))
#       xx <-as.numeric(apply(Mat*flags,1,sum))
#       init5=paste("uu=G",i,"(xx)",sep="")
#       eval(parse(text=init5))
#       S <- as.matrix(Legj(u = uu, m = M))
#       LP <- apply(S, FUN = "mean", 2)
#       init9=paste("LPsel<-denoise(LP,",ns[i-2],",method=\"BIC\") ",sep="")
#       eval(parse(text=init9))
#       init6=paste("Msel_eta.",eta,"_",i,"[b]=sum(LPsel!=0)",sep="")
#       init7=paste("D_eta.",eta,"_",i,"[b]=",ns[i-2],"*","sum(LPsel^2)",sep="")
#       init8=paste("T1_eta.",eta,"_",i,"[b]=",ns[i-2],"*","max(LP^2)",sep="")
#       init10=paste("LRT_eta.",eta,"_",i,"[b]=-2*(sum(log(bkg",i,"(xx)))+optimize(f=ll_bs_new",i,",interval=c(0,1), x=xx)$objective)",sep="")
#       eval(parse(text=init6))
#       eval(parse(text=init7))
#       eval(parse(text=init8))
#       eval(parse(text=init10))
#       print(c(i,eta,b))
#     }
#   }
# }

# The simulation results are summarized in table frame_n3,...,frame_n9 as: (take region 3 as an example):
frame_n3=data.frame(Msel_eta.0_3=Msel_eta.0_3,D_eta.0_3=D_eta.0_3,T1_eta.0_3=T1_eta.0_3,
                    Msel_eta.0.005_3=Msel_eta.0.005_3,D_eta.0.005_3=D_eta.0.005_3,T1_eta.0.005_3=T1_eta.0.005_3,
                    Msel_eta.0.025_3=Msel_eta.0.025_3,D_eta.0.025_3=D_eta.0.025_3,T1_eta.0.025_3=T1_eta.0.025_3,
                    Msel_eta.0.05_3=Msel_eta.0.05_3,D_eta.0.05_3=D_eta.0.05_3,T1_eta.0.05_3=T1_eta.0.05_3,
                    Msel_eta.0.075_3=Msel_eta.0.075_3,D_eta.0.075_3=D_eta.0.075_3,T1_eta.0.075_3=T1_eta.0.075_3,
                    Msel_eta.0.1_3=Msel_eta.0.1_3,D_eta.0.1_3=D_eta.0.1_3,T1_eta.0.1_3=T1_eta.0.1_3,
                    Msel_eta.0.125_3=Msel_eta.0.125_3,D_eta.0.125_3=D_eta.0.125_3,T1_eta.0.125_3=T1_eta.0.125_3,
                    Msel_eta.0.15_3=Msel_eta.0.15_3,D_eta.0.15_3=D_eta.0.15_3,T1_eta.0.15_3=T1_eta.0.15_3,
                    Msel_eta.0.175_3=Msel_eta.0.175_3,D_eta.0.175_3=D_eta.0.175_3,T1_eta.0.175_3=T1_eta.0.175_3,
                    Msel_eta.0.2_3=Msel_eta.0.2_3,D_eta.0.2_3=D_eta.0.2_3,T1_eta.0.2_3=T1_eta.0.2_3,
                    Msel_eta.0.25_3=Msel_eta.0.25_3,D_eta.0.25_3=D_eta.0.25_3,T1_eta.0.25_3=T1_eta.0.25_3,
                    Msel_eta.0.3_3=Msel_eta.0.3_3,D_eta.0.3_3=D_eta.0.3_3,T1_eta.0.3_3=T1_eta.0.3_3,
                    Msel_eta.0.35_3=Msel_eta.0.35_3,D_eta.0.35_3=D_eta.0.35_3,T1_eta.0.35_3=T1_eta.0.35_3,
                    Msel_eta.0.4_3=Msel_eta.0.4_3,D_eta.0.4_3=D_eta.0.4_3,T1_eta.0.4_3=T1_eta.0.4_3
)


pvalbonf<-function(xdf){
  pval<-ifelse(sum(xdf)!=0,M*pchisq(xdf[1],xdf[2],lower.tail=F),1)
  return(pval)
}


pval_like2=c();pval_like=c()
for (i in 3:9){
  initp = paste("pval_like[i-2]=mean((pchisq(LRT_eta.0_",i,",df=1,lower.tail = FALSE)/2)<alpha)",sep="")
  eval(parse(text=initp))
  initp2 = paste("pval_like2[i-2]=mean((pchisq(LRT_eta.0_",i,",df=1,lower.tail = FALSE)/2)<alphas)",sep="")
  eval(parse(text=initp2))
  initpo=paste("Power_like_",i,"=Power_like_sidak_",i,"=c()",sep="")
  eval(parse(text=initpo))
  for (j in 1:length(etas[1:18][-1])){
    initpo1 = paste("Power_like_",i,"[",j,"] = mean((pchisq(LRT_eta.",etas[j],"_",i,",df=1,lower.tail = FALSE)/2)<alpha)",sep="")
    eval(parse(text=initpo1))
    initpo2 = paste("Power_like_sidak_",i,"[",j,"] = mean((pchisq(LRT_eta.",etas[j],"_",i,",df=1,lower.tail = FALSE)/2)<alphas)",sep="")
    eval(parse(text=initpo2))
  }
}

alpha=0.05
sidak=function(x,M=9){return(1-(1-x)^(1/M))}
alphas=sidak(alpha,7)
pval_like2=c()
for (i in 3:9){
  call0=paste("frame_n",i,"=data.frame(Msel_eta.0_",i,"=Msel_eta.0_",i,",D_eta.0_",i,"=D_eta.0_",i,",T1_eta.0_",i,"=T1_eta.0_",i,",Msel_eta.0.005_",i,"=Msel_eta.0.005_",i,",D_eta.0.005_",i,"=D_eta.0.005_",i,",T1_eta.0.005_",i,"=T1_eta.0.005_",i,",Msel_eta.0.025_",i,"=Msel_eta.0.025_",i,",D_eta.0.025_",i,"=D_eta.0.025_",i,",T1_eta.0.025_",i,"=T1_eta.0.025_",i,",Msel_eta.0.05_",i,"=Msel_eta.0.05_",i,",D_eta.0.05_",i,"=D_eta.0.05_",i,",T1_eta.0.05_",i,"=T1_eta.0.05_",i,",Msel_eta.0.075_",i,"=Msel_eta.0.075_",i,",D_eta.0.075_",i,"=D_eta.0.075_",i,",T1_eta.0.075_",i,"=T1_eta.0.075_",i,",Msel_eta.0.1_",i,"=Msel_eta.0.1_",i,",D_eta.0.1_",i,"=D_eta.0.1_",i,",T1_eta.0.1_",i,"=T1_eta.0.1_",i,",Msel_eta.0.125_",i,"=Msel_eta.0.125_",i,",D_eta.0.125_",i,"=D_eta.0.125_",i,",T1_eta.0.125_",i,"=T1_eta.0.125_",i,",Msel_eta.0.15_",i,"=Msel_eta.0.15_",i,",D_eta.0.15_",i,"=D_eta.0.15_",i,",T1_eta.0.15_",i,"=T1_eta.0.15_",i,",Msel_eta.0.175_",i,"=Msel_eta.0.175_",i,",D_eta.0.175_",i,"=D_eta.0.175_",i,",T1_eta.0.175_",i,"=T1_eta.0.175_",i,",Msel_eta.0.2_",i,"=Msel_eta.0.2_",i,",D_eta.0.2_",i,"=D_eta.0.2_",i,",T1_eta.0.2_",i,"=T1_eta.0.2_",i,",Msel_eta.0.225_",i,"=Msel_eta.0.225_",i,",D_eta.0.225_",i,"=D_eta.0.225_",i,",T1_eta.0.225_",i,"=T1_eta.0.225_",i,",Msel_eta.0.25_",i,"=Msel_eta.0.25_",i,",D_eta.0.25_",i,"=D_eta.0.25_",i,",T1_eta.0.25_",i,"=T1_eta.0.25_",i,",Msel_eta.0.275_",i,"=Msel_eta.0.275_",i,",D_eta.0.275_",i,"=D_eta.0.275_",i,",T1_eta.0.275_",i,"=T1_eta.0.275_",i,",Msel_eta.0.3_",i,"=Msel_eta.0.3_",i,",D_eta.0.3_",i,"=D_eta.0.3_",i,",T1_eta.0.3_",i,"=T1_eta.0.3_",i,",Msel_eta.0.325_",i,"=Msel_eta.0.325_",i,",D_eta.0.325_",i,"=D_eta.0.325_",i,",T1_eta.0.325_",i,"=T1_eta.0.325_",i,",Msel_eta.0.35_",i,"=Msel_eta.0.35_",i,",D_eta.0.35_",i,"=D_eta.0.35_",i,",T1_eta.0.35_",i,"=T1_eta.0.35_",i,",Msel_eta.0.375_",i,"=Msel_eta.0.375_",i,",D_eta.0.375_",i,"=D_eta.0.375_",i,",T1_eta.0.375_",i,"=T1_eta.0.375_",i,",Msel_eta.0.4_",i,"=Msel_eta.0.4_",i,",D_eta.0.4_",i,"=D_eta.0.4_",i,",T1_eta.0.4_",i,"=T1_eta.0.4_",i,",Msel_eta.0.45_",i,"=Msel_eta.0.45_",i,",D_eta.0.45_",i,"=D_eta.0.45_",i,",T1_eta.0.45_",i,"=T1_eta.0.45_",i,",Msel_eta.0.5_",i,"=Msel_eta.0.5_",i,",D_eta.0.5_",i,"=D_eta.0.5_",i,",T1_eta.0.5_",i,"=T1_eta.0.5_",i,")",sep="")
  call1 = paste("Power_T1_",i,"<-MCerror_power_T1_",i,"<-Power_Bonf_",i,"<-MCerror_power_Bonf_",i,"<-Power_naive_",i,"<-MCerror_power_naive_",i,"<-c()",sep="")
  call2 = paste("Power_T1_sidak",i,"<-MCerror_power_T1_sidak",i,"<-Power_Bonf_sidak",i,"<-MCerror_power_Bonf_sidak",i,"<-Power_naive_sidak",i,"<-MCerror_power_naive_sidak",i,"<-c()",sep="")
  eval(parse(text=call0))
  eval(parse(text=call1))
  eval(parse(text=call2))
  
  initp = paste("pval_like=mean((pchisq(LRT_eta.0_",i,",df=1,lower.tail = FALSE)/2)<alpha)",sep="")
  eval(parse(text=initp))
  initp2 = paste("pval_like2=mean((pchisq(LRT_eta.0_",i,",df=1,lower.tail = FALSE)/2)<alphas)",sep="")
  eval(parse(text=initp2))
  initpo=paste("Power_like_",i,"=Power_like_sidak_",i,"=c()",sep="")
  eval(parse(text=initpo))
  for (j in 1:length(etas[1:14][-1])){
    initpo1 = paste("Power_like_",i,"[",j,"] = mean((pchisq(LRT_eta.",etas[j],"_",i,",df=1,lower.tail = FALSE)/2)<alpha)",sep="")
    eval(parse(text=initpo1))
    initpo2 = paste("Power_like_sidak_",i,"[",j,"] = mean((pchisq(LRT_eta.",etas[j],"_",i,",df=1,lower.tail = FALSE)/2)<alphas)",sep="")
    eval(parse(text=initpo2))
  }
}


for (j in 3:9){
  for(i in 1:17){
    call11 =  paste("T1<-frame_n",j,"[,3*i+3]",sep="")
    eval(parse(text=call11))
    #T1<-frame_n3[,3*i+3]
    pval_T1<-1-(pchisq(T1,1))^M
    #Power_T1_10n[i]<-mean(pval_T1<alpha)
    call2 = paste("Power_T1_",j,"[i]<-mean(pval_T1<alpha)",sep="")
    call2s = paste("Power_T1_sidak",j,"[i]<-mean(pval_T1<alphas)",sep="")
    eval(parse(text=call2))
    eval(parse(text=call2s))
    #MCerror_power_T1_10n[i]<-sqrt(Power_T1_10n[i]*(1- Power_T1_10n[i])/B)
    call3 = paste("MCerror_power_T1_",j,"[i]<-sqrt(Power_T1_",j,"[i]*(1- Power_T1_",j,"[i])/B)",sep ="")
    call3s = paste("MCerror_power_T1_sidak",j,"[i]<-sqrt(Power_T1_sidak",j,"[i]*(1- Power_T1_sidak",j,"[i])/B)",sep ="")
    eval(parse(text=call3))
    eval(parse(text=call3s))
    
    #Msel<-frame_n3[,3*i+1] # Deviance
    call4 = paste("Msel = frame_n",j,"[,3*i+1]",sep="")
    eval(parse(text=call4))
    
    #Dev<-frame_n3[,3*i+2] # Mselected or Degree of Freedom.
    call5 = paste("Dev = frame_n",j,"[,3*i+2]",sep="")
    eval(parse(text=call5))
    
    pval_Bonf<-apply(cbind(Dev,Msel),1,pvalbonf)
    #Power_Bonf_10n[i]<-mean(pval_Bonf<alpha)
    call6 =  paste("Power_Bonf_",j,"[i]<-mean(pval_Bonf<alpha)",sep="")
    call6s =  paste("Power_Bonf_sidak",j,"[i]<-mean(pval_Bonf<alphas)",sep="")
    eval(parse(text=call6))
    eval(parse(text=call6s))
    
    #MCerror_power_Bonf_10n[i]<-sqrt( Power_Bonf_10n[i]*(1- Power_Bonf_10n[i])/B)
    call7  =  paste("MCerror_power_Bonf_",j,"[i]<-sqrt( Power_Bonf_",j,"[i]*(1- Power_Bonf_",j,"[i])/B)",sep="")
    call7s =  paste("MCerror_power_Bonf_sidak",j,"[i]<-sqrt( Power_Bonf_sidak",j,"[i]*(1 - Power_Bonf_sidak",j,"[i])/B)",sep="")
    eval(parse(text=call7))
    eval(parse(text=call7s))
    
    pval_naive<-pchisq(Dev,M,lower.tail = F)
    call8 = paste("Power_naive_",j,"[i]<-mean(pval_naive<alpha)",sep="")
    call8s = paste("Power_naive_sidak",j,"[i]<-mean(pval_naive<alphas)",sep="")
    call9 = paste("MCerror_power_naive_",j,"[i]<-sqrt(Power_naive_",j,"[i]*(1-Power_naive_",j,"[i])/B)",sep="")
    call9s = paste("MCerror_power_naive_sidak",j,"[i]<-sqrt(Power_naive_sidak",j,"[i]*(1-Power_naive_sidak",j,"[i])/B)",sep="")
    
    eval(parse(text=call8))
    eval(parse(text=call8s))
    eval(parse(text=call9))
    eval(parse(text=call9s))
    #Power_naive_10n[i]<-mean(pval_naive<alpha)
    #MCerror_power_naive_10n[i]<-sqrt( Power_naive_10n[i]*(1- Power_naive_10n[i])/B)
    print(i)
  }
}


# Power plot for region 3
par(mar=c(4.5,5,3,0.5))
plot(etas_new[-1][1:9],Power_naive_3[1:9],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main=expression("n"[3]*"=730"),cex.main=2.5,lwd=2.5,font.main=3)
points(etas_new[-1][1:9],Power_Bonf_3[1:9],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_3[1:9],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_3[1:9],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

points(etas_new[-1][1:9],Power_naive_sidak3[1:9],pch=2,col="dodgerblue",cex=2,type="o",lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_Bonf_sidak3[1:9],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_sidak3[1:9],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_sidak_3[1:9],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood
Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2.5,cex=2.5, col = "black", col.lab = "black")
Axis(side = 1, at = c(0.005,0.025,0.05, 0.075, 0.1, 0.125,0.15, 0.175, 0.2), cex.axis=2.5,cex=2.5,col = "black", col.lab = "black")
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
# round(c(0.062,0.063,0.084,0.073,0.0875,0.092)*730,2)
# round(c(0.096,0.098,0.117,0.106,0.122,0.122)*730,2)
abline(v=0.072,col="grey")
# c(0.041,0.067,0.054,0.073)*730

# Power plot for region 4
plot(etas_new[-1][1:13],Power_naive_4[1:13],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main=expression("n"[4]*"=247"),cex.main=2.5,lwd=2.5)
points(etas_new[-1][1:13],Power_Bonf_4[1:13],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:13],Power_T1_4[1:13],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:13],Power_like_4[1:13],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

points(etas_new[-1][1:13],Power_naive_sidak4[1:13],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:13],Power_Bonf_sidak4[1:13],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:13],Power_T1_sidak4[1:13],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:13],Power_like_sidak_4[1:13],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood
Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2,cex=2, col = "black", col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:13], cex.axis=2,cex=2, col = "black", col.lab = "black")
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
# c(0.117,0.127,0.156,0.138,0.172,0.173)*247
# c(0.181,0.195,0.217,0.20,0.242,0.233)*247

#c(0.081,0.131,0.107,0.16)*247

# Power plot for region 5
plot(etas_new[-1][1:9],Power_naive_5[1:9],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main="n5=471",cex.main=2.5,lwd=2.5)
points(etas_new[-1][1:9],Power_Bonf_5[1:9],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_5[1:9],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_5[1:9],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

points(etas_new[-1][1:9],Power_naive_sidak5[1:9],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_Bonf_sidak5[1:9],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_sidak5[1:9],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_sidak_5[1:9],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood

Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2,cex=2, col = "black",  col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:9], cex.axis=2,cex=2, col = "black",  col.lab = "black")
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
round(c(0.075,0.078,0.101,0.088,0.107,0.111)*471,2)
round(c(0.117,0.121,0.139,0.125,0.147,0.149)*471,2)
c(0.051,0.075,0.064,0.093)*471

# Power plot for region 6
plot(etas_new[-1][1:9],Power_naive_6[1:9],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main="n6=390",cex.main=2.5,lwd=2.5
)
points(etas_new[-1][1:9],Power_Bonf_6[1:9],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_6[1:9],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_6[1:9],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

#legend("topleft",legend=c(expression(Bonferroni), expression(T1), expression(Naive)),pch=c(15,16,17),
#       col=c("tomato3","darkgreen","dodgerblue"),cex=1)
points(etas_new[-1][1:9],Power_naive_sidak6[1:9],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_Bonf_sidak6[1:9],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:9],Power_T1_sidak6[1:9],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:9],Power_like_sidak_6[1:9],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood

Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2.5,cex=2., col = "black", col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:9], cex.axis=2,cex=2, col = "black", col.lab = "black")

# This is how we calculate upper limits for region 6. 
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
c(0.084,0.088,0.113,0.1,0.121,0.125)*390
c(0.131,0.138,0.156,0.144,0.171,0.168)*390
c(0.058,0.089,0.072,0.101)*390

# Power plot for region 7 in Figure 4. 
plot(etas_new[-1][1:17],Power_naive_7[1:17],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main="n7=179",cex.main=2.5,lwd=2.5
)
points(etas_new[-1][1:17],Power_Bonf_7[1:17],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_7[1:17],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
#legend("topleft",legend=c(expression(Bonferroni), expression(T1), expression(Naive)),pch=c(15,16,17),
#       col=c("tomato3","darkgreen","dodgerblue"),cex=1)
points(etas_new[-1][1:17],Power_like_7[1:17],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

points(etas_new[-1][1:17],Power_naive_sidak7[1:17],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_Bonf_sidak7[1:17],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_sidak7[1:17],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:17],Power_like_sidak_7[1:17],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood
Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2,cex=2, col = "black",  col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:17], cex.axis=2,cex=2, col = "black",  col.lab = "black")

# This is how we calculate upper limits for region 7 .
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
round(c(0.15,0.168,0.205,0.183,0.231,0.23)*179,2)
round(c(0.235,0.256,0.28,0.26,0.317,0.3)*179,2)
abline(v=0.201)
round(c(0.1,0.166,0.135,0.201)*179,2)

# Power plot for region 8 in Figure 4. 
plot(etas_new[-1][1:17],Power_naive_8[1:17],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main="n8=145",cex.main=2.5,lwd=2.5)
points(etas_new[-1][1:17],Power_Bonf_8[1:17],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_8[1:17],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:17],Power_like_8[1:17],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

#legend("topleft",legend=c(expression(Bonferroni), expression(T1), expression(Naive)),pch=c(15,16,17),
#       col=c("tomato3","darkgreen","dodgerblue"),cex=1)

points(etas_new[-1][1:17],Power_naive_sidak8[1:17],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Naive
points(etas_new[-1][1:17],Power_Bonf_sidak8[1:17],pch=0,col="tomato3",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_sidak8[1:17],pch=1,col="darkgreen",cex=2,type="o", lty=2,lwd=2.5) # T1
points(etas_new[-1][1:17],Power_like_sidak_8[1:17],pch=5,col="grey1",cex=2,type="o", lty=2,lwd=2.5) # likelihood
Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2,cex=2, col = "black",  col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:17], cex.axis=2,cex=2, col = "black",  col.lab = "black")

# This is how we calculate upper limits for region 8.
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
abline(v=0.25,col="grey")
round(c(0.185,0.206,0.259,0.229,0.275,0.291)*145,2)
round(c(0.288,0.314,0.35,0.321,0.375,0.384)*145,2)
round(c(0.123,0.209,0.171,0.25)*145,2)

# Power plot for region 9 in Figure 4. 
plot(etas_new[-1][1:17],Power_naive_9[1:17],pch=17,col="dodgerblue",xlab=expression(eta),ylim=c(0,1),
     ylab="Power",cex=2.5,type="o",cex.axis=2.5,cex.lab=2.5,xaxt='n',yaxt='n',
     main="n9=390",cex.main=2.5,lwd=2.5)
points(etas_new[-1][1:17],Power_Bonf_9[1:17],pch=15,col="tomato3",cex=2,type="o",lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_9[1:17],pch=16,col="darkgreen",cex=2,type="o",lwd=2.5) # T1
points(etas_new[-1][1:17],Power_like_9[1:17],pch=23,col="grey1",bg="grey1",cex=2,type="o",lwd=2.5) # likelihood

points(etas_new[-1][1:17],Power_naive_sidak9[1:17],pch=2,col="dodgerblue",cex=2,type="o", lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_Bonf_sidak9[1:17],pch=0,col="tomato3",cex=2,type="o",lty=2,lwd=2.5) # Bon
points(etas_new[-1][1:17],Power_T1_sidak9[1:17],pch=1,col="darkgreen",cex=2,type="o",lty=2,lwd=2.5) # T1
points(etas_new[-1][1:17],Power_like_sidak_9[1:17],pch=5,col="grey1",cex=2,type="o",lty=2,lwd=2.5) # likelihood
Axis(side = 2, at = c(0.1,0.3,0.5,0.7,0.9), cex.axis=2,cex=2, col = "black",  col.lab = "black")
Axis(side = 1, at = etas_new[-1][1:17], cex.axis=2,cex=2, col = "black",  col.lab = "black")

# This is how we calculate upper limits for region 9.
abline(a=0.5,b=0,col="grey");abline(a=0.9,b=0,col="grey")
abline(v=0.197,col="grey")
c(0.15,0.17,0.204,0.183,0.234,0.229)*390
c(0.235,0.258,0.282,0.26,0.321,0.303)*390
c(0.097,0.163,0.133,0.197)*390


# Legend plot in Figure 4. 
par(mar=c(4.5,5,3,0.5))
plot(1,ylim=c(0,0.9),xaxt='n',yaxt='n',cex.main=2.5,ylab="",xlab="")
legend("center",legend=c(expression(Bonferroni), expression(K), expression(Naive),expression(LRT),
                         expression("Bonferroni with Sidak's correction"), expression("K with Sidak's correction"), 
                         expression("Naive with Sidak's correction"),expression("LRT with Sidak's correction")),pch=c(15,16,17,23,0,1,2,5),pt.bg=c("tomato3","darkgreen","dodgerblue","grey1"),
       col=c("tomato3","darkgreen","dodgerblue","grey1"),cex=2.4,lty=c(1,1,1,1,2,2,2,2), lwd=2.5,
       bty = "n")

TypeIerror_T1 <- TypeIerror_Bonferroni<-TypeIerror_naive<-MCerror_T1<-MCerror_Bonferroni<-MCerror_naive<-c()

# The p-values can be calculated as follows: 
for (i in 3:9){
  #K
  call1=paste("One_minus_pval_T1_eta0<-(pchisq(frame_n",i,"[,3],1))^M",sep="")
  eval(parse(text=call1))
  TypeIerror_T1[i-2]<-mean(One_minus_pval_T1_eta0>1-alpha)
  MCerror_T1[i-2] <- sqrt(TypeIerror_T1[i-2]*(1-TypeIerror_T1[i-2])/B)
  
  #Bonferroni
  call2=paste("pval_Bonf_eta0<-apply(cbind(frame_n",i,"[,2],frame_n",i,"[,1]),1,pvalbonf)",sep="")
  eval(parse(text=call2))
  TypeIerror_Bonferroni[i-2]<-mean(pval_Bonf_eta0<alpha)
  MCerror_Bonferroni[i-2] <- sqrt(TypeIerror_Bonferroni[i-2]*(1-TypeIerror_Bonferroni[i-2])/B)
  
  #Naive correction
  call3=paste("pval_naive_eta0<-pchisq(frame_n",i,"[,3],M,lower.tail = F)",sep="")
  eval(parse(text=call3))
  TypeIerror_naive[i-2]<-mean(pval_naive_eta0<alpha)
  MCerror_naive[i-2] <- sqrt(TypeIerror_naive[i-2]*(1-TypeIerror_naive[i-2])/B)
}

TypeIerror_T1;TypeIerror_Bonferroni;TypeIerror_naive
MCerror_T1;MCerror_Bonferroni;MCerror_naive


# 14.908 & 16.93 in section 5.3. 
# Constrction of upper limits and power curves can be got follows exactly the same way as above. 
ll_bs_new4_new<-function(pars,x){
  -sum(log((1-pars[1])*bkg4(x)+pars[1]*fe4_norm(x,rnorm(1,14.908,0.005))))
}

ll_bs_new5_new<-function(pars,x){
  -sum(log((1-pars[1])*bkg5(x)+pars[1]*fe5_norm(x,rnorm(1,16.93,0.005))))
}
set.seed(13)
MLE4_new = optimize(f=ll_bs_new4_new, interval=c(0,1), x=bs4)
MLE5_new = optimize(f=ll_bs_new5_new, interval=c(0,1), x=bs5)
(pval_like4n_new=pchisq(-2*(sum(log(bkg4(bs4)))+ MLE4_new$objective),df=1,lower.tail = FALSE)/2) #
(pval_like5n_new=pchisq(-2*(sum(log(bkg5(bs5)))+ MLE5_new$objective),df=1,lower.tail = FALSE)/2) #

