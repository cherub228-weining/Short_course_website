# Run from the repository root. Shared estimation functions live in one file.
source("code/R/exercise_setup.R")

# LAB A: fit and inspect the impact response before looping over horizons.
impact <- estimate_horizon(0)
print(impact$AIC); cat("Selected J:",impact$J,"; HAC lags:",impact$L,"\n")
print(infer_contrast(impact))

# LAB B: plot pointwise intervals and the wider state-grid simultaneous band.
# A fixed grid [-.8,.8] is used for comparability across horizons in this data.
fits <- lapply(0:8,estimate_horizon)
grid <- seq(-.8,.8,length.out=161)
curves <- do.call(rbind,lapply(fits,infer_curve,grid=grid))
contrasts <- do.call(rbind,lapply(fits,infer_contrast))
aic <- do.call(rbind,lapply(fits,function(f)cbind(h=f$h,f$AIC)))
write.csv(curves,file.path(out,"paper_curves.csv"),row.names=FALSE)
write.csv(contrasts,file.path(out,"paper_contrasts.csv"),row.names=FALSE)
write.csv(aic,file.path(out,"AIC_by_horizon.csv"),row.names=FALSE)
pdf(file.path(out,"pointwise_and_uniform_h0.pdf"),width=8,height=5,family="Helvetica")
s <- subset(curves,h==0)
plot(s$z,s$estimate,type="n",ylim=range(s$band_lower,s$band_upper),
 xlab="Pre-shock state z",ylab="Response per unit shock",main="Synthetic panel: horizon 0")
polygon(c(s$z,rev(s$z)),c(s$band_lower,rev(s$band_upper)),col="#007D8222",border=NA)
lines(s$z,s$point_lower,col="#007D82",lty=3)
lines(s$z,s$point_upper,col="#007D82",lty=3)
lines(s$z,s$estimate,col="#007D82",lwd=2)
lines(grid,1-grid^2,lty=2)
legend("bottom",c("Estimate","Known truth","Pointwise 95% limits","State-grid 95% band (shading)"),
 col=c("#007D82","black","#007D82","#007D8244"),lty=c(1,2,3,1),lwd=c(2,1,1,8),bty="n",cex=.8)
dev.off()

# LAB C: sensitivity to basis size and HAC bandwidth. Compare output at h=4.
baseline <- infer_contrast(fits[[5]])
fixed7 <- infer_contrast(estimate_horizon(4,Jset=7))
longerHAC <- infer_contrast(estimate_horizon(4,bandwidth=8))
sens <- rbind(cbind(specification="AIC / baseline HAC",baseline),
              cbind(specification="Fixed J=7",fixed7),
              cbind(specification="AIC / HAC 8",longerHAC))
print(sens); write.csv(sens,file.path(out,"sensitivity_h4.csv"),row.names=FALSE)
cat("Finished. Inspect the identification assumptions before interpreting these estimates causally.\n")


# Optional course extensions. These are not manuscript Monte Carlo results.
endog <- simulate_panel(state="endogenous")
endog_fit <- estimate_horizon(0,Jset=4:10,data=endog)
write.csv(infer_contrast(endog_fit),file.path(out,"endogenous_state_contrast.csv"),row.names=FALSE)
write.csv(nonlinear_shift_demo(),file.path(out,"nonlinear_shift.csv"),row.names=FALSE)
