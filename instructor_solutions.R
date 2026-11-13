# Run from the course root. Numerical answers for the 30-minute practical.
source("code/R/exercise_setup.R")
impact <- estimate_horizon(0,Jset=4:10)
grid <- seq(-.8,.8,length.out=161)
curve <- infer_curve(impact,grid)
contrast <- infer_contrast(impact)
fixed <- estimate_horizon(0,Jset=7)
m <- lm(y ~ x + I(x*z) + ylag,data=d)
glin <- coef(m)["x"]+coef(m)["I(x * z)"]*grid
core <- rbind(cbind(specification="AIC J=4:10",J=impact$J,contrast),
              cbind(specification="Fixed J=7",J=7,infer_contrast(fixed)))
write.csv(core,file.path(out,"core_answers.csv"),row.names=FALSE)
write.csv(curve,file.path(out,"core_curve.csv"),row.names=FALSE)
write.csv(data.frame(z=grid,linear=glin,truth=1-grid^2),
          file.path(out,"linear_benchmark.csv"),row.names=FALSE)
write.csv(data.frame(rows=nrow(d),units=length(unique(d$id)),
  dates=length(unique(d$t)),J=impact$J,L=impact$L,
  true_contrast=.75^2,true_h4_contrast=.6^4*.75^2,
  linear_intercept=coef(m)["x"],linear_slope=coef(m)["I(x * z)"]),
  file.path(out,"core_summary.csv"),row.names=FALSE)
pdf(file.path(out,"core_comparison.pdf"),width=8,height=5)
plot(grid,curve$estimate,type="n",ylim=range(curve$band_lower,curve$band_upper,glin),
     xlab="Pre-shock state",ylab="Impact response per unit shock",
     main="Synthetic course exercise: impact response")
polygon(c(grid,rev(grid)),c(curve$band_lower,rev(curve$band_upper)),
        col="#007D8225",border=NA)
lines(grid,curve$estimate,col="#007D82",lwd=2)
lines(grid,1-grid^2,lty=2,lwd=2)
lines(grid,glin,col="#C66B35",lwd=2)
legend("bottom",c("Spline estimate","Known truth","Linear interaction",
                  "95% state-grid band (shading)"),
       col=c("#007D82","black","#C66B35","#007D8244"),
       lty=c(1,2,1,1),lwd=c(2,2,2,8),bty="n",cex=.8)
dev.off()
print(core)
print(curve[which.min(abs(grid)),])

