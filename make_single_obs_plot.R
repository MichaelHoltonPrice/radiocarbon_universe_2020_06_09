uncalYear <- 1100 # Years BP
uncalSig  <- 20   # 20 year measurement uncertainty

phi_m <- exp(-uncalYear/ 8033)
sig_m <- uncalSig * exp(-uncalYear / 8033) / 8033

tau_min <- 850
tau_max <- 1050
dtau    <- 5
tau <- seq(tau_min,tau_max,by=dtau)
calibDf = baydem::bd_load_calib_curve("intcal13")
M <- baydem::bd_calc_meas_matrix(tau, phi_m, sig_m, calibDf, T, F)

fprior <- rep(1/(tau_max-tau_min),length(tau))
fpost  <- M[1,] / (sum(M[1,]*dtau))
pdf('single_obs_inf_plot1.pdf',width=12,height=8)
  plot(tau,fprior,type='l',lwd=3,xaxt='n',ylim=c(0,max(fprior,fpost)),xlab='Calendar Year [AD]',ylab='Density')
dev.off()

pdf('single_obs_inf_plot2.pdf',width=12,height=8)
  baydem::bd_vis_calib_curve(tau_min, tau_max, calibDf, xlab = "Calendar Year [AD]", ylab = "Fraction Modern", invertCol = "gray80")
dev.off()

pdf('single_obs_inf_plot3.pdf',width=12,height=8)
  plot(tau,fprior,type='l',lwd=3,xaxt='n',ylim=c(0,max(fprior,fpost)),xlab='Calendar Year [AD]',ylab='Density',col='grey')
  lines(tau,fpost,lwd=3)
dev.off()



