# CEClust: executed small-budget workflow.
# The budgets illustrate the API; they are not a validation of search quality.
library(CEClust)
set.seed(12)
X <- cbind(x=c(rnorm(30,-2,.5),rnorm(30,2,.5)),y=rnorm(60,sd=.5))
g <- CECfitBoundGrid(Z=X,C_grid=c(.2,.5,1),lambda_grid=c(.6,1,1.5),
  optimization="em_bank",penalty="none",k0=3,B=3,algo_seed=14,n_cores=1,
  optimizer_control=list(Rmax=3,Ninit=2,Cmax=1,max_iter=30),
  rsi_policy="all",stab_algo_threshold=.85,rsi_threshold=.95,sat_threshold=0)
print(CECsummariseGrid(g))
CECplotGrid(g)
s <- g$summary
q <- s[s$keep,,drop=FALSE]
q <- q[!duplicated(q$partition_id),,drop=FALSE]
print(q[,c("C","lambda","REO","partition_id")])
for(i in seq_len(nrow(q))) {
  CECplotPartition(g,C=q$C[i],lambda=q$lambda[i],x_var="x",y_var="y")
}
g_bic <- CECrescoreGrid(g,penalty="bic")
CECplotGrid(g_bic)
print(cbind(CEC=g$summary$REO,CEC_BIC=g_bic$summary$REO))
stopifnot(identical(g$banks,g_bic$banks))
# For persistent results, explicitly choose separate paths:
# saveRDS(g,"cec_none.rds"); saveRDS(g_bic,"cec_bic.rds")
# vignette("intro",package="CEClust")

# Extend the window in all four directions and refine relevant boundaries.
# Small caps keep this example short; inspect completed and restrictions.
g_adaptive <- CECfitBoundGrid(Z=X, grid_strategy="adaptive",
  optimization="em_bank", penalty="bic", k0=3, B=3, algo_seed=14,
  optimizer_control=list(Rmax=3, Ninit=2, Cmax=10, max_iter=30),
  adaptive_control=list(max_cells=300, max_rounds=3, max_extensions=2))
CECplotGrid(g_adaptive)
print(g_adaptive$adaptive[c("completed", "stop_reason", "restrictions")])
