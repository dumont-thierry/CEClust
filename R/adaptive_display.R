CECadaptive_cell <- function(x,C,lambda) {
  if(length(C)!=1L || length(lambda)!=1L || !is.finite(C) || !is.finite(lambda) || C<=0 || lambda<=0)
    stop('Supply one finite positive C and lambda.')
  s <- x$summary;u <- log(C)
  exact <- which(abs(log(s$C)-u)<1e-12 & abs(s$lambda-lambda)<1e-12)
  if(length(exact))return(exact[1])
  at <- which(u>=s$u0-1e-12 & u<=s$u1+1e-12 & lambda>=s$v0-1e-12 & lambda<=s$v1+1e-12)
  if(!length(at))stop('The requested point is outside the explored adaptive window.',call.=FALSE)
  at[which.min((log(s$C[at])-u)^2+(s$lambda[at]-lambda)^2)]
}

CECplot_adaptive <- function(grid_obj,selected_C=NULL,selected_lambda=NULL,
  stab_algo_threshold=grid_obj$meta$stab_algo_threshold,rsi_threshold=grid_obj$meta$rsi_threshold,
  sat_threshold=grid_obj$meta$sat_threshold,sensibility=.95,grid_summary=NULL,
  main='REO on the adaptive C x lambda diagram',reo_col=cec_grid_color_palette(100),
  bad_fill=grDevices::adjustcolor('grey70',alpha.f=.72),
  bad_hatch_col=grDevices::adjustcolor('grey25',alpha.f=.75),draw_zones=TRUE,show_legend=TRUE) {
  s <- if(is.null(grid_summary))cec_grid_apply_thresholds(grid_obj,stab_algo_threshold,rsi_threshold,sat_threshold) else grid_summary
  old <- graphics::par(c('mar','xpd'));on.exit(graphics::par(old),add=TRUE)
  graphics::par(mar=c(4.2,4.3,2.7,if(show_legend)7 else 1),xpd=FALSE)
  graphics::plot(NA,NA,xlim=range(c(s$u0,s$u1)),ylim=range(c(s$v0,s$v1)),xaxt='n',
    xlab='C (log scale)',ylab=expression(lambda),main=main)
  ticks <- pretty(range(c(s$u0,s$u1)),n=6)
  graphics::axis(1,at=ticks,labels=signif(exp(ticks),3))
  rr <- cec_grid_reo_range(s)
  ix <- pmax(1L,pmin(length(reo_col),1L+floor((length(reo_col)-1L)*(s$REO-rr[1])/diff(rr))))
  col <- ifelse(s$keep,reo_col[ix],bad_fill)
  unmeasured <- !s$measured;col[unmeasured] <- '#E6EBF0'
  quadrant <- s$measurement_status=='quadrant_unmeasured';col[quadrant] <- '#EAC985'
  graphics::rect(s$u0,s$v0,s$u1,s$v1,col=col,border='white',lwd=.35)
  fail <- s$measured & !s$keep_stability
  if(any(fail))graphics::rect(s$u0[fail],s$v0[fail],s$u1[fail],s$v1[fail],
    density=10,angle=45,col=bad_hatch_col,border=NA)
  if(any(s$exclude_saturation)) {
    at <- which(s$exclude_saturation)
    graphics::rect(s$u0[at],s$v0[at],s$u1[at],s$v1[at],density=10,angle=-45,col=bad_hatch_col,border=NA)
  }
  if(draw_zones) {
    e <- grid_obj$adaptive$mesh$edges
    for(k in seq_len(nrow(e))) {
      i <- e[k,1];j <- e[k,2]
      if(!s$keep[i] || !s$keep[j] || s$partition_id[i]==s$partition_id[j])next
      if(partition_similarity(grid_obj$parts[[s$partition_id[i]]]$labels,
          grid_obj$parts[[s$partition_id[j]]]$labels,'match')>=sensibility)next
      a <- s[i,];b <- s[j,]
      if(abs(a$u1-b$u0)<1e-10 || abs(b$u1-a$u0)<1e-10) {
        u <- if(abs(a$u1-b$u0)<1e-10)a$u1 else b$u1
        graphics::segments(u,max(a$v0,b$v0),u,min(a$v1,b$v1),lwd=1.2)
      } else {
        v <- if(abs(a$v1-b$v0)<1e-10)a$v1 else b$v1
        graphics::segments(max(a$u0,b$u0),v,min(a$u1,b$u1),v,lwd=1.2)
      }
    }
  }
  if(!is.null(selected_C)&&!is.null(selected_lambda)) {
    i <- CECadaptive_cell(grid_obj,selected_C,selected_lambda)
    graphics::rect(s$u0[i],s$v0[i],s$u1[i],s$v1[i],border='black',lwd=2)
  }
  if(show_legend) {
    graphics::par(xpd=NA)
    orders <- sort(unique(s$REO[s$keep]));colors <- reo_col[pmax(1L,pmin(length(reo_col),
      1L+floor((length(reo_col)-1L)*(orders-rr[1])/diff(rr))))]
    limits <- graphics::par('usr')
    graphics::legend(x=limits[2]+.02*diff(limits[1:2]),y=limits[4],xjust=0,yjust=1,bty='n',cex=.8,
      legend=c(paste('REO',orders),'Not retained','Not measured',if(any(quadrant))'Quadrant (not measured)'),
      fill=c(colors,bad_fill,'#E6EBF0',if(any(quadrant))'#EAC985'))
  }
  invisible(s)
}
