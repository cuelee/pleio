circlize_plot =  function(factors, ldf, gwas_beta, eta, eta_ci , lmp_input, eta_col = c('#E77E23','#309F86'), link_col = c('#F51929','#2341F6'), link_hex = 3, size_scale = 0.6, study_palette = default_palette){
  
  circos.clear()
  
  rainbow_palete_hex = c('#F44336','#2196F3','#8BC34A','#FF5722','#9C27B0','#00BCD4','#FFEB3B','#9E9E9E','#3F51B5','#4CAF50','#FF9800','#E91E63','#03A9F4','#CDDC39','#795548','#673AB7','#009688','#FFC107','#607D8B')
  default_alpha = format(as.hexmode(floor(0.7*255)),width = 2, upper.case=T)
  default_palette = paste(rainbow_palete_hex, default_alpha, sep='')
  
  col_mat = matrix(study_palette[1:length(factors)], nrow = 1)
  colnames(col_mat) = factors
  
  p=rep(0,length(factors))
  for (i in 1:length(factors)){
    p[i] = 10^(-lmp_input$snp_info[[factors[i]]]$y)
  }
  names(p) = factors
  
  circos.par(start.degree = 0, canvas.xlim = c(-1/size_scale,1/size_scale), gap.degree = 2, cell.padding = c(0,0,0,0), canvas.ylim = c(-1.5,1.5), points.overflow.warning=F)
  circos.initialize(factors = factors, xlim = c(-1, 1))
  
  ##  Track 1
  cex_text = 0.6
  circos.track(factors = factors, x = rep(-0.3,length(factors)), y = rep(0, length(factors)), ylim = c(-0.01,0.01),
               bg.col = NA,
               bg.border = NA,
               track.height = 0.05,
               panel.fun = function(x, y) {
                 circos.text(x = x, y = CELL_META$cell.ylim[2] + uy(2, "mm"), 
                             labels = sapply(CELL_META$sector.index, function(x) paste('',x,sep='')), facing='clockwise', niceFacing =T, adj = c(0,0), cex = cex_text, font = 4 )
                 circos.text(x = x+0.4, y = CELL_META$cell.ylim[2] + uy(2, "mm"), 
                             labels = paste('P =', formatC(p[CELL_META$sector.index],format='E',digit=2),sep='') , facing='clockwise', niceFacing =T,  adj = c(0,0), cex = cex_text, font = 3)
                 
                 circos.text(x = x+0.8, y = CELL_META$cell.ylim[2] + uy(2, "mm"), 
                             labels = paste('BETA: ',formatC(as.numeric(gwas_beta[CELL_META$sector.index]),format='E',digit=1,drop0trailing=T), sep = ''), facing='clockwise', niceFacing =T,  adj = c(0,0), cex = cex_text, font = 3)
                 
                 circos.rect(xleft = -1, xright = 1, ybottom = -0.005, ytop = 0.005, border = NA, col= col_mat[1,CELL_META$sector.index])
               })
  
  ##  Track 2
  ymin = 0
  ymax = 0
  for(fac in factors){
    ymax = max(c(ymax, lmp_input$points[[fac]][,'y']))
  }
  circos.track(ylim = c(ymin, ymax),
               x = rep(0,length(factors)),
               y = rep(0,length(factors)),
               bg.col = NA,
               bg.border = col_mat[, factors],
               bg.lty= 4,
               bg.lwd= 0.5,
               track.height = 0.2,
               panel.fun = function(x, y)
               {
                 circos.segments(x0=-0.9, y0 = 4, x1=0.9, y1=4, lty = 4, lwd = 0.1, col = '#FF000088')
               }
  )
  
  for (fac in factors){
    circos.points(x = lmp_input$points[[fac]][,'x'],
                  y = lmp_input$points[[fac]][,'y'],
                  sector.index = fac,
                  track.index = 2,
                  pch = '.',
                  cex = 0.4,
                  #lwd = 0.3,
                  col = col_mat[, fac]
    )
    circos.points(x = lmp_input$snp_info[[fac]]$x,
                  y = lmp_input$snp_info[[fac]]$y,
                  sector.index = fac,
                  track.index = 2,
                  pch = 20,
                  cex = 0.4,
                  col = '#000000AA'
    )
  }
  
  ##  Track 3
  yl = 1
  yp = 0.4
  ci_lwd=0.5
  circos.par("track.height" = 0.08)
  xl = eta_ci[[1]]/max(abs(eta))
  xu = eta_ci[[2]]/max(abs(eta))
  circos.track(ylim = c(-yl, yl), 
               x = 0.8 * eta/(max(abs(eta))), 
               y=rep(0,length(eta)), 
               bg.col = NA, 
               bg.border = col_mat[, factors],
               bg.lty= 4,
               bg.lwd= 0.5, 
               track.height = 0.08,
               panel.fun = function(x, y)
               {
                 circos.segments(x0=0, y0 = -yl*0.8, x1=0, y1=yl*0.8, lty = 4, lwd = 0.5, col = '#000000AA')
                 if (x>0){
                   circos.rect(xleft = min(x,0), xright = max(x,0), ybottom = -yl*yp, ytop = yl*yp, border = NA, col= eta_col[1])
                 }
                 else if (x<0){
                   circos.rect(xleft = min(x,0), xright = max(x,0), ybottom = -yl*yp, ytop = yl*yp, border = NA, col= eta_col[2])
                 }
                 
                 circos.segments(x0=max(-1, xl[1,CELL_META$sector.index]), y0=0, x1 = min(1, xu[1,CELL_META$sector.index]), y1=0,lwd = ci_lwd, col = '#000000AA')
                 circos.segments(x0=max(-1, xl[1,CELL_META$sector.index]), y0 = - yl*yp*0.4, x1=max(-1, xl[1,CELL_META$sector.index]), y1 = yl*yp*0.6,lwd = ci_lwd, col = '#000000AA')
                 circos.segments(x0= min(1,xu[1,CELL_META$sector.index]), y0 = - yl*yp*0.4, x1= min(1,xu[1,CELL_META$sector.index]), y1 = yl*yp*0.6,lwd = ci_lwd, col = '#000000AA')
               })
  
  ##  Track links
  hr = 0.5
  for(i in 1:nrow(ldf)){
    sign = as.numeric(as.matrix(ldf[i,'f_sign']))
    sector.index1=as.character(ldf[i,'f_factor'])
    point1=as.numeric(as.matrix(ldf[i,c('f_rstart','f_rend')]))
    sector.index2=as.character(ldf[i,'t_factor'])
    point2=as.numeric(as.matrix(ldf[i,c('t_rstart','t_rend')]))
    hex = format(as.hexmode(max(floor((abs(min(sign,1)))^(link_hex)*255-1),0)),width=2, upper.case=T)
    if(sign>0){
      circos.link(sector.index1, point1, sector.index2, point2, col=paste(link_col[1], hex,sep=''), h.ratio=hr)
    }
    if(sign<0){
      circos.link(sector.index1, point1, sector.index2, point2, col=paste(link_col[2], hex,sep=''), h.ratio=hr) 
    }
  }
  circos.clear()
  
  ## Legend
  x_legend = rev(seq(-1,1,0.1))
  y_legend = c((paste(link_col[1],format(as.hexmode(sapply(floor(abs(x_legend[x_legend>0])^link_hex*255-1), function(x) max(x,0))),width=2,upper.case=T),sep='')),'#FFFFFFFF', rev(paste(link_col[2],format(as.hexmode(sapply(floor(abs(x_legend[x_legend>0])^link_hex*255-1),function(x) max(x,0))),width=2, upper.case=T),sep='')))
  
  col_fun_link = colorRamp2(x_legend, sapply(y_legend, EHtoSH))
  
  lgd_links = Legend(at = rev(seq(-1,1,0.5)), labels = rev(seq(-1,1,0.5)), col_fun = col_fun_link, title_position = "topleft", title = bquote(r[g]),labels_gp = gpar(fontsize = 5), direction = "vertical", grid_width = unit(0.3,'cm'))
  
  #draw(lgd_links, x = unit(0.03, "npc"), y = unit(0.2, "npc"), just = c("left", "top"))
  
  col_fun_beta = function(xs){
    res = rep(0,length(xs))
    for (i in 1:length(xs)){
      if (xs[i] >= 0) {
        res[i] = eta_col[1]
      }
      else if(xs[i]<0){
        res[i] = eta_col[2]}
    }
    return(res)
  }
  
  x_legend = eta/max(abs(eta))*0.8
  
  if(min(eta)>0){
    at_eta = c(formatC(as.numeric(max(eta)),format='E',digit=2) , 0) 
  }
  else if(max(eta)<0){
    at_eta = c(formatC(as.numeric(min(eta)),format='E',digit=2), 0) 
  }
  else {
    at_eta = c(formatC(as.numeric(min(eta)),format='E',digit=2), 0,formatC(as.numeric(max(eta)),format='E',digit=2)) 
  }
  
  lgd_eta = Legend(at = round(as.numeric(at_eta),3), labels = at_eta, col_fun = col_fun_beta, title = bquote(eta), labels_gp = gpar(fontsize = 5), grid_width = unit(0.3,'cm'))
  
  lgd_list_horizontal = packLegend(lgd_links, lgd_eta, direction = "horizontal")
  
  draw(lgd_list_horizontal, x = unit(0.09, "npc"), y = unit(0.99, "npc") , just = c("top"))
}

EHtoSH = function(hex_col){
  white = '#FFFFFF'
  out_col = function(f){
    format(as.hexmode(floor(f)),width=2,upper.case=T)
  }
  R3 = as.hexmode(substr(hex_col,2,3))
  G3 = as.hexmode(substr(hex_col,4,5))
  B3 = as.hexmode(substr(hex_col,6,7))
  A1 = as.hexmode(substr(hex_col,8,9))
  R2 = as.hexmode(substr(white,2,3))
  G2 = as.hexmode(substr(white,4,5))
  B2 = as.hexmode(substr(white,6,7))
  #if(A1 == '0'){return(toupper(paste('#',R3,G3,B3,sep='')))}
  if(A1 == as.hexmode('ff')){return(toupper(paste('#',R2,G2,B2,sep='')))}
  op = (as.integer(A1))/255
  R1 = (as.integer(R3)*(op) + as.integer(R2)*(1-op))
  G1 = (as.integer(G3)*(op) + as.integer(G2)*(1-op))
  B1 = (as.integer(B3)*(op) + as.integer(B2)*(1-op))
  return(toupper( paste('#',out_col(R1), out_col(G1),  out_col(B1),sep='')))
} 


generate_link_in = function(m, h2){
  f = rownames(m)
  
  fs = (apply(abs(m), MARGIN = 1, sum)-1)
  h2_mod = h2/max(h2)*0.8
  low = -1*h2_mod
  high = h2_mod
  
  bin = start = end = sign = NULL;
  for (bs in f){
    s = low[bs]
    for (ps in f){
      if(bs==ps) {next}
      else{
        e = high[bs]
        si = m[bs,ps]
        bin = c(bin,paste(bs,ps,sep='-'))
        start = c(start, s)
        end = c(end, e)
        sign = c(sign,si)
      }
    } 
  }
  p = cbind(start,end,sign)
  row.names(p) =bin
  remove(bin,start,end,e,s,bs,ps,sign,si)
  
  lm = done = NULL;
  for (ffac in factors){
    for (tfac in factors[factors %notin% done]){
      if (ffac == tfac ) {next}
      else{
        lm = rbind(lm, c(ffac, p[paste(ffac,tfac,sep='-'),c('start','end','sign')], tfac, p[paste(tfac,ffac,sep='-'),c('start','end','sign')])) 
      }
    }
    done = c(done,ffac)
  }
  colnames(lm) = c('f_factor','f_rstart','f_rend','f_sign','t_factor','t_rstart','t_rend','t_sign')
  ldf = data.frame(lm)
  
  return(ldf)
}


corrMatOrder <- function(
  corr,
  order = c("AOE", "FPC", "hclust", "alphabet"),
  hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single",
                    "average", "mcquitty", "median", "centroid") )
{
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)
  
  switch(order,
         AOE = reorder_using_aoe(corr),
         FPC = reorder_using_fpc(corr),
         hclust = reorder_using_hclust(corr, hclust.method),
         alphabet = sort(rownames(corr))
  )
}

#' Reorder the variables using the angular order of the eigenvectors.
reorder_using_aoe <- function(corr) {
  x.eigen <- eigen(corr)$vectors[, 1:2]
  e1 <- x.eigen[, 1]
  e2 <- x.eigen[, 2]
  alpha <- ifelse(e1 > 0, atan(e2 / e1), atan(e2 / e1) + pi)
  order(alpha) # returned vector
}

#' Reorder the variables using the first principal component.
reorder_using_fpc <- function(corr) {
  x.eigen <- eigen(corr)$vectors[, 1:2]
  e1 <- x.eigen[, 1]
  order(e1) # returned vector
}

#' Reorder the variables using hclust (Hierarchical Clustering).
reorder_using_hclust <- function(corr, hclust.method) {
  hc <- hclust(as.dist(1 - corr), method = hclust.method)
  order.dendrogram(as.dendrogram(hc)) # returned vector
}

# save.image(file = '/Users/cuelee/Dropbox/github/pleio/pleiotropy_plot/R/pleiotropy_plot.Rdata')
