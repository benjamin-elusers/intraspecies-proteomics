workspace.size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x)object.size(get(x))))
  class(ws) <- "object_size"
  ws
}

stat_box_data <- function(y, upper_limit = max(y) +sd(y)) { #can be modified if needed#
  return(
    data.frame(
      y = 1.15 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round((mean(y)), 6), '\n')
    )
  )
}

asinh_trans <- scales::trans_new(
  "inverse_hyperbolic_sine",
  transform = function(x) {asinh(x)},
  inverse = function(x) {sinh(x)}
)

tri.to.squ<-function(x)
{
  rn <- row.names(x)
  cn <- colnames(x)
  an <- unique(c(cn,rn))
  myval <-  x[!is.na(x)]
  mymat <-  matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
  }
  return(mymat)
}

give.n <- function(x){
  return(c(y = max(x)+sd(x), label = length(x)))
  # experiment with the multiplier to find the perfect position
}
get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
dist_point_line <- function(a, slope, intercept) {
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  return((det(m))/sqrt(sum(v1*v1)))
}


library(ggplot2)
library(gridExtra)
library(ggsci)
library(ggthemes)

theme_black = function(base_size = 12, base_family = "") {

  theme_grey(base_size = base_size, base_family = base_family) %+replace%

    theme(
      # Specify axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border =  element_blank(),
      panel.grid.major = element_line(color = "grey35",linetype = 'dashed'),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")

    )

}

mytheme = theme_clean(base_size = 14) +
  theme(line=element_line(size=1),
        panel.grid.major = element_line(size=1,lineend='round',linetype='22',color='gray50'),
        panel.grid.minor = element_line(size=0.5,lineend='round',linetype='12',color='gray50'),
        axis.ticks.length = unit(1.5,'mm'),
        axis.ticks = element_line(size=0.5),
        legend.position = 'none')



