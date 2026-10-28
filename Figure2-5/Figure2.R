# =====================================================================================================
# A. Fig. 2b, 2c, 2d, 2e, 2g, 2h  and  Fig. S6a, S6b, S6c, S6d   --  IMR90 chr2 MERFISH (Su et al. 2020)
# =====================================================================================================
#
# INPUTS:  
#   chromosome2.tsv, chromosome21.tsv  (Zenodo 3928890, IMR90 MERFISH), use the MERFISH_chr2_lite.csv here as a lite example here, which contains 1,2,3,198,546,556,1259,2368,2862 allel.
#   distances_2-100.txt, distances_chr2R{0-50,50-100,...}M_3000.txt (neighbour 250 kb distances),
#   distances_chr2L100-3000.txt  and  distances_chr2L50-3000.txt  (pair distances for >100 Mb and 50-100 Mb),

# KEY PARAMETERS in the code: d3li = 299 nm (LRC), 250 kb bins, genomic separation >= 100 Mb,
# sLRC = LRC present in > 23 traces; traces with > 10000 LRC removed.
# =====================================================================================================

####################################################################################################
# SETUP: load chr2 MERFISH traces, 250-kb bin pairs, distance table
####################################################################################################
library(ggplot2)
library(plotly)
library(dplyr)
library(circlize)
library(htmlwidgets)

IMR_2 <- read.delim("MERFISH_chr2_lite.csv")
print(colnames(IMR_2))
print(dim(IMR_2))

cells = unique(IMR_2$Chromosome.copy.number)
spots = unique(IMR_2$Genomic.coordinate)
print(paste(length(cells),length(spots)))

Co = IMR_2$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_2$Coordinate = Coordinate

distances_2 = read.table(paste0("distances_2-100.txt"),skip = 1,header = FALSE)
print(dim(distances_2))
distances_2 = round(distances_2[,1],0)[-1]
print(length(distances_2))

spots2 = unique(Coordinate)
distances_tmp = matrix(0,ncol = length(spots2),nrow = length(spots2))
colnames(distances_tmp) = spots2
rownames(distances_tmp) = spots2
distances_tmp[lower.tri(distances_tmp)] <- distances_2[1:(length(distances_2)/100)]

st = 0
ed = 0
dst = 0
for (jj in 1:length(colnames(distances_tmp))) {
  st = c(st,rep(colnames(distances_tmp)[jj],nrow(distances_tmp)-jj))
  ed = c(ed, rownames(distances_tmp)[-c(1:jj)])
  # dst = c(dst,distances_tmp[(jj+1):nrow(distances_tmp),jj])
}
st = st[-1]
ed = ed[-1]
# dst = dst[-1]
lner = as.numeric(ed)-as.numeric(st)

tt1 = which(lner>=100000000)
tt2 = which((lner>50000000)&(lner<100000000))


####################################################################################################
# Fig. 2b  3D distance of neighbour (250 kb) / 50-100 Mb / >100 Mb locus pairs, 1st quartile of neighbours
####################################################################################################
## compared between distances 50-100 and >100 MB #####
distances_2m100 = read.table('distances_chr2L100-3000.txt',skip = 1)
distances_2m1000 = data.frame(distances_2m100)
colnames(distances_2m1000) = c('distance','homolog')

distances_2m50 = read.table('distances_chr2L50-3000.txt',skip = 1)
distances_2m500 = data.frame(distances_2m50)
colnames(distances_2m500) = c('distance','homolog')


# Combine them
color_A <- "orange"  # Soft blue
color_B <- "#9b59b6"  # Coral red
color_C <- "darkolivegreen3"  # Coral red
A = as.numeric(distances_2m1000$distance)
# calculate l between 100 MB bins
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between pair of loci >100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_A, border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

B = as.numeric(distances_2m500$distance)
# calculate l between 50 MB bins
quartiles <- quantile(B, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(B, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between pair of loci 50-100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_B, border = "white")
box(lwd = 3)
lines(density(B,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

C = as.vector(neighbor250k)
# calculate l between neibhor bins
quartiles <- quantile(C, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(C, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between neighbor loci,chr2", 
     xlab = "3D dist. (nm)", col = color_C, border = "white")
box(lwd = 3)
lines(density(C,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)
# Calculate the range for both A and B
# x_range <- range(c(A, B),na.rm = T)
# Calculate densities
dens_A <- density(A,na.rm = T)
dens_B <- density(B,na.rm = T)
dens_C <- density(C,na.rm = T)
# Calculate the y-axis limit
y_max <- max(c(dens_A$y, dens_B$y,dens_C$y))
# Create an empty plot
plot(0, 0, type = "n", xlim = c(0,6000), ylim = c(0, y_max * 1.1),
     main = "Overlaid Histograms with Density Curves",
     xlab = "3D dist. (nm)", ylab = "Density")
box(lwd = 3)
# Add subtle grid lines for better readability
grid(col = "gray90", lty = "dotted")
hist(A, freq = FALSE, col = adjustcolor(color_A, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
hist(B, freq = FALSE, col = adjustcolor(color_B, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
hist(C, freq = FALSE, col = adjustcolor(color_C, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
lines(dens_A, col = color_A, lwd = 2)
lines(dens_B, col = color_B, lwd = 2)
lines(dens_C, col = color_C, lwd = 2)
abline(v = quantile(C,c(0.25,0.5,0.75),na.rm=T)[1], col = "darkblue", lty = 2, lwd = 2) 
# Add a legend
legend("topright", legend = c("bin pairs >100M","bin pairs 50-100M", "bin pairs = 250K","1st QTL of the neibhor"),
       col = c(adjustcolor(color_A, alpha.f = 0.5), adjustcolor(color_B, alpha.f = 0.5), adjustcolor(color_C, alpha.f = 0.5),adjustcolor("darkblue", alpha.f = 0.5)),
       lty = c(1,1,1,2), lwd = 2,
       border = NA, bty = "n")

####################################################################################################
# Figure 2C, relationship between convex hull volume v.s. LRC counts ####
####################################################################################################
convex_hull_volume <- function(points) {
  library(geometry)
  hull <- convhulln(points, output.options = "FA")
  return(hull$vol)
}

distances_2m100 = read.table('distances_chr2L100-3000.txt',skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

chr2 = 242200000
rsl = 50
d3li = 299

## distribution of distances of loci pairs ####
conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

chv = rep(NA,length(cells))
for (ii in 1:length(cells)) {
  if (ii%%500==0) {
    print(ii)
  }
  tmp = which(IMR_2$Chromosome.copy.number ==ii)
  set = as.matrix(IMR_2[tmp,1:3])
  set = na.omit(set)
  
  vol <- convex_hull_volume(set)
  chv[ii] = vol
}
rcr = cbind(chv, colSums(conformation_matrix))
colnames(rcr) = c('chv','numbLRC')
plot(rcr[,1],rcr[,2],col = 'darkgrey',pch = '.',cex = 2)
plot(log10(rcr[,1]),log10(rcr[,2]),col = 'darkgrey',pch = '.')
abline(v = 2.5)

# contours plot
data <- data.frame(x = log10(rcr[,1]), y = log10(rcr[,2]))
data %>% ggplot(aes(x, y))+
  geom_point(colour = "brown", size = 0.1)+
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "darkblue",alpha =0.3,bins = 20)+
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  theme_classic()+
  labs(title = "relationship between # of LRC and compaction",
       x = "convex hull volume",
       y = "# LRC/chr.")

## number of count of LRC in individual chromosome ####
hist(log10(colSums(conformation_matrix)), probability = TRUE, breaks = 200,
     main = "count of LRC in individual chr.", 
     xlab = "count of LRC", col = "lightgreen", border = "white",xaxt = "n")  # Suppress default x-axis

##  distribution of individual LRC ####
rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))
print(hist(distances_stat2[,4]))

tgt = which(rowSums(conformation_matrix)>24)
print(length(tgt))
data_df =distances_matrix[tgt,]
data_df2 = data.frame(value=as.vector(data_df), pairs =rep(1:nrow(data_df)+length(tgt),3000))
rdm = 1:nrow(distances_stat2)
rdm = rdm[-tgt]
rdm = sample(rdm,length(tgt))
data_rdm =distances_matrix[rdm,]
data_rdm2 = data.frame(value=as.vector(data_rdm), pairs =rep(1:nrow(data_rdm),3000))
colors = c(rep('random',nrow(data_rdm2)),rep('targets',nrow(data_df2)))
data_df3 = rbind(data_rdm2,data_df2)
data_df3$colors = factor(colors,levels = c('targets','random'))
print(length(tgt))

custom_colors <- c("targets" = "purple", "random" = "lightblue")

# ggplot(data_df3, aes(x = log2(value), group = pairs,color = colors)) +
#   geom_line(stat = "density", size = 0.1, alpha = 0.5) +
#   scale_color_manual(values = custom_colors) +
#   theme_minimal() +
#   labs(title = "Distribution of 1000 Vectors",
#        x = "Value",
#        y = "Density") +
#   xlim(5,13) +  # Limit x-axis from 1 to 100
#   geom_vline(xintercept = log2(299), linetype = "dashed", color = "red")

data_df3t = data_df3[which(data_df3$colors=='random'),]
ggplot(data_df3t, aes(x = (value), group = pairs,color = colors)) +
  geom_line(stat = "density", size = 0.1, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "3D dist. in each chr.",
       x = "3D dist. (nm)",
       y = "Density") +
  xlim(-10,6000) +  # Limit x-axis from 1 to 100
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 299, linetype = "dashed", color = "red")
# ggsave("random.png",width = 947/72, height = 405/72, units = "in", dpi = 300,  bg = "transparent")

data_df3t = data_df3[which(data_df3$colors=='targets'),]
ggplot(data_df3t, aes(x = log2(value), group = pairs,color = colors)) +
  geom_line(stat = "density", size = 0.1, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "3D dist. in each chr.",
       x = "3D dist. (nm)",
       y = "Density") +
  xlim(5,13) +  # Limit x-axis from 1 to 100
  geom_vline(xintercept = 5.1, linetype = "dashed", color = "red")
# ggsave("targets.png",width = 947/72, height = 405/72, units = "in", dpi = 72,  bg = "transparent")
####################################################################################################
# Figure 2E, distribution of number of cells that a LRC exist in #####
####################################################################################################
distances_2m100 = read.table('distances_chr2L100-3000.txt',skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))

data = rowSums(conformation_matrix)
a = boxplot(data,notch = T,col = 'orange',outpch = 19,outcex = 0.5,outcol = 'darkorange3',horizontal = T)
hist(data, breaks = 200, prob = TRUE, col = "orange", xlab = "Values", main = " ")
box(lwd = 3)
abline(v = a$stats[5],lty = 3, col = 'orchid', lwd = 2)

####################################################################################################
# Fig. 2c (linkage map + 3D) for chr2 traces 2368, 2862, 1259; Fig. S6c,d for traces 556 and 546
####################################################################################################
# interactive plots, linkage maps, 3D, networks ####

chr2 = 242200000
rsl = 50
d3li = 299
distances_2m100 = read.table('distances_chr2L100-3000.txt',skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

chrom <- (1:(round((chr2+2000000)/2000000)+1))*2-2
start <- rep(0,(round((chr2+2000000)/2000000)+1))
end <- rep(2000000/(rsl*1000),(round((chr2+2000000)/2000000)+1))
chr_sizes_df <- data.frame(chrom,start,end)

smpl = sample(1:3000,6)
for (ss in smpl) {
  print(ss)
  print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
  distances_stat_tmp = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(distances_matrix[,ss])))
  tmp = rep(0,nrow(distances_stat_tmp))
  tmp[which(as.numeric(as.character(distances_stat_tmp[,4]))<d3li)] = 1
  distances_stat_tmp= cbind(distances_stat_tmp,tmp)
  distances_stat_tmp = distances_stat_tmp[which(as.numeric(as.character(distances_stat_tmp[,4]))<d3li),]
  print(dim(distances_stat_tmp))
  
  ## linkage map ####
  circos.clear()
  circos.par(cell.padding=c(0,0,0,0),
             track.margin=c(0,0.05),
             start.degree = 90,
             gap.degree = 0.1,
             clock.wise = TRUE)
  circos.initialize(factors = chr_sizes_df$chrom,
                    xlim = cbind(chr_sizes_df$start, chr_sizes_df$end))
  
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, facing = 'clockwise')
  },bg.col="#d3e5ed", bg.border=TRUE, track.height=0.1)
  circos.rect(
    sector.index = 244,
    # track.index	= 1,
    xleft = 1,
    xright = 40,
    col = 'white',
    border = NA,
    ytop = 1.5,
    ybottom=-0.5
  )
  circos.rect(
    sector.index = 242,
    # track.index	= 1,
    xleft = 0,
    xright = 40,
    col = 'white',
    border = NA,
    ytop = 1.5,
    ybottom=-0.5
  )
  title(paste0('linkage map: ',ss), cex.main = 1.5)
  locs = seq(0,(chr2+2000000),2000000)
  for (ii in 1:(length(locs)-1)) {
    loc_st = locs[ii]
    loc_ed = locs[ii+1]
    locst = seq(loc_st,loc_ed,rsl*1000)
    tmp1 = which((distances_stat_tmp[,1]>=loc_st) & (distances_stat_tmp[,1]<loc_ed))
    if (length(tmp1)>0) {
      tmp1 = matrix(as.matrix(distances_stat_tmp[tmp1,]),ncol=5)
      for (jj in (ii+1):122) {
        # for (jj in 121:121) {  # create the small links for 240-242 MB locus
        loc_st = locs[jj]
        loc_ed = locs[jj+1]
        loced = seq(loc_st,loc_ed,rsl*1000)
        tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
        if(length(tmp2)!=0){
          tmp2 = matrix(tmp1[tmp2,],ncol=5)
          for (kk in 1:length(locst)) {
            tmpp = which((tmp2[,1]-1)==locst[kk])
            if(length(tmpp)!=0){
              tmpp = matrix(tmp2[tmpp,],ncol=5)
              for (kkk in 1:length(loced)) {
                tmppp = which((tmpp[,2]-1)==loced[kkk])
                if(length(tmppp)!=0){
                  lwdd = (tmpp[tmppp,4]/1000)
                  if(tmppp!=0){
                    circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'blue',lwd = lwdd)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  ## 3D plot ####
  
  # Create sample data (replace this with your actual data)
  tmp = which(IMR_2$Chromosome.copy.number==ss)
  all_colors = colorspace::rainbow_hcl(length(tmp)+length(tmp))
  set.seed(123)
  data <- data.frame(
    index = 1:length(tmp),
    x = IMR_2$X.nm.[tmp],
    y = IMR_2$Y.nm.[tmp],
    z = IMR_2$Z.nm.[tmp],
    cord = (IMR_2$Coordinate[tmp]-1)/1000000
  )
  
  # Create a group column in the data
  data$group <- "Other"
  distance_group_indices <- unique(c(
    which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp[,1]),
    which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp[,2])
  ))
  data$group[distance_group_indices] <- "sLRC"
  all_colors[distance_group_indices] = rainbow(length(tmp)+length(tmp))[distance_group_indices]
  
  ## Create the base plot for html ####
  plot <- plot_ly()
  
  # Add all points
  plot <- add_trace(plot,
                    data = subset(data), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers+lines",
                    marker = list(size = 7, color = all_colors, symbol = "circle"),
                    line = list(color = "grey", width = 1),
                    name = "Other",
                    visible = TRUE,
                    text = data$cord,  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  
  # Add sLRC points
  plot <- add_trace(plot,
                    data = subset(data, group == "sLRC"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 7, color = all_colors[distance_group_indices], symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "sLRC",
                    legendgroup = "sLRC",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="sLRC")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  # Add DistanceGroup links
  for (pp in 1:nrow(distances_stat_tmp)) {
    pp1 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp[pp,1])
    pp2 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp[pp,2])
    
    plot <- add_trace(plot,
                      x = c(data$x[pp1], data$x[pp2]),
                      y = c(data$y[pp1], data$y[pp2]),
                      z = c(data$z[pp1], data$z[pp2]),
                      type = "scatter3d",
                      mode = "lines",
                      line = list(color = "black", width = 3),
                      name = "sLRC",
                      legendgroup = "sLRC",
                      showlegend = F,
                      visible = TRUE)
  }
  
  # Add scale bar
  scale_length <- 1000
  min_x <- min(data$x)
  min_y <- min(data$y)
  min_z <- min(data$z)
  
  plot <- add_trace(plot,
                    x = c(min_x, min_x + scale_length),
                    y = c(min_y, min_y),
                    z = c(min_z, min_z),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  plot <- add_trace(plot,
                    x = c(min_x, min_x),
                    y = c(min_y, min_y + scale_length),
                    z = c(min_z, min_z),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  plot <- add_trace(plot,
                    x = c(min_x, min_x),
                    y = c(min_y, min_y),
                    z = c(min_z, min_z + scale_length),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  
  # Set up layout
  plot <- plot %>% layout(
    scene = list(
      aspectmode = "cube",
      xaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = ""),
      yaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = ""),
      zaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = "")
    ),
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        xanchor = "center",
        yanchor = "top",
        pad = list(r = 0, t = 10, b = 10),
        x = 0.5,
        y = 1.1,
        buttons = list(
          list(method = "update",
               args = list(visible = list(TRUE, TRUE)),
               label = "All"),
          list(method = "update",
               args = list(visible = list(TRUE, FALSE)),
               label = "Distance Group"),
          list(method = "update",
               args = list(visible = list(FALSE, TRUE)),
               label = "Other Group")
        )
      )
    )
  )
  print(plot)
  
  plot <- plot %>% layout(
    scene = list(
      camera = list(
        up = list(x = 0, y = 0, z = 1),
        center = list(x = 0, y = 0, z = 0),
        eye = list(x = 1.25, y = 1.25, z = 1.25)
      ),
      dragmode = "turntable"
    )
  )%>%onRender("
    function(el, x) {
      var id = el.getAttribute('id');
      var gd = document.getElementById(id);
      
      var time = 0;
      function rotate() {
        time += 0.01;
        
        // Create a spiral rotation
        var radius = 1.25 + 0.25 * Math.sin(time * 0.5);  // Varying radius
        Plotly.relayout(gd, {
          'scene.camera.eye.x': radius * Math.cos(time),
          'scene.camera.eye.y': radius * Math.sin(time),
          'scene.camera.eye.z': 1.25 + 1.25 * Math.cos(time)  // Up and down motion
        });
        requestAnimationFrame(rotate);
      }
      rotate();
    }
  ")
  print(plot)
  # Save the plot as an interactive HTML file
  # htmlwidgets::saveWidget(plot, paste0("interactive_",ss,".html"))
  # print( paste0("interactive_",ss,".html"))
  # 
  # ## network map ####
  # links = data.frame(from = (distances_stat_tmp[,1]-1)/1000000, 
  #                    to = (distances_stat_tmp[,2]-1)/1000000, 
  #                    weight = 100/distances_stat_tmp[,4],
  #                    type = rep('lks',nrow(distances_stat_tmp)))
  # nodes = data.frame(id = unique(c((distances_stat_tmp[,1]-1)/1000000,(distances_stat_tmp[,2]-1)/1000000)),
  #                    audience.size = rep(1,length(unique(c(distances_stat_tmp[,1],distances_stat_tmp[,2])))))
  # print(head(distances_stat_tmp))
  # print(head(links))
  # print(head(nodes))
  # 
  # net <- graph_from_data_frame(links, nodes, directed=T)
  # tkid <- tkplot(net) 
  # V(net)$color <- 'white'
  # V(net)$frame.color<- 'purple'
  # # Compute node degrees (#links) and use that to set node size:
  # deg <- igraph::degree(net, mode="all")
  # V(net)$size <- deg/3
  # # V(net)$label <- ""
  # # V(net)$size <- V(net)$weight*
  # # Set edge width based on weight:
  # E(net)$width <- E(net)$weight
  # l = layout.fruchterman.reingold
  # plot(net,vertex.label.cex = deg/50, edge.arrow.size=.04, edge.color="orange",layout=l)
  # 
  # cut.off <- 1
  # net.sp <- igraph::delete_vertices(net, V(net)[size<cut.off])
  # deg = deg[which(V(net)$size>=cut.off)]
  # print(c(length(V(net.sp)$size), length(deg)))
  # cut.off <- mean(links$weight)
  # net.sp <- igraph::delete_edges(net.sp, E(net.sp)[weight<cut.off])
  # 
  # print(c(length(V(net.sp)$size), length(deg)))
  # 
  # l = layout.fruchterman.reingold
  # plot(net.sp,vertex.label.cex = deg/30,edge.arrow.size=.04, edge.color="orange",layout=l)
  # 
  # tkid <- tkplot(net) #tkid is the id of the tkplot that will open
  # l <- tk_coords(tkid) # grab the coordinates from tkplot
  # deg <- igraph::degree(net, mode="all")
  # V(net)$size <- deg
  # vertex_colors = match(V(net),seq(0,242,0.25))
  # vertex_colors <- rainbow(length(seq(0,242,0.25)*2))[vertex_colors]
  # 
  # plot(net, layout=l,vertex.label.cex = 0.5,
  #      edge.arrow.size=.04, edge.color="orange",
  #      # mark.groups=list(c(1), c(length(V(net.sp)))),
  #      mark.col= vertex_colors,
  #      # mark.border=NA
  #      # vertex.label.color = vertex_colors
  # )
}

####################################################################################################
# Fig. 2h  overlay of MERFISH sLRC (DNA FISH) and Hi-C LRC (CLOD outliers) on chr2
# NOTE: a MERFISH sLRC and a Hi-C LRC count as overlaid if both ends are within 250 kb.
####################################################################################################
# Overlay of DNA FISH and HiC ####
library(circlize)
chr2 = 242200000
chr_l = chr2
chr = 2
rsl = 100
d3li = 299

## input HiC ####
flss = "IMR90"
fls = paste0("outliers4_",flss,'-',rsl,"k-chr",chr,".csv")
print(fls)
hic_IQRs = read.delim(fls,header = T,sep = ',')
head(hic_IQRs)

distances_stat_hic = data.frame(hic_IQRs[,2:5])
tmp = which((distances_stat_hic$loc2-distances_stat_hic$loc1)>100000000)
distances_stat_hic = distances_stat_hic[tmp,]
tmp = which((distances_stat_hic[,4])>2)
distances_stat_hic = distances_stat_hic[tmp,]
print(dim(distances_stat_hic))

## input DNA FISH ####
IMR_2 <- read.delim("chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))

cells = unique(IMR_2$Chromosome.copy.number)
spots = unique(IMR_2$Genomic.coordinate)
print(paste(length(cells),length(spots)))

Co = IMR_2$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_2$Coordinate = Coordinate

distances_2 = read.table(paste0("distances_2-100.txt"),skip = 1,header = FALSE)
dim(distances_2)
distances_2 = round(distances_2[,1],0)[-1]
print(length(distances_2))

spots2 = unique(Coordinate)
distances_tmp = matrix(0,ncol = length(spots2),nrow = length(spots2))
colnames(distances_tmp) = spots2
rownames(distances_tmp) = spots2
distances_tmp[lower.tri(distances_tmp)] <- distances_2[1:(length(distances_2)/100)]

st = 0
ed = 0
dst = 0
for (jj in 1:length(colnames(distances_tmp))) {
  st = c(st,rep(colnames(distances_tmp)[jj],nrow(distances_tmp)-jj))
  ed = c(ed, rownames(distances_tmp)[-c(1:jj)])
  # dst = c(dst,distances_tmp[(jj+1):nrow(distances_tmp),jj])
}
st = st[-1]
ed = ed[-1]
# dst = dst[-1]
lner = as.numeric(ed)-as.numeric(st)

tt1 = which(lner>=100000000)
tt2 = which((lner>50000000)&(lner<100000000))

distances_2m100 = read.table(paste0('distances_chr2L100-3000.txt'),skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))

tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]

## overlaid DNA FISH and HiC ####
head(distances_stat22)
head(distances_stat_hic)
distances_stat22[,1] = distances_stat22[,1]-1
distances_stat22[,2] = distances_stat22[,2]-1

label_hic = rep(0,nrow(distances_stat_hic))
for (ii in 1:nrow(distances_stat_hic)) {
  tmp = as.numeric(distances_stat_hic[ii,])
  st1 = which(abs((spots2-1)-tmp[1])<=250000)
  ed1 = which(abs((spots2-1)-tmp[2])<=250000)
  if (length(st1)==0 | length(ed1)==0 ) {
    next
  }
  st2 = which(abs(distances_stat22[,1]-tmp[1])<=250000)
  ed2 = which(abs(distances_stat22[,2]-tmp[2])<=250000)
  p2= intersect(st2,ed2)
  if(length(p2)!=0){
    label_hic[ii] = 2
  }else{
    label_hic[ii] = 1
    print(ii)
  }
}

label_hic = rep(0,nrow(distances_stat22))
for (ii in 1:nrow(distances_stat22)) {
  tmp = as.numeric(distances_stat22[ii,])
  st2 = which(abs(distances_stat_hic[,1]-tmp[1])<250000)
  ed2 = which(abs(distances_stat_hic[,2]-tmp[2])<250000)
  p2= intersect(st2,ed2)
  if(length(p2)!=0){
    label_hic[ii] = 2
  }else{
    label_hic[ii] = 1
    print(ii)
  }
}

# plot basic ring ####
rsl = 50
chrom <- (1:(round((chr_l+2000000)/2000000)+1))*2-2
start <- rep(0,(round((chr_l+2000000)/2000000)+1))
end <- rep(2000000/(rsl*1000),(round((chr_l+2000000)/2000000)+1))
chr_sizes_df <- data.frame(chrom,start,end)

circos.clear()
circos.par(cell.padding=c(0,0,0,0),
           track.margin=c(0,0.05),
           start.degree = 90,
           gap.degree = 0.1,
           clock.wise = TRUE)
circos.initialize(factors = chr_sizes_df$chrom,
                  xlim = cbind(chr_sizes_df$start, chr_sizes_df$end))
transparent_color <- adjustcolor('lightblue', alpha.f = 0.4)
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex = 1, facing = 'clockwise')
},bg.col=transparent_color, bg.border=TRUE, track.height=0.1)
circos.rect(
  sector.index = 244,
  # track.index	= 1,
  xleft = 1,
  xright = 8,
  col = 'white',
  border = NA,
  ytop = 1.5,
  ybottom=-0.5
)

title(paste(flss,'Overlay'), cex.main = 1.5)
#plot MERFISH 
locs = seq(0,(chr2+2000000),2000000)+1
for (ii in 1:(length(locs)-1)) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,rsl*1000)
  tmp1 = which((distances_stat22[,1]>=loc_st) & (distances_stat22[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(distances_stat22[tmp1,]),ncol=4)
    for (jj in (ii+1):((round((chr2)/2000000)))) {
      # for (jj in 122:122) {  # create the small links for 240-242 MB locus
      loc_st = locs[jj]
      loc_ed = locs[jj+1]
      loced = seq(loc_st,loc_ed,rsl*1000)
      tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<loc_ed))
      if(length(tmp2)!=0){
        tmp2 = matrix(tmp1[tmp2,],ncol=4)
        for (kk in 1:length(locst)) {
          tmpp = which((tmp2[,1])==locst[kk])
          if(length(tmpp)!=0){
            tmpp = matrix(tmp2[tmpp,],ncol=4)
            for (kkk in 1:length(loced)) {
              tmppp = which((tmpp[,2])==loced[kkk])
              if(length(tmppp)!=0){
                lwdd = (tmpp[tmppp,4]/1000)
                if(length(tmppp)!=0){
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'deepskyblue3',lwd = lwdd)
                }
              }
            }
          }
        }
      }
    }
  }
}

#plot HiC  
distances_stat_hic2 = distances_stat_hic[which(label_hic==0),]
for (ii in 1:(round((chr_l+2000000)/2000000))) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,rsl*1000)
  tmp1 = which((distances_stat_hic2[,1]>=loc_st) & (distances_stat_hic2[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(distances_stat_hic2[tmp1,]),ncol=4)
    for (jj in (ii+1):((round((chr_l)/2000000)))) {# create the small links for all locus but 242MB
      # for (jj in 122:122) {
      loc_st = locs[jj]
      loc_ed = locs[jj+1]
      loced = seq(loc_st,loc_ed,rsl*1000)
      tmp2 = which((tmp1[,2]>=loc_st) & ((tmp1[,2])<loc_ed))
      if(length(tmp2)!=0){
        tmp2 = matrix(tmp1[tmp2,],ncol=4)
        for (kk in 1:length(locst)) {
          tmpp = which((tmp2[,1])==locst[kk])
          if(length(tmpp)!=0){
            tmpp = matrix(tmp2[tmpp,],ncol=4)
            for (kkk in 1:length(loced)) {
              tmppp = which((tmpp[,2])==loced[kkk])
              if(length(tmppp)!=0){
                tmppp = as.numeric(tmpp[tmppp,3])
                if(tmppp!=0){
                  link_thickness <- tmppp/20
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'darkgrey', lwd = 0.1)
                }
              }
            }
          }
        }
      }
    }  
  }
}

#plot HiC unique line
distances_stat_hic2 = distances_stat_hic[which(label_hic==1),]
for (ii in 1:(round((chr_l+2000000)/2000000))) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,rsl*1000)
  tmp1 = which((distances_stat_hic2[,1]>=loc_st) & (distances_stat_hic2[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(distances_stat_hic2[tmp1,]),ncol=4)
    for (jj in (ii+1):((round((chr_l)/2000000)))) {# create the small links for all locus but 242MB
      # for (jj in 122:122) {
      loc_st = locs[jj]
      loc_ed = locs[jj+1]
      loced = seq(loc_st,loc_ed,rsl*1000)
      tmp2 = which((tmp1[,2]>=loc_st) & ((tmp1[,2])<loc_ed))
      if(length(tmp2)!=0){
        tmp2 = matrix(tmp1[tmp2,],ncol=4)
        for (kk in 1:length(locst)) {
          tmpp = which((tmp2[,1])==locst[kk])
          if(length(tmpp)!=0){
            tmpp = matrix(tmp2[tmpp,],ncol=4)
            for (kkk in 1:length(loced)) {
              tmppp = which((tmpp[,2])==loced[kkk])
              if(length(tmppp)!=0){
                tmppp = as.numeric(tmpp[tmppp,3])
                if(tmppp!=0){
                  link_thickness <- tmppp/20
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'darkolivegreen4', lwd = 0.1)
                }
              }
            }
          }
        }
      }
    }  
  }
}

#plot overlaid
distances_stat_hic2 = distances_stat_hic[which(label_hic==2),]
for (ii in 1:(round((chr_l+2000000)/2000000))) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,rsl*1000)
  tmp1 = which((distances_stat_hic2[,1]>=loc_st) & (distances_stat_hic2[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(distances_stat_hic2[tmp1,]),ncol=4)
    for (jj in (ii+1):((round((chr_l)/2000000)))) {# create the small links for all locus but 242MB
      # for (jj in 122:122) {
      loc_st = locs[jj]
      loc_ed = locs[jj+1]
      loced = seq(loc_st,loc_ed,rsl*1000)
      tmp2 = which((tmp1[,2]>=loc_st) & ((tmp1[,2])<loc_ed))
      if(length(tmp2)!=0){
        tmp2 = matrix(tmp1[tmp2,],ncol=4)
        for (kk in 1:length(locst)) {
          tmpp = which((tmp2[,1])==locst[kk])
          if(length(tmpp)!=0){
            tmpp = matrix(tmp2[tmpp,],ncol=4)
            for (kkk in 1:length(loced)) {
              tmppp = which((tmpp[,2])==loced[kkk])
              if(length(tmppp)!=0){
                tmppp = as.numeric(tmpp[tmppp,3])
                if(tmppp!=0){
                  link_thickness <- tmppp/20
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'brown', lwd = 1)
                }
              }
            }
          }
        }
      }
    }  
  }
}

####################################################################################################
# Fig. S6a  chr21:32.45-33.35 Mb, distance between adjacent 50-kb loci (median line)
####################################################################################################
# genomic distances (l) between pair of bins####
## Fig. S6B ####
IMR_21 <- read.delim("chromosome21.tsv")
print(head(IMR_21))
print(colnames(IMR_21))
print(dim(IMR_21))

spots = unique(IMR_21$Genomic.coordinate)

Co = IMR_21$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_21$Coordinate = Coordinate

# chr21 32.45 MB - 33.35 MB, Fig. S1B
tmp = which(Coordinate>32450000&Coordinate<33350000)
IMR_211 = IMR_21[tmp,]

neighbor50k2 = 0
for (ii in 2:nrow(IMR_211)) {
  if (ii%%100000==0) {
    print(ii)
  }
  if ((IMR_211$Coordinate[ii]-IMR_211$Coordinate[ii-1])==50000) {
    tmp = ((IMR_211$Z.nm.[ii-1]-IMR_211$Z.nm.[ii])^2+(IMR_211$X.nm.[ii-1]-IMR_211$X.nm.[ii])^2+(IMR_211$Y.nm.[ii-1]-IMR_211$Y.nm.[ii])^2)^0.5
    neighbor50k2 = c(neighbor50k2,as.integer(tmp))
  }
}
neighbor50k2 = neighbor50k2[-1]
print(median(neighbor50k2,na.rm = T))

A = as.vector(neighbor50k2)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
print(quartiles)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 400, xlim = c(1,600), 
     main = "chr21:32.45 MB-33.35 MB", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c( "Median"),
       col = c( "brown3"),
       lty = c( 2), lwd = 2)

####################################################################################################
# Fig. S6b  chr21 (all) and chr2 (250 kb) neighbor-locus distance distributions
####################################################################################################
neighbor50k = 0
for (ii in 2:nrow(IMR_21)) {
  if (ii%%100000==0) {
    print(ii)
  }
  if (IMR_21$Coordinate[ii]-IMR_21$Coordinate[ii-1]==50000) {
    tmp = ((IMR_21$Z.nm.[ii-1]-IMR_21$Z.nm.[ii])^2+(IMR_21$X.nm.[ii-1]-IMR_21$X.nm.[ii])^2+(IMR_21$Y.nm.[ii-1]-IMR_21$Y.nm.[ii])^2)^0.5
    # if (is.na(tmp)) {
    #   print(ii)
    # }
    neighbor50k = c(neighbor50k,as.integer(tmp))
  }
}
neighbor50k = neighbor50k[-1]
print(median(neighbor50k,na.rm = T))

A = as.vector(neighbor50k)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
print(quartiles)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "chr21", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

IMR_2 <- read.delim("chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))

cells = unique(IMR_2$Chromosome.copy.number)
spots = unique(IMR_2$Genomic.coordinate)
print(paste(length(cells),length(spots)))

Co = IMR_2$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_2$Coordinate = Coordinate

## the 3D distances between pairs of neibhor 250k ####
neighbor250k = matrix(0,nrow = 1000, ncol = 3000)
aa = 1
for (ll in 1:5) {
  tt = which((spots2>(ll-1)*50000000)&(spots2<(ll*50000000+1)))
  distances_2m_tt = read.table(paste0('distances_chr2R',(ll-1)*50,'-',(ll*50),'M_3000.txt'))
  print(dim(distances_2m_tt))
  colnames(distances_2m_tt) = distances_2m_tt[1,]
  distances_2m_tt = t(distances_2m_tt[-1,])
  print(dim(distances_2m_tt))
  if((dim(distances_2m_tt)[2]/3000)!=dim(distances_2m_tt)[1]){
    print('WRONG!!!')
    break
  }
  bb = aa+dim(distances_2m_tt)[1]-2
  for (cells in 1:3000) {
    start_col <- (cells - 1) * dim(distances_2m_tt)[1] + 1
    end_col <- (cells) * dim(distances_2m_tt)[1]
    xx<- distances_2m_tt[, start_col:end_col]
    neighbor250k[aa:bb,cells] = xx[row(xx) == col(xx) - 1]
  }
  aa = bb+1
  print(c(bb,aa))
}
print(aa)
neighbor250k = neighbor250k[1:bb,]
# # double check
# neighbor250k2 = 0
# for (ii in 914186:nrow(IMR_2)) {
#   if (ii%%100000==0) {
#     print(ii)
#   }
#   if (IMR_2$Coordinate[ii]-IMR_2$Coordinate[ii-1]==250000) {
#     tmp = ((IMR_2$Z.nm.[ii-1]-IMR_2$Z.nm.[ii])^2+(IMR_2$X.nm.[ii-1]-IMR_2$X.nm.[ii])^2+(IMR_2$Y.nm.[ii-1]-IMR_2$Y.nm.[ii])^2)^0.5
#     neighbor250k2 = c(neighbor250k2,as.integer(tmp))
#   }
# }
# neighbor250k2 = neighbor250k2[-1]

print(median(neighbor250k,na.rm = T))
A = as.vector(neighbor250k)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "chr2", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)
