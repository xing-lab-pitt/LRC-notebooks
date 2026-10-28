# =====================================================================================================
# A. Fig. 2b, 2c, 2d, 2e, 2g, 2h  and  Fig. S6a, S6b, S6c, S6d   --  IMR90 chr2 MERFISH (Su et al. 2020)
# =====================================================================================================
#
# INPUTS:  
#   chromosome2.tsv, chromosome21.tsv  (Zenodo 3928890, IMR90 MERFISH),
#   distances_2-100.txt, distances_chr2R{0-50,50-100,...}M_3000.txt (neighbour 250 kb distances),
#   distances_chr2L100-3000.txt  and  distances_chr2L50-3000.txt  (pair distances for >100 Mb and 50-100 Mb),
#   outliers4_IMR90-100k-chr2.csv  (Hi-C CLOD output used in Fig. 2h; produced by the CLOD code,
#   micro-c_GSM7661135_IMR90.chr2.output.2k.txt'
#   genomic-scale-with transcription and nuclear bodies.tsv     3 Mb tracing data (Su et al. 2020 Cell; chr1-7) 
#   41586_2025_8838_MOESM5_ESM.csv                              probe table of the NMuMG DNA seqFISH+ paper 
#   NMuMG_rep1_pos0..pos9.csv                                   DNA seqFISH+ spots (Zenodo 7693825)
# =====================================================================================================

############################################################
# Fig. 2b  3D distance of neighbor (250 kb) / 50-100 Mb / >100 Mb locus pairs, 1st quartile of neighbors ####
############################################################
## compared with distances 50-100 and >100 MB #####
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
# Fig. 2c; ####
# Fig. S6c,d f####
####################################################################################################
## interactive plots, linkage maps, 3D, networks ####

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

smpl = c(2368,2862,1259)   # representative examples for Fig. 2c. and Fig. S6c,d 
for (ss in smpl) {
  print(ss)
  # print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
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
  ## # network map ####
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


############################################################
# Fig. 2d (convex hull volume vs #LRC per trace), ####
# Fig. 2g (distance distributions sLRC vs random pairs),#### 
# Fig. 2e (histogram + boxplot of #traces per LRC)  ####
####################################################################################################
# Figure 2d, relationship between convex hull volume v.s. LRC counts ####
convex_hull_volume <- function(points) {
  library(geometry)
  hull <- convhulln(points, output.options = "FA")
  return(hull$vol)
}

distances_2m100 = read.table(paste0('distances_chr2L100-3000.txt'),skip = 1)
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

#  Fig. 2g distribution of individual LRC ####
rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))
hist(distances_stat2[,4])

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

# Figure 2e, distribution of number of cells that a LRC exists in #####
distances_2m100 = read.table(paste0('distances_chr2L100-3000.txt'),skip = 1)
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
# Fig. 2h  overlay of MERFISH sLRC (DNA FISH) and Hi-C LRC (CLOD outliers) on chr2 ####
####################################################################################################
# Overlay of DNA FISH and HiC ####
library(circlize)
chr2 = 242200000
chr_l = chr2
chr = 2
rsl = 100
d3li = 299

## input HiC ####
fls = paste0("outliers4_IMR90-",rsl,"k-chr",chr,".csv")
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
IMR_2 <- read.delim("Verify/chromosome2.tsv")
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

distances_2 = read.table("Verify/distances_2-100.txt",skip = 1,header = FALSE)
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

distances_2m100 = read.table('distances_chr2L100-3000.txt',skip = 1)
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
#plot DNA FISH 
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
# Fig. 2i  3D structure + LRC linkage map of single NMuMG chr1 traces (trace # 22, 20, 10 in the figure) ####
####################################################################################################
## plot ####
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
sample_cell = c(22,20,10)   # example
for (ff in 1:length(file_names)) {
  file_name = file_names[ff]
  fls = paste0(file_names[ff],'.csv')
  print(fls)
  seqFISH_DNA <- read.csv(fls)
  split_list <- strsplit(seqFISH_DNA$name, "-")
  loci <- do.call(rbind, split_list)
  loci[, 2] <- as.numeric(loci[, 2])
  seqFISH_DNA2 = data.frame(
    cbind(
      seqFISH_DNA$cellID,seqFISH_DNA$chrom,loci[,2],
      seqFISH_DNA$x_um, seqFISH_DNA$y_um,  seqFISH_DNA$z_um,
      seqFISH_DNA$dbscan_ldp_nbr_allele
    )
  )
  colnames(seqFISH_DNA2) = c('cell', 'chr', 'loci','x_um','y_um','z_um', 'allele')
  cells = seqFISH_DNA2$cell%>%unique()%>%as.character()%>%as.numeric()%>%sort()
  print(paste('#cell:',length(cells)))
  
  for (cc in 'chr1') {   # example
    probes_tmp = which(probes$chrom==cc)
    probes_tmp = probes[probes_tmp,c(1,2,5)]
    
    seqFISH_DNA2_tmp = seqFISH_DNA2[which(seqFISH_DNA2$chr==cc),]
    print(cc)
    print(dim(seqFISH_DNA2_tmp))
    spots2 = unique(seqFISH_DNA2_tmp$loci)%>%as.character()%>%as.numeric()%>%sort()
    print(length(spots2))
    
    for (allel in c(1,2)) {
      fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
      print(fls)
      seqFISH_dist <- read.csv(fls)
      print(dim(seqFISH_dist))
      
      for (cell in sample_cell) {
        tmp = which(seqFISH_DNA2_tmp$cell==cell)
        seqFISH_DNA2_tmp3 = seqFISH_DNA2_tmp[tmp,]
        
        tmp = which(seqFISH_DNA2_tmp3$allele==allel)
        seqFISH_DNA2_tmp2 = seqFISH_DNA2_tmp3[tmp,]
        plot(seqFISH_DNA2_tmp2$z_um)
        tmp = match(probes_tmp$probe_ID,seqFISH_DNA2_tmp2$loci)
        all_colors = colorspace::rainbow_hcl(length(tmp)+length(tmp))
        data <- data.frame(
          index = 1:length(tmp),
          x = seqFISH_DNA2_tmp2$x_um[tmp]%>%as.character()%>%as.numeric(),
          y = seqFISH_DNA2_tmp2$y_um[tmp]%>%as.character()%>%as.numeric(),
          z = seqFISH_DNA2_tmp2$z_um[tmp]%>%as.character()%>%as.numeric(),
          cord = (probes_tmp$start%>%as.character()%>%as.numeric())/1000000
        )
        
        seqFISH_dist_LRC = seqFISH_dist[,c(1:2,cell+2)]
        LRC_tmp = which(seqFISH_dist_LRC[,3]<299)
        seqFISH_dist_LRC = seqFISH_dist_LRC[LRC_tmp,]
        # Create a group column in the data
        data$group <- "Other"
        distance_group_indices <- unique(c(
          which(data$cord %in% (seqFISH_dist_LRC[,1]/1000)),
          which(data$cord %in% (seqFISH_dist_LRC[,2]/1000))
        ))
        data$group[distance_group_indices] <- "LRC"
        all_colors[distance_group_indices] = rainbow(length(tmp)+length(tmp))[distance_group_indices]
        
        ## Create the base plot for html ####
        plot <- plot_ly()
        
        # Add all points
        plot <- add_trace(plot,
                          data = subset(data),
                          x = ~x, y = ~y, z = ~z,
                          type = "scatter3d", mode = "markers+lines",
                          marker = list(size = 3, color = all_colors, symbol = "circle"),
                          line = list(color = "grey", width = 1),
                          name = "Other",
                          visible = TRUE,
                          text = data$cord,  # Use the index as the text
                          hoverinfo = "text",    # Show only the text (index) on hover
                          hoverlabel = list(namelength = 0))
        
        
        # Add LRC points
        plot <- add_trace(plot,
                          data = subset(data, group == "LRC"),
                          x = ~x, y = ~y, z = ~z,
                          type = "scatter3d", mode = "markers",
                          marker = list(size = 5, color = all_colors[distance_group_indices], symbol = "circle"),
                          # line = list(color = "grey", width = 1),
                          name = "LRC",
                          legendgroup = "LRC",
                          showlegend = T,
                          visible = TRUE,
                          text = data$cord[which(data$group=="LRC")],  # Use the index as the text
                          hoverinfo = "text",    # Show only the text (index) on hover
                          hoverlabel = list(namelength = 0))
        
        # Add DistanceGroup links
        for (pp in 1:nrow(seqFISH_dist_LRC)) {
          pp1 <- which(data$cord == (seqFISH_dist_LRC[pp,1]/1000))
          pp2 <- which(data$cord == (seqFISH_dist_LRC[pp,2]/1000))
          
          plot <- add_trace(plot,
                            x = c(data$x[pp1], data$x[pp2]),
                            y = c(data$y[pp1], data$y[pp2]),
                            z = c(data$z[pp1], data$z[pp2]),
                            type = "scatter3d",
                            mode = "lines",
                            line = list(color = "black", width = 3),
                            name = "LRC",
                            legendgroup = "LRC",
                            showlegend = F,
                            visible = TRUE)
        }
        
        print(plot)
        
        # pie plot
        distances_stat_tmp = seqFISH_dist_LRC
        distances_stat_tmp[,2] = seqFISH_dist_LRC[,1]*1000
        distances_stat_tmp[,1] = seqFISH_dist_LRC[,2]*1000
        
        ### linkage map ####
        rsl = 25
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
        
        circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
          chr=CELL_META$sector.index
          xlim=CELL_META$xlim
          ylim=CELL_META$ylim
          circos.text(mean(xlim), mean(ylim), chr)
        },bg.col="#d3e5ed", bg.border=TRUE, track.height=0.1)
        
        title(paste0('linkage map: ',cell), cex.main = 1.5)
        locs = seq(0,(chr_l+2000000),2000000)
        for (ii in 1:(length(locs)-1)) {
          loc_st = locs[ii]
          loc_ed = locs[ii+1]
          locst = seq(loc_st,loc_ed,rsl*1000)
          tmp1 = which((distances_stat_tmp[,1]>=loc_st) & (distances_stat_tmp[,1]<loc_ed))
          if (length(tmp1)>0) {
            tmp1 = matrix(as.matrix(distances_stat_tmp[tmp1,]),ncol=3)
            for (jj in (ii+1):99) {
              # for (jj in 121:121) {  # create the small links for 240-242 MB locus
              loc_st = locs[jj]
              loc_ed = locs[jj+1]
              loced = seq(loc_st,loc_ed,rsl*1000)
              tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
              if(length(tmp2)!=0){
                tmp2 = matrix(tmp1[tmp2,],ncol=3)
                for (kk in 1:length(locst)) {
                  tmpp = which((tmp2[,1])==locst[kk])
                  if(length(tmpp)!=0){
                    tmpp = matrix(tmp2[tmpp,],ncol=3)
                    for (kkk in 1:length(loced)) {
                      tmppp = which((tmpp[,2])==loced[kkk])
                      if(length(tmppp)!=0){
                        lwdd = 1
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
      }
    }
  }
}

# genomic distances (l) between pair of bins####
# Fig. S6a & b ####
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

# chr21 32.45 MB - 33.35 MB
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

# Fig. S6b  chr21 (all) and chr2 (250 kb) neighbor-locus distance distributions ####

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

IMR_2 <- read.delim("Verify/chromosome2.tsv")
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

# the 3D distances between pairs of neighbor 250k ####
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

# =====================================================================================================
# B. Fig. S6e   micro-C  vs MERFISH
# =====================================================================================================
# INPUTS: micro-c_GSM7661135_IMR90.chr2.output.2k.txt (GEO GSM7661135), 
#         chromosome2.tsv, distances_2-100.txt, distances_chr2L100-3000.txt
# =====================================================================================================

library(geometry)
library(dplyr) 
library(ggplot2)
library(circlize)
options(scipen = 999)  
library(gridExtra)
library(plotly)

####################################################################################################
# Fig. S6e micro-C LRC (linkage map), MERFISH sLRC, overlay and Venn diagram ####
####################################################################################################
## overlay MERFISH and micro-C ####
### input micro-C data and extract >100 Mb & > median of n #####
ff = 'micro-c_GSM7661135_IMR90.chr2.output.2k.txt'
col = '#EEB2D1'
chr_hic = read.delim(paste0(ff), header=T)
print(ff)
print(dim(chr_hic))
chr_hic = data.frame(chr = chr_hic$chrom1,bin_position1 = chr_hic$start1, bin_position2 = chr_hic$start2, count = chr_hic$count)
chr_hic_chr2 = chr_hic[which(chr_hic$chr=='chr2'),]
chr_hic_chr2$dist = chr_hic_chr2$bin_position2-chr_hic_chr2$bin_position1
chr_hic_chr2 = na.omit(chr_hic_chr2)
tmp = which(chr_hic_chr2$dist==0)
chr_hic_chr2 = chr_hic_chr2[-tmp,]
# tmp = which(chr_hic_chr2$dist<=256000)
# md = median((chr_hic_chr2$count[tmp]))
# plot(log10(chr_hic_chr2$dist), log10(chr_hic_chr2$count), pch = '.', col = 'grey',
#      xlim = c(log10(100000-1000),log10(243000000)))
# tmp = which(chr_hic_chr2$count>100)
# points(log10(chr_hic_chr2$dist)[tmp], log10(chr_hic_chr2$count)[tmp], pch = '.',col = 'orange')
# abline(h = 1, col = 'blue',lty = 2)
# abline(v = log10(100000), col = 'blue',lty = 2)
# abline(v = log10(1000000), col = 'blue',lty = 2)

### group 2000 bp to 250 kb
loci_bin = seq(0,244000000,250000)
chr_hic_chr2$posA_bin <- loci_bin[findInterval(chr_hic_chr2$bin_position1, loci_bin)]
chr_hic_chr2$posB_bin <- loci_bin[findInterval(chr_hic_chr2$bin_position2, loci_bin)]

result_df_dplyr <- chr_hic_chr2 %>%
  group_by(posA_bin, posB_bin) %>%
  summarise(sums = sum(count), .groups = "drop")
result_df_dplyr = data.frame(result_df_dplyr)
result_df_dplyr$dist = result_df_dplyr$posB_bin-result_df_dplyr$posA_bin
result_df_dplyr = result_df_dplyr[which(result_df_dplyr$dist!=0),]
print(dim(result_df_dplyr))
print(head(result_df_dplyr))

plot(log10(result_df_dplyr$dist), log10(result_df_dplyr$sums), pch = '.', col = 'grey',
     xlim = c(log10(100000-1000),log10(243000000)))
tmp = which((result_df_dplyr$sums>20)& (result_df_dplyr$dist>100000000))
points(log10(result_df_dplyr$dist)[tmp], log10(result_df_dplyr$sums)[tmp], pch = '.',col = 'orange')
abline(h = 1, col = 'blue',lty = 2)
abline(v = log10(100000), col = 'blue',lty = 2)
abline(v = log10(1000000), col = 'blue',lty = 2)

LRC = result_df_dplyr[which(result_df_dplyr$posB_bin-result_df_dplyr$posA_bin>100000000),]
LRC$dist = LRC$posB_bin-LRC$posA_bin
print(dim(LRC))
print(head(LRC))

microc_cutoff = 20   # strict standard (right panel); use 10 for the loose standard (left panel)
distances_stat2_microc = LRC[which((LRC$sums>microc_cutoff)),]

chr_l = 242000000
rsl = 50
chrom <- (1:(round((chr_l+2000000)/2000000)+1))*2-2
start <- rep(0,(round((chr_l+2000000)/2000000)+1))
end <- rep(2000000/(rsl*1000),(round((chr_l+2000000)/2000000)+1))
chr_sizes_df <- data.frame(chrom,start,end)
# output_file = paste0(ff, '-chr', cc, '_100k.png')
# png(output_file, width = 800, height = 800)
col = '#F3CCE2'#EEB2D1
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
  circos.text(mean(xlim), mean(ylim), chr,
              facing = "clockwise")
},bg.col=col, bg.border=TRUE, track.height=0.1)
title(paste0('micro-c, IMR90'), cex.main = 1)

locs = seq(0,(chr_l+4000000),2000000)
for (ii in 1:(round((chr_l)/2000000))) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,100*1000)
  tmp1 = which((distances_stat2_microc[,1]>=loc_st) & (distances_stat2_microc[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(distances_stat2_microc[tmp1,]),ncol=4)
    for (jj in (ii+1):((round((chr_l+2000000)/2000000)))) {
      # for (jj in 122:122) {# create the small links for 242MB loci
      loc_st = locs[jj]
      loc_ed = locs[jj+1]
      loced = seq(loc_st,loc_ed,100*1000)
      tmp2 = which((tmp1[,2]>=loc_st) & ((tmp1[,2])<loc_ed))
      if(length(tmp2)!=0){
        tmp2 = matrix(tmp1[tmp2,],ncol=4)
        for (kk in 1:(length(locst)-1)) {
          tmpp = which((tmp2[,1]>=locst[kk])&(tmp2[,1]<locst[kk+1]))
          if(length(tmpp)!=0){
            tmpp = matrix(tmp2[tmpp,],ncol=4)
            for (kkk in 1:(length(loced)-1)) {
              tmppp = which((tmpp[,2]>=loced[kkk])&(tmpp[,2]<loced[kkk+1]))
              if(length(tmppp)!=0){
                tmppp = tmpp[tmppp,3]
                if(tmppp!=0){
                  link_thickness <- tmppp/100
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = col, lwd = link_thickness)
                }
              }
            }
          }
        }
      }
    }
  }
}

### input MERFISH ####
IMR_2 <- read.delim("chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))
d3li= 299

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

fils = seq(100,3000,100)
ff = 1
distances_2 = read.table(paste0("distances_2-",fils[ff],".txt"),skip = 1,header = FALSE)
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
distances_stat2_MERFISH = distances_stat2[tmp,]
distances_stat2_MERFISH[,1] = distances_stat2_MERFISH[,1]-1
distances_stat2_MERFISH[,2] = distances_stat2_MERFISH[,2]-1

### plot overlaid ####

# 1. micro-C processing
microc_pairs <- distances_stat2_microc[, 1:2] %>%
  as.data.frame() %>%
  rename(pos1 = posA_bin, pos2 = posB_bin) %>%
  mutate(
    pos1 = as.numeric(pos1),
    pos2 = as.numeric(pos2)
  ) %>%
  distinct()

# 2. MERFISH processing
merfish_pairs <- distances_stat2_MERFISH[, 1:2] %>%
  as.data.frame() %>%
  rename(pos1 = V1, pos2 = V2) %>%
  mutate(
    pos1 = as.numeric(pos1),
    pos2 = as.numeric(pos2)
  ) %>%
  distinct()
common_pairs <- inner_join(microc_pairs, merfish_pairs, by = c("pos1", "pos2"))
microc_only_pairs <- anti_join(microc_pairs, merfish_pairs, by = c("pos1", "pos2"))
merfish_only_pairs <- anti_join(merfish_pairs, microc_pairs, by = c("pos1", "pos2"))

distances_stat_hic2 = cbind(common_pairs,1,1)
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
                tmppp = as.numeric(tmpp[tmppp,4])
                if(tmppp!=0){
                  link_thickness <- tmppp/10
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'bisque4', lwd = 1)
                }
              }
            }
          }
        }
      }
    }  
  }
}

microc_set = paste(microc_pairs[,1], microc_pairs[,2])
merfish_set = paste(merfish_pairs[,1], merfish_pairs[,2])
venn_data <- list(
  MicroC = microc_set,
  MERFISH = merfish_set
)
library(VennDiagram)
venn_plot <- draw.pairwise.venn(
  area1 = length(microc_set),          
  area2 = length(merfish_set),         
  cross.area = length(intersect(microc_set, merfish_set)),  
  category = c("MicroC", "MERFISH"),   
  fill = c("lightgreen", "lightblue"),
  alpha = 0.6,                         
  lwd = 2,                              
  lty = "solid",                      
  fontface = "bold",                   
  fontfamily = "sans",                 
  cat.cex = 1.2,                      
  cex = 1.1,                          
  cat.pos = c(30, 330),               
  cat.dist = 0.07,                   
  main = "MicroC vs MERFISH",         
  # print.mode = 'none'
  main.cex = 1.5                      
)

# Function to calculate point-to-triangle distance
point_to_triangle_distance <- function(point, triangle_vertices) {
  # Calculate normal vector of the triangle
  v1 <- triangle_vertices[2,] - triangle_vertices[1,]
  v2 <- triangle_vertices[3,] - triangle_vertices[1,]
  normal <- crossprod(matrix(c(v1, v2), nrow=2))
  normal <- normal / sqrt(sum(normal^2))
  
  # Calculate distance from point to plane of triangle
  distance <- abs(sum((point - triangle_vertices[1,]) * normal))
  return(distance)
}
# Function to calculate Euclidean distance between two 3D points
calc_dist <- function(p1, p2) {
  if (any(is.na(p1)) || any(is.na(p2))) {
    return(NA)
  }
  sqrt(sum((p1 - p2)^2))
}
# Function to create distance matrix from xyz coordinates
create_dist_matrix <- function(coords) {
  # Check if input is matrix
  if (!is.matrix(coords)) {
    stop("Input must be a matrix")
  }
  
  # Check if matrix has 3 columns
  if (ncol(coords) != 3) {
    stop("Matrix must have exactly 3 columns (x, y, z)")
  }
  
  n <- nrow(coords)
  dist_matrix <- matrix(NA, nrow = n, ncol = n)
  
  # Calculate pairwise distances
  for (i in 1:n) {
    for (j in i:n) {
      dist <- calc_dist(coords[i,], coords[j,])
      dist_matrix[i,j] <- dist
      dist_matrix[j,i] <- dist  # Make matrix symmetric
    }
  }
  
  return(dist_matrix)
}
#Extract distances where index difference > 16 ~ 48 Mb
extract_distant_pairs <- function(dist_matrix) {
  n <- nrow(dist_matrix)
  
  # Initialize vectors to store results
  distances <- c()
  
  # Extract distances where index difference > 16 ~ 48 Mb
  for (i in 1:(n-16)) {
    for (j in (i+16):n) {  
      distances <- c(distances, dist_matrix[i,j])
    }
  }
  
  # Combine into a data frame
  result <- distances
  return(result)
}

####################################################################################################
# Fig. S6f #####
####################################################################################################
## 3M--Distances among nuclear structure before and after binding 
## calculate distances 
# input original file
M3strut <- read.delim("genomic-scale-with transcription and nuclear bodies.tsv")
print(colnames(M3strut))
print(dim(M3strut))

# extract chr info
chr = rep(0,length(M3strut$genomic.coordinate))
for (ii in 1:nrow(M3strut)) {
  tmp = strsplit(M3strut$genomic.coordinate[ii],split = ':')[[1]]
  chr[ii] = tmp[1]
}
# extract tag loci info
coordinate_tmp = unique(M3strut$genomic.coordinate)
print(length(coordinate_tmp))
coordinates = matrix(0,nrow = nrow(M3strut), ncol = 2)
for (ii in 1:length(coordinate_tmp)) {
  co_tmp = which(M3strut$genomic.coordinate==coordinate_tmp[ii])
  tmp = strsplit(coordinate_tmp[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  coordinates[co_tmp,1] = as.numeric(tmp[1])
  coordinates[co_tmp,2] = as.numeric(tmp[2])
}
coordinates = coordinates[,1]

# check the length of each chr
chr_idx = unique(chr)
for (ii in chr_idx) {
  print(c(ii,min(as.numeric(as.character(coordinates[which(chr ==ii)]))),max(as.numeric(as.character(coordinates[which(chr ==ii)]))),length(which(chr ==ii))))
}

# Calculate the neighbor distance and distances > 50 Mb for chr 1-7
## start calculation ####
neighbor3M_list <- vector("list", 7)
names(neighbor3M_list) <- paste0("chr", 1:7)
LongDist_3M_list <- vector("list", 7)
names(LongDist_3M_list) <- paste0("chr", 1:7)

# Loop through 7 chromosomes
for (chrr in 1:7) {
  cat("Processing chromosome", chrr, "\n")
  
  # Get chromosome-specific data
  chr_n = which(chr ==paste0('chr',chrr))
  M3strut_chr = M3strut[chr_n,]
  coordinate_chr = seq(0,as.character(coordinates[chr_n])%>%as.numeric()%>%max(),3000000)+2950000
  
  print(colnames(M3strut_chr))
  print(dim(M3strut_chr))
  print(length(coordinate_chr))
  
  cells_chr = unique(M3strut_chr$cell.number)
  idx1 = 0
  idx2 = 0
  for (i in 1:(length(coordinate_chr)-16)) {
    for (j in (i+16):length(coordinate_chr)) {
      idx1 <- c(idx1, coordinate_chr[i])
      idx2 <- c(idx2, coordinate_chr[j])
    }
  }
  idx_chr = cbind(idx1,idx2)
  idx_chr = idx_chr[-1,]/1000000+0.05
  print(dim(idx_chr))
  
  # Initialize matrices for this chromosome
  # Note: Different chromosomes will have different numbers of regions (rows)
  neighbor3M_chr <- matrix(0, nrow = length(coordinate_chr)-1, ncol = length(cells_chr) * 2)
  LongDist_3M_chr <- matrix(0, nrow = nrow(idx_chr), ncol = length(cells_chr) * 2)
  
  # Process each cell for this chromosome
  for (ss in 1:length(cells_chr)) {
    if (ss %% 500 == 0) {
      cat("  Cell", ss, "of", length(cells_chr), "\n")
    }
    
    # Process homolog 1
    cells_tmp <- which(M3strut_chr$cell.number == cells_chr[ss] &
                         M3strut_chr$homolog.number == 1)
    M3strut_tmp <- M3strut_chr[cells_tmp, ]
    coordinates_tmp = rep(0,nrow(M3strut_tmp))
    for (ii in 1:length(coordinates_tmp)) {
      tmp = strsplit(M3strut_tmp$genomic.coordinate[ii],split =':')[[1]]
      tmp = strsplit(tmp[2],split ='-')[[1]]
      coordinates_tmp[ii] = as.numeric(tmp[1])
    }
    # print(length(coordinates_tmp))
    tmp = match(coordinate_chr,coordinates_tmp)
    xx <- as.matrix(M3strut_tmp[tmp, 1:3])%>%create_dist_matrix()
    neighbor3M_chr[, (ss * 2 - 1)] <- xx[row(xx) == col(xx) - 1]
    LongDist_3M_chr[, (ss * 2 - 1)] <- extract_distant_pairs(xx)
    
    # Process homolog 2
    cells_tmp <- which(M3strut_chr$cell.number == cells_chr[ss] &
                         M3strut_chr$homolog.number == 2)
    M3strut_tmp <- M3strut_chr[cells_tmp, ]
    coordinates_tmp = rep(0,nrow(M3strut_tmp))
    for (ii in 1:length(coordinates_tmp)) {
      tmp = strsplit(M3strut_tmp$genomic.coordinate[ii],split =':')[[1]]
      tmp = strsplit(tmp[2],split ='-')[[1]]
      coordinates_tmp[ii] = as.numeric(tmp[1])
    }
    # print(length(coordinates_tmp))
    tmp = match(coordinate_chr,coordinates_tmp)
    xx <- as.matrix(M3strut_tmp[tmp, 1:3])%>%create_dist_matrix()
    neighbor3M_chr[, (ss * 2)] <- xx[row(xx) == col(xx) - 1]
    LongDist_3M_chr[, (ss * 2)] <- extract_distant_pairs(xx)
  }
  
  # Add idx to LongDist_3M_chr
  LongDist_3M_chr <- cbind(idx_chr, LongDist_3M_chr)
  
  # Store results for this chromosome with chromosome label
  neighbor3M_list[[paste0("chr", chrr)]] <- neighbor3M_chr
  LongDist_3M_list[[paste0("chr", chrr)]] <- LongDist_3M_chr
  
  cat("Completed chromosome", chrr, "- neighbor3M dims:", dim(neighbor3M_chr),
      "LongDist_3M dims:", dim(LongDist_3M_chr), "\n")
}

# The final results are stored as lists, not combined matrices
neighbor3M <- neighbor3M_list  # List of 7 matrices
LongDist_3M <- LongDist_3M_list  # List of 7 matrices

chr_length = c(249000000,242000000,198000000,190000000,182000000,171000000,159000000)
d3li = 920
d3li = 299
for(chrr in 1:7) {
  distdiff = LongDist_3M[[chrr]][,2]-LongDist_3M[[chrr]][,1]
  tmp = which(distdiff>100)
  LongDist_3M1 = LongDist_3M[[chrr]][tmp,c(1:2)]
  LongDist_3M2 = LongDist_3M[[chrr]][tmp,-c(1:2)]
  print(dim(LongDist_3M2))
  
  LongDist_3M1 <- LongDist_3M1[!apply(is.na(LongDist_3M2), 1, all), ]
  LongDist_3M2 <- LongDist_3M2[!apply(is.na(LongDist_3M2), 1, all), ]
  print(dim(LongDist_3M1))
  print(dim(LongDist_3M2))
  
  LongDist_3M2[which(LongDist_3M2>d3li,arr.ind = T)] = 0
  LongDist_3M2[which(LongDist_3M2!=0, arr.ind = T)] = 1
  data = rowSums(LongDist_3M2,na.rm = T)/ncol(LongDist_3M2)
  
  par(mar = c(2,2,1,1))
  a = boxplot(data,notch = T,col = 'orange',outpch = 19,outcex = 0.5,outcol = 'darkorange3',horizontal = T,plot = T)
  hist(data, breaks = 100, prob = TRUE, col = "orange", xlab = "Values", main = paste('chr',chrr),add = T)
  # box(lwd = 3)
  # abline(v = a$stats[5],lty = 3, col = 'orchid', lwd = 2)
  stableline = a$stats[5]
  # stableline = nrow(LongDist_3M2)*0.4
  # stableline =46
  tmp = which(data>stableline)
  print(length(tmp))
  distances_stat22 = cbind(LongDist_3M1[tmp,]*1000000,data[tmp])
  ### linkage map of sLRC in 3M ####
  chr2 = chr_length[chrr]
  rsl = 1000
  chrom <- (1:(round((chr2+2000000)/2000000)+1))*2-2
  start <- rep(1,(round((chr2+2000000)/2000000)+1))
  end <- rep(3,(round((chr2+2000000)/2000000)+1))
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
    circos.text(mean(xlim), mean(ylim), chr)
  },bg.col=transparent_color, bg.border=TRUE, track.height=0.1)
  # circos.rect(
  #   sector.index = 244,
  #   # track.index	= 1,
  #   xleft = 1,
  #   xright = 2,
  #   col = 'white',
  #   border = NA,
  #   ytop = 1.5,
  #   ybottom=-0.5
  # )
  title(paste('sLRC 250k->3M chr:',chrr), cex.main = 1.5)
  #plot DNA FISH 
  locs = seq(0,(chr2+2000000),2000000)
  for (ii in 1:(length(locs)-1)) {
    loc_st = locs[ii]
    loc_ed = locs[ii+1]
    locst = seq(loc_st,loc_ed,rsl*1000)
    tmp1 = which((distances_stat22[,1]>=loc_st) & (distances_stat22[,1]<loc_ed))
    if (length(tmp1)>0) {
      tmp1 = matrix(as.matrix(distances_stat22[tmp1,]),ncol=3)
      for (jj in (ii+1):((round((chr2)/2000000)))) {
        # for (jj in 122:122) {  # create the small links for 240-242 MB locus
        loc_st = locs[jj]
        loc_ed = locs[jj+1]
        loced = seq(loc_st,loc_ed,rsl*1000)
        tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<loc_ed))
        if(length(tmp2)!=0){
          tmp2 = matrix(tmp1[tmp2,],ncol=3)
          for (kk in 1:length(locst)) {
            tmpp = which((tmp2[,1])==locst[kk])
            if(length(tmpp)!=0){
              tmpp = matrix(tmp2[tmpp,],ncol=3)
              for (kkk in 1:length(loced)) {
                tmppp = which((tmpp[,2])==loced[kkk])
                if(length(tmppp)!=0){
                  lwdd = tmpp[tmppp,3]/50
                  # lwdd = 1
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
}

####################################################################################################
# Fig. S6g, h,i####
####################################################################################################
## DNA seqFISH data clean up####
probes = read.csv('41586_2025_8838_MOESM5_ESM.csv')
probes_ID = strsplit(probes$name, split = '-')
probes_ID = do.call(rbind,probes_ID)
probes = cbind(probes_ID,probes$name,probes$chrom,probes$start)
colnames(probes) = c('chr','probe_ID','names','chrom','start')
probes = data.frame(probes)

file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)   # all FOVs pos0-pos9
for (file_name in file_names) {
  fls = paste0('seqFISH/NMuMG/',file_name,'.csv')
  print(fls)
  seqFISH_DNA <- read.csv(fls)
  print(dim(seqFISH_DNA))
  # remove the duplicated dots (all fov, cellID, loci, location and allele are the same)
  seqFISH_DNA_tmp = cbind(seqFISH_DNA$fov, seqFISH_DNA$cellID, seqFISH_DNA$name, seqFISH_DNA$dbscan_ldp_nbr_allele)
  duplicated_all <- duplicated(seqFISH_DNA_tmp) | duplicated(seqFISH_DNA_tmp, fromLast = TRUE)
  print(length(which(duplicated_all)))
  seqFISH_DNA = seqFISH_DNA[-which(duplicated_all),]
  # This returns a list where each element is a vector of the split parts
  split_list <- strsplit(seqFISH_DNA$name, "-")
  loci <- do.call(rbind, split_list)
  loci[, 2] <- as.numeric(loci[, 2])
  seqFISH_DNA2 = data.frame(
    cbind(
      seqFISH_DNA$cellID,seqFISH_DNA$chrom,loci[,2],
      seqFISH_DNA$x_um, seqFISH_DNA$y_um,  seqFISH_DNA$z_um,
      seqFISH_DNA$dbscan_ldp_nbr_allele
    )
  )
  colnames(seqFISH_DNA2) = c('cell', 'chr', 'loci','x_um','y_um','z_um', 'allele')
  for (cc in paste0('chr',1:1)) {
    seqFISH_DNA2_tmp3 = seqFISH_DNA2[which(seqFISH_DNA2$chr==cc),]
    print(cc)
    print(dim(seqFISH_DNA2_tmp3))
    # alleles = sort(unique(seqFISH_DNA2_tmp3$allele))
    # alleles = alleles[which(alleles!='-1')]
    for (allel in c(1,2)) {
      seqFISH_DNA2_tmp = seqFISH_DNA2_tmp3[which(seqFISH_DNA2_tmp3$allele==allel),]
      print(allel)
      print(dim(seqFISH_DNA2_tmp))
      
      cells = seqFISH_DNA2_tmp$cell%>%unique()%>%as.character()%>%as.numeric()%>%sort()
      spots2 = unique(seqFISH_DNA2_tmp$loci)%>%as.character()%>%as.numeric()%>%sort()
      print(c(length(cells),length(spots2)))
      if (length(cells)==1) {
        next
      }
      for (cell in cells[1:1]) {
        tmp = which(seqFISH_DNA2_tmp$cell==cell)
        seqFISH_DNA2_tmp2 = seqFISH_DNA2_tmp[tmp,]
        # print(dim(seqFISH_DNA2_tmp2))
        coords <- cbind(as.numeric(seqFISH_DNA2_tmp2$x_um),as.numeric(seqFISH_DNA2_tmp2$y_um),as.numeric(seqFISH_DNA2_tmp2$z_um))
        tmp2 = match(spots2,seqFISH_DNA2_tmp2$loci%>%as.character()%>%as.numeric())
        coords2 =  matrix(NA, ncol = 3, nrow = length(tmp2))
        coords2[which(!is.na(tmp2)),] = coords[tmp2[which(!is.na(tmp2))],]
        
        dist_matrix <- as.matrix(dist(coords2))
        
        rownames(dist_matrix) <- spots2
        colnames(dist_matrix) <- spots2
        combinations <- which(lower.tri(dist_matrix, diag = FALSE), arr.ind = TRUE)
        
        # Create the three-column result
        result <- data.frame(
          ID1 = spots2[combinations[, "row"]],
          ID2 = spots2[combinations[, "col"]],
          Distance = dist_matrix[lower.tri(dist_matrix)]
        )
        
        # Sort the result (optional)
        chr_dist <- result[order(result$ID1, result$ID2), ]
        chr_dist$Distance = round(as.numeric(chr_dist$Distance*1000),2)
      }
      for (cell in cells[2:length(cells)]) {
        tmp = which(seqFISH_DNA2_tmp$cell==cell)
        seqFISH_DNA2_tmp2 = seqFISH_DNA2_tmp[tmp,]
        # print(dim(seqFISH_DNA2_tmp2))
        coords <- cbind(as.numeric(seqFISH_DNA2_tmp2$x_um),as.numeric(seqFISH_DNA2_tmp2$y_um),as.numeric(seqFISH_DNA2_tmp2$z_um))
        tmp2 = match(spots2,seqFISH_DNA2_tmp2$loci%>%as.character()%>%as.numeric())
        coords2 =  matrix(NA, ncol = 3, nrow = length(tmp2))
        coords2[which(!is.na(tmp2)),] = coords[tmp2[which(!is.na(tmp2))],]
        
        dist_matrix <- as.matrix(dist(coords2))
        
        rownames(dist_matrix) <- spots2
        colnames(dist_matrix) <- spots2
        combinations <- which(lower.tri(dist_matrix, diag = FALSE), arr.ind = TRUE)
        
        # Create the three-column result
        result <- data.frame(
          ID1 = spots2[combinations[, "row"]],
          ID2 = spots2[combinations[, "col"]],
          Distance = dist_matrix[lower.tri(dist_matrix)]
        )
        
        # Sort the result (optional)
        result <- result[order(result$ID1, result$ID2), ]
        chr_dist = cbind(chr_dist,round(as.numeric(result$Distance%>%as.character())*1000,2))
      }
      colnames(chr_dist) = c('loc1', 'loc2',paste0('cell-',cells))
      write.csv(chr_dist, file = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele',allel, ".csv"), row.names = FALSE)
    }
  }
}

file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
for (cc in paste0('chr',c(1))) {
  for (file_name in file_names) {   # original: file_names[-1] (pos0 skipped; its split files are read later, so they must come from an earlier run)
    for (allel in c(1,2)) {
      fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele',allel, ".csv")
      print(fls)
      seqFISH_DNA <- read.csv(fls)
      loci = cbind(paste0(cc,'-',seqFISH_DNA$loc1),paste0(cc,'-',seqFISH_DNA$loc2))
      
      tmp = match(loci[,1],probes$name)
      loci_st = probes$start[tmp]
      tmp= match(loci[,2],probes$name)
      loci_ed = probes$start[tmp]
      
      seqFISH_DNA$loc1 = as.numeric(as.character(loci_st))/1000
      seqFISH_DNA$loc2 = as.numeric(as.character(loci_ed))/1000
      
      genomic_dist = abs(as.numeric(as.character(loci_st))-as.numeric(as.character(loci_ed)))
      tmp =  which(genomic_dist==25000)
      print(length(tmp))
      write.csv(seqFISH_DNA[tmp,], file = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele', allel,"_neighbor.csv"), row.names = FALSE)
      tmp =  which((genomic_dist>25000)&(genomic_dist<=50000000))
      print(length(tmp))
      if (length(tmp)>0) {
        tmp = sample(tmp,length(tmp)/2)
        write.csv(seqFISH_DNA[tmp,], file = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele', allel,"_near.csv"), row.names = FALSE)
      }
      tmp =  which((genomic_dist>50000000)&(genomic_dist<=100000000))
      print(length(tmp))
      if (length(tmp)>0) {
        write.csv(seqFISH_DNA[tmp,], file = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele', allel,"_50Mb.csv"), row.names = FALSE)
      }
      tmp =  which(genomic_dist>100000000)
      print(length(tmp))
      if (length(tmp)>0) {
        write.csv(seqFISH_DNA[tmp,], file = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc,'-allele', allel, "_100Mb.csv"), row.names = FALSE)
      }
    }
  }
}

cc = 'chr1'
## basci dist. info ####
### neighbors 3D distance ####
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
for (file_name in file_names[1]) {
  for (allel in c(1)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_neighbor.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    # print(seqFISH_DNA[1:5,1:5])
    loci_pair = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = na.omit(unlist(seqFISH_DNA[,3:ncol(seqFISH_DNA)]))
    LRC = which(seqFISH_DNA[,3:ncol(seqFISH_DNA)]<299,arr.ind = T)
    loci_pair_LRC = loci_pair[LRC[,1],]
    tab <- table(c(loci_pair_LRC[,1], loci_pair_LRC[,2]))
    print(head(sort(tab, decreasing = TRUE), 10))
    seqFISH_DNA_nbr = seqFISH_DNA2
    print(length(seqFISH_DNA_nbr))
  }
}
for (file_name in file_names[1]) {
  for (allel in c(2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_neighbor.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    # print(seqFISH_DNA[1:5,1:5])
    loci_pair = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = na.omit(unlist(seqFISH_DNA[,3:ncol(seqFISH_DNA)]))
    LRC = which(seqFISH_DNA[,3:ncol(seqFISH_DNA)]<299,arr.ind = T)
    loci_pair_LRC = loci_pair[LRC[,1],]
    tab <- table(c(loci_pair_LRC[,1], loci_pair_LRC[,2]))
    print(head(sort(tab, decreasing = TRUE), 10))
    seqFISH_DNA_nbr = seqFISH_DNA2
    print(length(seqFISH_DNA_nbr))
  }
}
for (file_name in file_names[-c(1,8)]) {
  for (allel in c(1,2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_neighbor.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    # print(seqFISH_DNA[1:5,1:5])
    loci_pair = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = na.omit(unlist(seqFISH_DNA[,3:ncol(seqFISH_DNA)]))
    LRC = which(seqFISH_DNA[,3:ncol(seqFISH_DNA)]<299,arr.ind = T)
    loci_pair_LRC = loci_pair[LRC[,1],]
    tab <- table(c(loci_pair_LRC[,1], loci_pair_LRC[,2]))
    print(head(sort(tab, decreasing = TRUE), 10))
    seqFISH_DNA_nbr = c(seqFISH_DNA_nbr,seqFISH_DNA2)
    print(length(seqFISH_DNA_nbr))
  }
}

### 250 kb 3D distance ####
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
for (file_name in file_names[1]) {
  for (allel in c(1)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_near.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    
    loci_pair_tmp = seqFISH_DNA[,1:2]
    loci_pair_tmp = which(abs(loci_pair_tmp[,1]-loci_pair_tmp[,2])==250)
    print(length(loci_pair_tmp))
    seqFISH_DNA_250k = na.omit(unlist(seqFISH_DNA[loci_pair_tmp,3:ncol(seqFISH_DNA)]))
    print(length(seqFISH_DNA_250k))
  }
}
for (file_name in file_names[1]) {
  for (allel in c(2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_near.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    
    loci_pair_tmp = seqFISH_DNA[,1:2]
    loci_pair_tmp = which(abs(loci_pair_tmp[,1]-loci_pair_tmp[,2])==250)
    print(length(loci_pair_tmp))
    
    seqFISH_DNA2 = na.omit(unlist(seqFISH_DNA[loci_pair_tmp,3:ncol(seqFISH_DNA)]))
    seqFISH_DNA_250k = c(seqFISH_DNA_250k,seqFISH_DNA2)
    print(length(seqFISH_DNA_250k))
  }
}
for (file_name in file_names[-c(1,8)]) {
  for (allel in c(1, 2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-',cc, '-allele', allel, "_near.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    
    loci_pair_tmp = seqFISH_DNA[,1:2]
    loci_pair_tmp = which((loci_pair_tmp[,1]-loci_pair_tmp[,2])==250)
    print(length(loci_pair_tmp))
    seqFISH_DNA_250k = c(seqFISH_DNA_250k,na.omit(unlist(seqFISH_DNA[loci_pair_tmp,3:ncol(seqFISH_DNA)])))
    print(length(seqFISH_DNA_250k))
  }
}

A = as.vector(seqFISH_DNA_250k)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
print(quartiles)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 300, 
     main = fls, 
     xlim = c(0,5000),
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2) 
abline(v = quartiles[1], col ="orange", lty = 2, lwd = 2)  
abline(v = 299, col ="darkgreen", lty = 2, lwd = 2)
legend("topright", 
       legend = c( "299 nm", "1st QTL", "Median"),
       col = c("darkgreen", "orange", "brown3"),
       lty = c( 2), lwd = 2)

### 100 Mb 3D distance ####
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
for (file_name in file_names[1]) {
  for (allel in c(1)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    locs = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = seqFISH_DNA[,3:ncol(seqFISH_DNA)]
    seqFISH_DNA2 = unlist(seqFISH_DNA2)
    seqFISH_DNA2 = na.omit(seqFISH_DNA2)
    
    seqFISH_DNA_100 = seqFISH_DNA2
    print(length(seqFISH_DNA_100))
  }
}
for (file_name in file_names[1]) {
  for (allel in c(2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    locs = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = seqFISH_DNA[,3:ncol(seqFISH_DNA)]
    seqFISH_DNA2 = unlist(seqFISH_DNA2)
    seqFISH_DNA2 = na.omit(seqFISH_DNA2)
    
    seqFISH_DNA_100 = c(seqFISH_DNA_100, seqFISH_DNA2)
    print(length(seqFISH_DNA_100))
  }
}
for (file_name in file_names[-c(1,8)]) {
  for (allel in c(1, 2)){
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    locs = seqFISH_DNA[,1:2]
    seqFISH_DNA2 = seqFISH_DNA[,3:ncol(seqFISH_DNA)]
    seqFISH_DNA2 = unlist(seqFISH_DNA2)
    seqFISH_DNA2 = na.omit(seqFISH_DNA2)
    
    seqFISH_DNA_100 = c(seqFISH_DNA_100, seqFISH_DNA2)
    print(length(seqFISH_DNA_100))
  }
}

color_A <- "orange"  # Soft blue
color_B <- "#9b59b6"  # Coral red
color_C <- "darkolivegreen3"  # Coral red
A = as.numeric(seqFISH_DNA_100)
# calculate l between 100 MB bins
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
hist(A, probability = TRUE, breaks = 200, xlim = c(1,8000), 
     main = "3D dist. between pair of loci >100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_A, border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

B = as.numeric(seqFISH_DNA_250k)
# calculate l between 50 MB bins
quartiles <- quantile(B, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(B, probability = TRUE, breaks = 200, xlim = c(1,8000), 
     main = "3D dist. between pair of loci 50-100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_B, border = "white")
box(lwd = 3)
lines(density(B,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

C = as.vector(seqFISH_DNA_nbr)
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
abline(v = quantile(B,c(0.25,0.5,0.75),na.rm=T)[1], col = "darkblue", lty = 2, lwd = 2) 
abline(v = quantile(C,c(0.25,0.5,0.75),na.rm=T)[2], col = "darkgreen", lty = 2, lwd = 2) 
abline(v = 299, col = "red", lty = 2, lwd = 2)
# Add a legend
legend("topright", legend = c("bin pairs >100M","bin pairs = 250k", "bin pairs = 25K","1st QTL of 250 kb pairs","median of 25 kb pairs","299 nm"),
       col = c(adjustcolor(color_A, alpha.f = 0.5), adjustcolor(color_B, alpha.f = 0.5), adjustcolor(color_C, alpha.f = 0.5),adjustcolor("darkblue", alpha.f = 0.5),adjustcolor("darkgreen", alpha.f = 0.5),adjustcolor("red", alpha.f = 0.5)),
       lty = c(1,1,1,2), lwd = 2,
       border = NA, bty = "n")

####################################################################################################
# Fig. S6h  ####
####################################################################################################
### LRC grouped to 250 kb resolution ####
loci_bin = seq(0,244000000,250000)

sLRC = matrix(0,nrow = 1, ncol = 4)
colnames(sLRC) = c('loc1_bin','loc2_bin','LRC','not_na')
cell_number = 0
loci_number = 0
for (file_name in file_names[-c(8)]) {
  for (allel in c(1,2)) {
    fls = paste0("seqFISH/NMuMG/dist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    print(dim(seqFISH_DNA))
    cell_number = cell_number+ncol(seqFISH_DNA)-2
    loci_number = c(loci_number, as.numeric(table(which(!is.na(seqFISH_DNA[,-c(1:2)]),arr.ind = T)[,2])))
    
    locs = seqFISH_DNA[,1:2]
    locs[,1] = seqFISH_DNA[,2]
    locs[,2] = seqFISH_DNA[,1]
    seqFISH_DNA2 = seqFISH_DNA[,3:ncol(seqFISH_DNA)]
    tmp = which(seqFISH_DNA2<299,arr.ind = T)
    seqFISH_DNA2_LRC = matrix(0,nrow = nrow(seqFISH_DNA2), ncol = ncol(seqFISH_DNA2))
    seqFISH_DNA2_LRC[tmp] = 1
    seqFISH_DNA2_LRC = cbind(locs,seqFISH_DNA2_LRC)
    
    seqFISH_DNA2_LRC$posA_bin <- loci_bin[findInterval(seqFISH_DNA2_LRC$loc1*1000, loci_bin)]
    seqFISH_DNA2_LRC$posB_bin <- loci_bin[findInterval(seqFISH_DNA2_LRC$loc2*1000, loci_bin)]
    result_df_dplyr <- seqFISH_DNA2_LRC %>%
      group_by(posA_bin, posB_bin) %>%  # 分组列正常使用
      summarise(
        across(
          .cols = -c(loc1, loc2),       # 只排除loc1/loc2，其余列（包括cell.*）求和
          .fns = ~sum(.x, na.rm = TRUE),# 忽略NA，避免求和失败
          .names = "{.col}_sum"         # 求和列名：原列名+_sum
        ),
        .groups = "drop"  
      ) %>%
      as.data.frame()
    seqFISH_DNA2_LRC = result_df_dplyr[,-c(1:2)]
    
    tmp = which(!is.na(seqFISH_DNA2),arr.ind = T)
    seqFISH_DNA2_na = matrix(0,nrow = nrow(seqFISH_DNA2), ncol = ncol(seqFISH_DNA2))
    seqFISH_DNA2_na[tmp] = 1
    seqFISH_DNA2_na = cbind(locs,seqFISH_DNA2_na)
    
    seqFISH_DNA2_na$posA_bin <- loci_bin[findInterval(seqFISH_DNA2_na$loc1*1000, loci_bin)]
    seqFISH_DNA2_na$posB_bin <- loci_bin[findInterval(seqFISH_DNA2_na$loc2*1000, loci_bin)]
    result_df_dplyr <- seqFISH_DNA2_na %>%
      group_by(posA_bin, posB_bin) %>%  # 分组列正常使用
      summarise(
        across(
          .cols = -c(loc1, loc2),       # 只排除loc1/loc2，其余列（包括cell.*）求和
          .fns = ~sum(.x, na.rm = TRUE),# 忽略NA，避免求和失败
          .names = "{.col}_sum"         # 求和列名：原列名+_sum
        ),
        .groups = "drop"  
      ) %>%
      as.data.frame()
    seqFISH_DNA2_na = result_df_dplyr[,-c(1:2)]
    seqFISH_DNA2_na[which(seqFISH_DNA2_na!=0,arr.ind = T)]=1
    
    
    sLRC_tmp = cbind(result_df_dplyr$posA_bin, result_df_dplyr$posB_bin,rowSums(seqFISH_DNA2_LRC),rowSums(seqFISH_DNA2_na))
    colnames(sLRC_tmp) = c('loc1_bin','loc2_bin','LRC','not_na')
    sLRC = rbind(sLRC,sLRC_tmp)
  }
}
print(cell_number)
sLRC = sLRC[-1,]
sLRC2 = as.data.frame(sLRC, stringsAsFactors = FALSE)
result <- sLRC2 %>%
  group_by(loc1_bin, loc2_bin) %>%
  summarise(
    LRC = sum(LRC),  # aa列求和
    not_na = sum(not_na),  # bb列求和
    .groups = "drop"
  )
result = data.frame(result)
hist(result$LRC, probability = TRUE, 
     main = 'exsiting of LRC',
     xlab = "# cells LRC exsiting", col = "lightblue", border = "white")
hist(result$not_na, probability = TRUE, 
     main = '# detected cells',
     xlab = "# detected cells", col = "lightblue", border = "white")
plot(result$LRC, result$not_na,pch = '.', col = 'blue',cex = 3)
hist(result$not_na)

####################################################################################################
# Fig. S6i ####
####################################################################################################
tmp = which(result$LRC>10)
paired_loci2 = result[tmp,]
print(c(nrow(result),length(tmp)))

chr_l = 196000000 
rsl = 250
chrom <- (1:(round((chr_l+2000000)/2000000)+1))*2-2
start <- rep(1,(round((chr_l+2000000)/2000000)+1))
end <- rep(2000000/(rsl*1000),(round((chr_l+2000000)/2000000)+1))
chr_sizes_df <- data.frame(chrom,start,end)

col = 'lightblue'
circos.clear()
circos.par(cell.padding=c(0,0,0,0),
           track.margin=c(0,0.05),
           start.degree = 90,
           gap.degree = 0.0,
           clock.wise = TRUE)
circos.initialize(factors = chr_sizes_df$chrom,
                  xlim = cbind(chr_sizes_df$start, chr_sizes_df$end))
transparent_color <- adjustcolor(col, alpha.f = 0.4)
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr,cex = 0.5,facing = "clockwise")
},bg.col=transparent_color, bg.border=TRUE, track.height=0.1)
# circos.rect(
#   sector.index = 244,
#   # track.index	= 1,
#   xleft = 1,
#   xright = 8,
#   col = 'white',
#   border = NA,
#   ytop = 1.5,
#   ybottom=-0.5
# )
# circos.rect(
#   sector.index = 242,
#   # track.index	= 1,
#   xleft = 3,
#   xright = 8,
#   col = 'white',
#   border = NA,
#   ytop = 1.5,
#   ybottom=-0.5
# )

locs = seq(0,(chr_l+4000000),2000000)
for (ii in 1:(length(locs)-1)) {
  loc_st = locs[ii]
  loc_ed = locs[ii+1]
  locst = seq(loc_st,loc_ed,rsl*1000)
  tmp1 = which(((paired_loci2[,1])>=loc_st) & (paired_loci2[,1]<loc_ed))
  if (length(tmp1)>0) {
    tmp1 = matrix(as.matrix(paired_loci2[tmp1,]),ncol=4)
    for (jj in (ii+1):(length(locs)-1)) {
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
                lwdd = (tmpp[tmppp,3]/20)
                if(length(tmppp)!=0){
                  circos.link(ii*2-2, kk, jj*2-2, kkk, col = col,lwd = lwdd)
                }
              }
            }
          }
        }
      }
    }
  }
}
