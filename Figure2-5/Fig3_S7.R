# =====================================================================================================
# Fig. 3a, 3e, 3f, 3g, 3h, 3i   and   Fig. S7a, S7c, S7d, S7e, S7f   --  sLRC -> metaLRC -> multivalent binding model
# 
# INPUT: chromosome2.tsv, distances_chr2L100-3000.txt (same as Fig2_S6.R).
# =====================================================================================================

####################################################################################################
# SETUP (libraries, phi correlation function, sLRC definition)
####################################################################################################
library(ggplot2)
library(ggbreak)
library(readxl)
library(VennDiagram)
library("DESeq2")
library(factoextra)
library(rtracklayer)
library(hicrep)
library(heatmaply)
library(pheatmap)
library(plotly)
library("viridis")
library(ggExtra)
library(circlize)
library(igraph)


phi <- function (x) {
  ncx = nrow(x)
  ncy = nrow(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      x2 <- x[i,]
      y2 <- x[j,]
      aa = length(which(x2==1&y2==1))/ncol(x)
      dd = length(which(x2==0&y2==0))/ncol(x)
      bb = length(which(x2==1&y2==0))/ncol(x)
      cc = length(which(x2==0&y2==1))/ncol(x)
      phi = (aa*dd-bb*cc)/sqrt((aa+bb)*(cc+dd)*(aa+cc)*(bb+dd))
      r[i, j] <-   phi
    }
  }
  r <- r + t(r) - diag(diag(r))
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
} 

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

spots2 = unique(Coordinate)

st = 0
ed = 0
for (jj in 1:length((spots2))) {
  st = c(st,rep(spots2[jj],length(spots2)-jj))
  ed = c(ed, spots2[-c(1:jj)])
}
st = st[-1]
ed = ed[-1]
lner = as.numeric(ed)-as.numeric(st)
tt1 = which(lner>100000000)

distances_2m100 = read.table(paste0('distances_chr2L100-3000.txt'),skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

chr2 = 242200000
rsl = 50
d3li = 299

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))
colnames(distances_stat2) = c('start','end','genom.dist.','numblks')
distances_stat2 = data.frame(distances_stat2)
print(dim(distances_stat2))
print(head(distances_stat2))
print(head(IMR_2))

####################################################################################################
# Fig. S7a  genomic interval between neighboring ends of sLRC ####
####################################################################################################
# The genomic interval (l) between two neighboring ends of conserved LRCs (299nm, in >23 chr.)
tgt = which(rowSums(conformation_matrix)>23)
distances_stat2_tmp = distances_stat2[tgt,]
rmv = which(distances_stat2_tmp[,2]==242000001)
distances_stat2_tmp = distances_stat2_tmp[-rmv,]
intvl2 = rep(0,length(spots2))
for (ss in 1:length(spots2)) {
  tmp = which(distances_stat2_tmp[,1]==spots2[ss]|distances_stat2_tmp[,2]==spots2[ss])
  if (length(tmp)!=0) {
    intvl2[ss] = 1
  }
}

A = diff(spots2[which(intvl2==1)])
# A = intrvl[which(intrvl!=250000)]
hist(log10(A), probability = TRUE, 
     main = "interval-conserved", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)

####################################################################################################
# Fig. 3a  distribution of l_ab ####
####################################################################################################
## find the threshhood of genomic distances (l) for clustering small links, the basic idea is:
# for a pair of small links, x and y, they have 4 genomic distances, from the shortest to longest are l1, l2, l3, and l4. 
## plot the distribution of l2 ####
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
print(dim(distances_stat22))
secdist = 0
lk1 = 0
lk2 = 0
for (ii in 1:(nrow(distances_stat22)-1)) {
  print(ii)
  tmp1 = as.numeric(distances_stat22[ii,])
  tmp2 = distances_stat22[(ii+1):nrow(distances_stat22),]
  
  l1 = as.numeric(tmp2[,1])-tmp1[1]
  l2 = as.numeric(tmp2[,2])-tmp1[1]
  l3 = abs(as.numeric(tmp2[,1])-tmp1[2])
  l4 = as.numeric(tmp2[,2])-tmp1[2]
  
  result_matrix = cbind(l1, l2, l3, l4)
  result_matrix = t(apply(result_matrix, 1, sort))
  secdist = c(secdist,result_matrix[,2])
  
  lk1 = c(lk1,rep(ii,length((ii+1):nrow(distances_stat22))))
  lk2 = c(lk2, ((ii+1):nrow(distances_stat22)))
}

secdist = secdist[-1]
lk1 = lk1[-1]
lk2 = lk2[-1]

set.seed(123)
fit <- mixtools::normalmixEM(secdist, k = 3)  # k specifies the number of components in the mixture
# Plot the histogram of the data

hist(secdist, breaks = 100, prob = TRUE,
     main = "The 2nd smallest l between pair of LRC", 
     xlab = "l (MB)", col = "orange", border = "brown")
# Overlay the Gaussian mixture distribution curve
x <- seq(min(secdist), max(secdist), length.out = 1000)
density_curve <- function(x_vals) {
  sum(sapply(1:length(fit$lambda), function(i) {
    fit$lambda[i] * dnorm(x_vals, mean = fit$mu[i], sd = fit$sigma[i])
  }))
}
densities <- sapply(x, density_curve)  # Calculate densities at each point in x
lines(x, densities, col = "purple", lwd = 2)  # Plot the density curve for the Gaussian mixture model
box(lwd = 3)
abline(v = 12000000, col = 'blue',lwd = 2,lty = 2)

####################################################################################################
# Fig. S7c = histogram of #sLRC per metaLRC ####
####################################################################################################
## apply the threshold =1MB of genomic distances (l) for clustering ####
### remove 242 MB first ####
bin_thresh = 1000000
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
conformation_matrix22 = conformation_matrix[tmp,]
tmp = which(distances_stat22[,2]!=242000001)
distances_stat22 = distances_stat22[tmp,]
conformation_matrix22 = conformation_matrix22[tmp,]
print(dim(distances_stat22))

library(dbscan)
epss = 250001
clusters <- dbscan(distances_stat22[,1:2], 
                   eps = epss, 
                   minPts = 1)
print(length(table(clusters$cluster)))
# Visualize the results
plot(distances_stat22[,1:2], col = clusters$cluster + 1, pch = 16, 
     main = "DBSCAN Clustering with Custom Connectivity",
     xlab = "X coordinate", 
     ylab = "Y coordinate")
binned_links2 = (cbind(clusters$cluster,distances_stat22))
colnames(binned_links2)[1] = 'group'

cluster_names = names(table(clusters$cluster))
binned_links2no242 = matrix(0,nrow = length(cluster_names),ncol = 5)
colnames(binned_links2no242) = c('st1','st2','ed1','ed2','numblinks')
for (ii in 1:length(cluster_names)) {
  tmp = which(clusters$cluster==ii)
  distances_stat22_tmp = distances_stat22[tmp,]
  binned_links2no242[ii,1] = min(distances_stat22_tmp[,1]) 
  binned_links2no242[ii,2] = max(distances_stat22_tmp[,1]) 
  binned_links2no242[ii,3] = min(distances_stat22_tmp[,2]) 
  binned_links2no242[ii,4] = max(distances_stat22_tmp[,2]) 
  binned_links2no242[ii,5] = length(tmp)
}

### add 242 MB back ####
bin_thresh = 1000000
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
conformation_matrix22 = conformation_matrix[tmp,]

tmp = which(distances_stat22[,2]==242000001) #remove 242 temporarily
distances_stat22 = distances_stat22[tmp,]
conformation_matrix22 = conformation_matrix22[tmp,]

clusters <- dbscan(distances_stat22[,1:2], 
                   eps = epss, 
                   minPts = 1)
print(length(table(clusters$cluster)))
binned_links2_tmp = (cbind(clusters$cluster+length(unique(distances_stat22[,1])),distances_stat22))
colnames(binned_links2_tmp)[1] = c('group')
binned_links2 = rbind(binned_links2,binned_links2_tmp)
# Visualize the results
plot(distances_stat22[,1:2], col = clusters$cluster + 1, pch = 16, 
     main = "DBSCAN Clustering with Custom Connectivity",
     xlab = "X coordinate", 
     ylab = "Y coordinate")

cluster_names = names(table(clusters$cluster))
binned_links2w242 = matrix(0,nrow = length(cluster_names),ncol = 5)
colnames(binned_links2w242) = c('st1','st2','ed1','ed2','numblinks')
for (ii in 1:length(cluster_names)) {
  tmp = which(clusters$cluster==ii)
  distances_stat22_tmp = distances_stat22[tmp,]
  binned_links2w242[ii,1] = min(distances_stat22_tmp[,1]) 
  binned_links2w242[ii,2] = max(distances_stat22_tmp[,1]) 
  binned_links2w242[ii,3] = min(distances_stat22_tmp[,2]) 
  binned_links2w242[ii,4] = max(distances_stat22_tmp[,2]) 
  binned_links2w242[ii,5] = length(tmp)
}

binned_links3 = binned_links2no242[which(binned_links2no242[,5]>1),]
print(dim(binned_links3))
binned_links3_tmp = binned_links2w242[which(binned_links2w242[,5]>1),]
print(dim(binned_links3_tmp))
binned_links3 = rbind(binned_links3, binned_links3_tmp)
print(dim(binned_links3))

# binned_links3 = cbind(binned_links3,binned_links3[,2]-binned_links3[,1],binned_links3[,4]-binned_links3[,3])
# colnames(binned_links3) = c('st1','st2','ed1','ed2','numblinks','interval_st','interval_ed')
# write.table(binned_links3,'metaLRC.txt',quote = F)

hist(binned_links2no242[,5],breaks =300,xlim = c(0,40))
# numblink = 10
hist(binned_links3[,5],breaks =200,freq = FALSE, xlim = c(0,100),col ="#cc3333",boader = 'white')
# abline(v = numblink,lty = 2, col = 'red')
# binned_links3 = binned_links3[which(binned_links3[,5]>=numblink),]

####################################################################################################
# Fig. S7f (occurrence of sLRC vs metaLRC) and ####
# Fig. 3i (alpha/gamma vs fitted delta/gamma for all metaLRC with > 7 sLRC) ####
####################################################################################################
## Calculate correlation among LRC in meta LRC #####
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
conformation_matrix22 = conformation_matrix[tmp,]

tmp = which(distances_stat22[,2]!=242000001)
distances_stat22 = distances_stat22[tmp,]
conformation_matrix22 = conformation_matrix22[tmp,]

numbLRC =0
binned_links22 = binned_links3[which(binned_links3[,5]>numbLRC),]

distances_stat222 = matrix(0,nrow = 1, ncol = ncol(binned_links2)+4)
linknumb = 0
maxlink = 0
for (ii in 1:nrow(binned_links22)) {
  tmp = as.numeric(binned_links22[ii,])
  st_tmp = which((distances_stat22[,1]>=tmp[1])&(distances_stat22[,1]<=tmp[2]))
  ed_tmp = which((distances_stat22[,2]>=tmp[3])&(distances_stat22[,2]<=tmp[4]))
  link_tmp = intersect(st_tmp,ed_tmp)
  print(length(link_tmp))
  linknumb = c(linknumb,length(link_tmp))
  conformation_matrix22_tmp = conformation_matrix22[link_tmp,]
  rownames(conformation_matrix22_tmp) =  1:length(link_tmp)
  maxlink = c(maxlink,length(which(colSums(conformation_matrix22_tmp)!=0)))
} 
linknumb = linknumb[-1]
maxlink = maxlink[-1]

hist(rowSums(conformation_matrix), freq = FALSE, col = 'lightgrey', main = "LRC v.s. metaLRC", xlab = "# chr", xlim = c(0, 1000))
hist(maxlink, freq = FALSE, breaks = 200, col = rgb(1, 0.65, 0, 0.5), boader = 'brown', add = TRUE)
# # Add density curves
# lines(density(rowSums(conformation_matrix)), col = "black", lwd = 2)
# lines(density(maxlink), col = "brown", lwd = 2)
# Add a legend
legend("topright", legend = c("single LRC", "metaLRC"), 
       fill = c('lightgrey', 'orange'))
box(lwd = 3)

# physical model-bin model
# a, alpha, association rate without any others' help;
# r, gamma, dissociation rate;
# d, delta, association rate with cooperativity binding
# while steady state:
# p1 = N*a/r * p0
# p2 = (N-1)(a+d)/2r *p1
# p3 = (N-2)(a+d)/3r *p2
# pN = (a+d)/Nr *pN-1

# Function to calculate fitted values
calculate_fitted <- function(ard, n_real) {
  n_fitted <- numeric(length(n_real))
  n_fitted[1] <- n_real[1]  # First value remains the same
  for (i in 2:length(n_real)) {
    xx = ((length(n_real)-i+1)* n_fitted[i - 1])/i
    n_fitted[i] <- xx*(arr+ard)
  }
  return(n_fitted)
}

# Objective function to minimize
objective_function <- function(ard, n_real) {
  n_fitted <- calculate_fitted(ard, n_real)
  return(sum((n_fitted - n_real)^2))
}

numbLRC =7
binned_links22 = binned_links3[which(binned_links3[,5]>numbLRC),]
distances_stat222 = matrix(0,nrow = 1, ncol = ncol(binned_links2)+4)
p0 = 0
p1 = 0
ar = 0
adr = 0
linknumb = 0
for (ii in 1:nrow(binned_links22)) {
  tmp = as.numeric(binned_links22[ii,])
  st_tmp = which((distances_stat22[,1]>=tmp[1])&(distances_stat22[,1]<=tmp[2]))
  ed_tmp = which((distances_stat22[,2]>=tmp[3])&(distances_stat22[,2]<=tmp[4]))
  link_tmp = intersect(st_tmp,ed_tmp)
  print(length(link_tmp))
  linknumb = c(linknumb,length(link_tmp))
  conformation_matrix22_tmp = conformation_matrix22[link_tmp,]
  NN = dim(conformation_matrix22_tmp)[1]
  NS = dim(conformation_matrix22_tmp)[2]

  p00 = length(which(colSums(conformation_matrix22_tmp)==0))/NS
  p11 =  length(which(colSums(conformation_matrix22_tmp)==1))/NS
  arr = p11/(NN*p00)
  n_real = length(which(colSums(conformation_matrix22_tmp)==1))/NS
  for (jj in 2:NN) {
    adrr_tmp = length(which(colSums(conformation_matrix22_tmp)==jj))/NS
    n_real <- c(n_real,adrr_tmp)
  }
  # Use optim to find the optimal d
  result <- optim(par = 0, fn = objective_function, n_real = n_real, method = "BFGS")
  # Extract the optimal d
  optimal_d <- result$par
  
  p0 = c(p0,p00)
  p1 = c(p1,p11)
  ar = c(ar,arr)
  adr = c(adr,optimal_d)
}
p0 = p0[-1]
p1 = p1[-1]
ar = ar[-1]
adr = adr[-1]
plot(ar,adr,pch = 19,xlim = c(0.00,0.008),ylim = c(0.001,0.2),col = 'blue')
abline(coef = c(0,1),lty = 2, col = 'purple')

####################################################################################################
# Fig. 3f (bound-sLRC distribution vs independent / cooperative model), ####
# Fig. 3e (phi heatmap of sLRC in one metaLRC), ####
# Fig. S7d,e ####
####################################################################################################
binned_links22 = binned_links2no242[which((binned_links2no242[,5]>8)&(binned_links2no242[,5]<20)&(binned_links2no242[,4]-binned_links2no242[,3]>0)),]
for (ii in 1:nrow(binned_links22)) {
  tmp = as.numeric(binned_links22[ii,])
  st_tmp = which((distances_stat22[,1]>=tmp[1])&(distances_stat22[,1]<=tmp[2]))
  ed_tmp = which((distances_stat22[,2]>=tmp[3])&(distances_stat22[,2]<=tmp[4]))
  link_tmp = intersect(st_tmp,ed_tmp)
  print(length(link_tmp))
  # linknumb = c(linknumb,length(link_tmp))
  conformation_matrix22_tmp = conformation_matrix22[link_tmp,]
  NN = dim(conformation_matrix22_tmp)[1]
  NS = dim(conformation_matrix22_tmp)[2]
  
  p00 = length(which(colSums(conformation_matrix22_tmp)==0))/NS
  p11 =  length(which(colSums(conformation_matrix22_tmp)==1))/NS
  arr = p11/(NN*p00)
  adrr = 0
  pp = p11
  n_real = length(which(colSums(conformation_matrix22_tmp)==1))/NS
  for (jj in 2:NN) {
    adrr_tmp = length(which(colSums(conformation_matrix22_tmp)==jj))/NS
    n_real <- c(n_real,adrr_tmp)
  }
  # Use optim to find the optimal d
  result <- optim(par = 0, fn = objective_function, n_real = n_real, method = "BFGS")
  # Extract the optimal d
  optimal_d <- result$par
  
  # Calculate the fitted values using the optimal d
  n_fitted <- calculate_fitted(optimal_d, n_real)
  
  pd0 = c(p00,p11)
  for (jj in 2:NN) {
    tmp = ((NN-jj+1)*arr*pd0[jj])/jj
    pd0 = c(pd0,tmp)
  }
  
  # Create a data frame
  df <- data.frame(
    NumbBinds = c(0,1:NN),
    Real = c(p00,n_real),
    No_cor = pd0,
    With_cor = c(p00,n_fitted)
  )
  
  df = df[-1,]
  
  # Reshape the data for ggplot
  df_melted <- reshape2::melt(df, id.vars = "NumbBinds", variable.name = "Set", value.name = "Value")
  
  # Create the plot
  p <- ggplot(df_melted, aes(x = NumbBinds, y = Value, fill = Set)) +
    geom_bar(data = subset(df_melted, Set %in% c("Real", "No_cor")), 
             stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    geom_line(data = subset(df_melted, Set == "With_cor"), 
              aes(group = 1, color = "With_cor"), size = .8,linetype = 'dashed') +
    geom_point(data = subset(df_melted, Set == "With_cor"), 
               aes(color = "With_cor"), size = 1) +
    scale_fill_manual(values = c("Real" = "purple", "No_cor" = "darkgoldenrod3")) +
    scale_color_manual(values = c("With_cor" = "blue")) +
    labs(title = paste(ii,NN,':',(binned_links22[ii,1]-1)/1000000,(binned_links22[ii,2]-1)/1000000,"-",(binned_links22[ii,3]-1)/1000000,(binned_links22[ii,4]-1)/1000000),
         x = "Number of binds",
         y = "% of chr") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    theme(legend.position = "bottom") + 
    # scale_y_break(c(0.02, 0.8)) +
    guides(fill = guide_legend(title = "Bar Sets"),
           color = guide_legend(title = "Line Set"))
  
  # Display the plot
  print(p)
  
  # pheatmap(conformation_matrix22_tmp,color = c('white','blue'),main =paste(ii,tmp[1],tmp[2],tmp[3],tmp[4]))
  xxx = phi(((conformation_matrix22_tmp)))
  rownames(xxx) = 1:length(link_tmp)
  colnames(xxx) = 1:length(link_tmp)
  diag(xxx) = NA
  pheatmap(xxx,
           # color = colorRampPalette(c("blue", "white", "red"))(60),
           # breaks = seq(0, 0.6, length.out = 61),
           main =paste(ii,NN,':',(binned_links22[ii,1]-1)/1000000,(binned_links22[ii,2]-1)/1000000,"-",(binned_links22[ii,3]-1)/1000000,(binned_links22[ii,4]-1)/1000000))
}

####################################################################################################
# Fig. 3h  distance between the two ends of a metaLRC (minimal bounding spheres) vs number of bound sLRC per chr2 trace (violin) ####
####################################################################################################
## plot 3D of certain multivalaent binding ####
# Function to calculate the minimal bounding sphere
calculate_minimal_sphere <- function(points) {
  points = as.matrix(points,ncol = 4)
  points = matrix(points[,1:3],ncol = 3)
  # Calculate the centroid
  if (nrow(points)==1) {
    center = points
    radius = 0
  }else if (nrow(points)==2) {
    center <- colMeans(points)
    distances <-  sqrt(rowSums((t(points[1,]) - t(points[2,]))^2))
    radius <- (distances)
  }
  else{
    center <- colMeans(points)
    # Calculate distances from each point to the centroid
    distances <- sqrt(rowSums((t(points) - center)^2))
    # Find the point furthest from the centroid
    max_distance <- (max(distances))
    furthest_point <- points[which.max(distances), ]
    # Recalculate center as midpoint between centroid and furthest point
    center <- (center + furthest_point) / 2
    # Recalculate distances from each point to the new center
    distances <- sqrt(rowSums((t(points) - center)^2))
    # The radius is the maximum of these distances
    radius <- max(distances)
  }
  return(list(center = center, radius = radius))
}

rsl = 250
chrom <- (1:(round((chr2+2000000)/2000000)+1))*2-2
start <- rep(0,(round((chr2+2000000)/2000000)+1))
end <- rep(2000000/(rsl*1000),(round((chr2+2000000)/2000000)+1))
chr_sizes_df <- data.frame(chrom,start,end)

head(binned_links3)
head(distances_stat2)
print(dim(distances_stat2))

cells = unique(distances_2m[,2])
conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

tmp = which(distances_stat2$numblks>23)
distances_stat22 = distances_stat2[tmp,]
print(dim(distances_stat22))
conformation_matrix22 = conformation_matrix[tmp,]
print(dim(conformation_matrix22))

binned_conformation = matrix(0, nrow = nrow(binned_links3), ncol = ncol(conformation_matrix))
for (bl in 1:(nrow(binned_links3))) {
  print(bl)
  tmp1 = as.numeric(binned_links3[bl,])
  st_tmp1 = which((distances_stat22[,1]>=tmp1[1])&(distances_stat22[,1]<=tmp1[2]))
  ed_tmp1 = which((distances_stat22[,2]>=tmp1[3])&(distances_stat22[,2]<=tmp1[4]))
  link_tmp1 = intersect(st_tmp1,ed_tmp1)
  conformation_matrix_tmp1 = conformation_matrix22[link_tmp1,]
  binned_conformation[bl,] = colSums(conformation_matrix_tmp1)
}
binned_conformation2 = binned_conformation
binned_conformation2[which(binned_conformation2!=0,arr.ind = )] = 1
hist(rowSums(binned_conformation2))
binned_links33 = cbind(binned_links3,rowSums(binned_conformation2))
head(binned_links33)
plot(binned_links33[,5],binned_links33[,6])

# calculate distance distribution of metaLRC
# binned_links333 = binned_links33[which((binned_links33[,5]<10)&(binned_links33[,6]>150)),]# the chosen metaLRCs have at most 10 sLRC and show in at least 150 chr copies. 
print(dim(binned_links333))
for (bl in c(1)) {
  print(bl)
  tmp1 = as.numeric(binned_links333[bl,])
  st_tmp1 = which((distances_stat22[,1]>=tmp1[1])&(distances_stat22[,1]<=tmp1[2]))
  ed_tmp1 = which((distances_stat22[,2]>=tmp1[3])&(distances_stat22[,2]<=tmp1[4]))
  link_tmp1 = intersect(st_tmp1,ed_tmp1)
  conformation_matrix_tmp1 = conformation_matrix22[link_tmp1,]
  distances_stat22_tmp1 = distances_stat22[link_tmp1,]
  
  st_tmp2 = which((distances_stat2[,1]>=tmp1[1])&(distances_stat2[,1]<=tmp1[2]))
  ed_tmp2 = which((distances_stat2[,2]>=tmp1[3])&(distances_stat2[,2]<=tmp1[4]))
  link_tmp2 = intersect(st_tmp2,ed_tmp2)
  conformation_matrix_tmp2 = conformation_matrix[link_tmp2,]
  distances_stat22_tmp2 = distances_stat2[link_tmp2,]
  
  dist_locipair = matrix(0,nrow = 3000,ncol = 5)
  for (cc in 1:3000) {
    tmp = which(IMR_2$Chromosome.copy.number==cc)
    loci1 = which(IMR_2$Coordinate[tmp]%in%distances_stat22_tmp1[,1])
    dist_locipair[cc,1] = length(loci1)
    points1 <- data.frame(
      x = IMR_2$X.nm.[tmp][loci1],
      y = IMR_2$Y.nm.[tmp][loci1],
      z = IMR_2$Z.nm.[tmp][loci1],
      cord = IMR_2$Coordinate[tmp][loci1]
    )
    points1 = na.omit(points1)
    if (nrow(points1)>0) {
      result1 <- calculate_minimal_sphere(points1)
    }
    
    loci2 = which(IMR_2$Coordinate[tmp]%in%distances_stat22_tmp1[,2])
    dist_locipair[cc,2] = length(loci2)
    points2 <- data.frame(
      x = IMR_2$X.nm.[tmp][loci2],
      y = IMR_2$Y.nm.[tmp][loci2],
      z = IMR_2$Z.nm.[tmp][loci2],
      cord = IMR_2$Coordinate[tmp][loci2]
    )
    points2 = na.omit(points2)
    if (nrow(points2)>0) {
      result2 <- calculate_minimal_sphere(points2)
    }
    dist_locipair[cc,3] = sqrt(rowSums((t(as.numeric(result1$center)) - t(as.numeric(result2$center)))^2))
    dist_locipair[cc,4] = as.numeric(result1$radius)
    dist_locipair[cc,5] = as.numeric(result2$radius)
  }
  data = data.frame(sLRC_numb = colSums(conformation_matrix_tmp1),
                    LRC_numb = colSums(conformation_matrix_tmp2),
                    dist = dist_locipair[,3]
  )
}
for (bl in 2:(nrow(binned_links333))) {
  print(bl)
  tmp1 = as.numeric(binned_links333[bl,])
  st_tmp1 = which((distances_stat22[,1]>=tmp1[1])&(distances_stat22[,1]<=tmp1[2]))
  ed_tmp1 = which((distances_stat22[,2]>=tmp1[3])&(distances_stat22[,2]<=tmp1[4]))
  link_tmp1 = intersect(st_tmp1,ed_tmp1)
  conformation_matrix_tmp1 = conformation_matrix22[link_tmp1,]
  distances_stat22_tmp1 = distances_stat22[link_tmp1,]
  
  st_tmp2 = which((distances_stat2[,1]>=tmp1[1])&(distances_stat2[,1]<=tmp1[2]))
  ed_tmp2 = which((distances_stat2[,2]>=tmp1[3])&(distances_stat2[,2]<=tmp1[4]))
  link_tmp2 = intersect(st_tmp2,ed_tmp2)
  conformation_matrix_tmp2 = conformation_matrix[link_tmp2,]
  distances_stat22_tmp2 = distances_stat2[link_tmp2,]
  
  dist_locipair = matrix(0,nrow = 3000,ncol = 5)
  for (cc in 1:3000) {
    tmp = which(IMR_2$Chromosome.copy.number==cc)
    loci1 = which(IMR_2$Coordinate[tmp]%in%distances_stat22_tmp1[,1])
    dist_locipair[cc,1] = length(loci1)
    points1 <- data.frame(
      x = IMR_2$X.nm.[tmp][loci1],
      y = IMR_2$Y.nm.[tmp][loci1],
      z = IMR_2$Z.nm.[tmp][loci1],
      cord = IMR_2$Coordinate[tmp][loci1]
    )
    points1 = na.omit(points1)
    if (nrow(points1)>0) {
      result1 <- calculate_minimal_sphere(points1)
    }
    
    loci2 = which(IMR_2$Coordinate[tmp]%in%distances_stat22_tmp1[,2])
    dist_locipair[cc,2] = length(loci2)
    points2 <- data.frame(
      x = IMR_2$X.nm.[tmp][loci2],
      y = IMR_2$Y.nm.[tmp][loci2],
      z = IMR_2$Z.nm.[tmp][loci2],
      cord = IMR_2$Coordinate[tmp][loci2]
    )
    points2 = na.omit(points2)
    if (nrow(points2)>0) {
      result2 <- calculate_minimal_sphere(points2)
    }
    dist_locipair[cc,3] = sqrt(rowSums((t(as.numeric(result1$center)) - t(as.numeric(result2$center)))^2))
    dist_locipair[cc,4] = as.numeric(result1$radius)
    dist_locipair[cc,5] = as.numeric(result2$radius)
  }
  data_tmp = data.frame(sLRC_numb = colSums(conformation_matrix_tmp1),
                        LRC_numb = colSums(conformation_matrix_tmp2),
                        dist = dist_locipair[,3]
  )
  data = rbind(data_tmp,data)
}

data$sLRC_numb = factor(data$sLRC_numb, levels = c(0,1:10))
data$sLRC_numb = factor(data$sLRC_numb, levels = c(0,1:10))

plot(dist_locipair[,3],colSums(conformation_matrix_tmp2))
ggplot(data, aes(x = sLRC_numb, y = dist, fill = sLRC_numb)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "sLRC_numb v.s distances",
       x = "sLRC_numb",
       y = "distance") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.position = "none")

ggplot(data, aes(x = LRC_numb, y = dist, fill = LRC_numb)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "LRC_numb v.s distances",
       x = "LRC_numb",
       y = "distance") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.position = "none")

####################################################################################################
# Fig. 3g  3D structures (plotly, saved as metaLRC_<trace>-metaLRC_<bl>.html) of the example metaLRC ####
####################################################################################################
five_cols = c('#1d0ef7','#a90ef7','#db1501','#81b108','#008B8B','#FFCC00')
# Fig. 3g shows representative metaLRC (e.g.: 28-29.75 Mb & 216.25-216.5 Mb).
binned_links333 = binned_links33[which((binned_links33[,1]==28000001)&(binned_links33[,2]==29750001)&(binned_links33[,3]==216250001)&(binned_links33[,4]==216500001)),]
for (bl in 1:(nrow(binned_links333))) {
  print(bl)
  tmp1 = as.numeric(binned_links333[bl,])
  st_tmp1 = which((distances_stat22[,1]>=tmp1[1])&(distances_stat22[,1]<=tmp1[2]))
  ed_tmp1 = which((distances_stat22[,2]>=tmp1[3])&(distances_stat22[,2]<=tmp1[4]))
  link_tmp1 = intersect(st_tmp1,ed_tmp1)
  conformation_matrix_tmp1 = conformation_matrix22[link_tmp1,]
  distances_stat22_tmp1 = distances_stat22[link_tmp1,]
  
  st_tmp2 = which((distances_stat2[,1]>=tmp1[1])&(distances_stat2[,1]<=tmp1[2]))
  ed_tmp2 = which((distances_stat2[,2]>=tmp1[3])&(distances_stat2[,2]<=tmp1[4]))
  link_tmp2 = intersect(st_tmp2,ed_tmp2)
  conformation_matrix_tmp2 = distances_matrix[link_tmp2,]
  distances_stat22_tmp2 = distances_stat2[link_tmp2,]
  
  LC = unique(c(distances_stat22_tmp1[,1],distances_stat22_tmp1[,2]))
  smpl = 0
  for (ll in LC) {
    smpl_tmp = which((distances_stat22_tmp1[,1]==ll)|(distances_stat22_tmp1[,2]==ll))
    smpl_tmp = which(matrix(conformation_matrix_tmp1[smpl_tmp,],ncol = ncol(conformation_matrix_tmp1))>0,arr.ind = T)[,2]
    smpl_tmp = sample(smpl_tmp,2)
    smpl = c(smpl, smpl_tmp)
  }
  smpl = smpl[-1]
  # smpl = sample(smpl,10)2972,328,1043
  smpl = c(2040,2988,1976,2008,1726)
  print(smpl)
  for (ss in smpl) {
    print(ss)
    # print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
    distances_stat_tmp1 = cbind(distances_stat22_tmp1,as.numeric(as.character(conformation_matrix_tmp1[,ss])))
    print(distances_stat_tmp1)
    
    distances_stat_tmp2 = cbind(distances_stat22_tmp2,as.numeric(as.character(conformation_matrix_tmp2[,ss])))
    print(dim(distances_stat_tmp2))
    distances_stat_tmp3 = distances_stat_tmp2[which(distances_stat_tmp2[,5]<d3li),]
    print(distances_stat_tmp3)
    
    # # linkage map
    # circos.clear()
    # circos.par(cell.padding=c(0,0,0,0),
    #            track.margin=c(0,0.05),
    #            start.degree = 90,
    #            gap.degree = 0.1,
    #            clock.wise = TRUE)
    # circos.initialize(factors = chr_sizes_df$chrom,
    #                   xlim = cbind(chr_sizes_df$start, chr_sizes_df$end))
    # 
    # circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    #   chr=CELL_META$sector.index
    #   xlim=CELL_META$xlim
    #   ylim=CELL_META$ylim
    #   circos.text(mean(xlim), mean(ylim), chr, cex = 0.5)
    # },bg.col="#d3e5ed", bg.border=TRUE, track.height=0.1)
    # 
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
    #   xleft = 2,
    #   xright = 8,
    #   col = 'white',
    #   border = NA,
    #   ytop = 1.5,
    #   ybottom=-0.5
    # )
    # title(paste0('metaLRC: ',bl,'; linkage map: ',ss), cex.main = 1.5)
    # locs = seq(0,(chr2+2000000),2000000)
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp1[,1]>=loc_st) & (distances_stat_tmp1[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp1[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'grey',lwd = .2)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp2[,1]>=loc_st) & (distances_stat_tmp2[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp2[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'black',lwd = 0.4)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    distances_stat_tmp22 = distances_stat_tmp1[which(distances_stat_tmp1[,5]!=0),]
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp22[,1]>=loc_st) & (distances_stat_tmp22[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp22[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = five_cols[bl],lwd = 0.5)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    
    ## 3D plot ####
    
  smpl = c(1346,564,1541,2734) # for example
  print(smpl)
  for (ss in smpl) {
    print(ss)
    # print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
    distances_stat_tmp1 = cbind(distances_stat22_tmp1,as.numeric(as.character(conformation_matrix_tmp1[,ss])))
    print(distances_stat_tmp1)
    
    distances_stat_tmp2 = cbind(distances_stat22_tmp2,as.numeric(as.character(conformation_matrix_tmp2[,ss])))
    print(dim(distances_stat_tmp2))
    distances_stat_tmp3 = distances_stat_tmp2[which(distances_stat_tmp2[,5]<d3li),]
    print(distances_stat_tmp3)
    
    # # linkage map
    # circos.clear()
    # circos.par(cell.padding=c(0,0,0,0),
    #            track.margin=c(0,0.05),
    #            start.degree = 90,
    #            gap.degree = 0.1,
    #            clock.wise = TRUE)
    # circos.initialize(factors = chr_sizes_df$chrom,
    #                   xlim = cbind(chr_sizes_df$start, chr_sizes_df$end))
    # 
    # circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    #   chr=CELL_META$sector.index
    #   xlim=CELL_META$xlim
    #   ylim=CELL_META$ylim
    #   circos.text(mean(xlim), mean(ylim), chr,cex = 0.5)
    # },bg.col="#d3e5ed", bg.border=TRUE, track.height=0.1)
    # 
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
    #   xleft = 2,
    #   xright = 8,
    #   col = 'white',
    #   border = NA,
    #   ytop = 1.5,
    #   ybottom=-0.5
    # )
    # title(paste0('metaLRC: ',bl,'; linkage map: ',ss), cex.main = 1.5)
    # locs = seq(0,(chr2+2000000),2000000)
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp1[,1]>=loc_st) & (distances_stat_tmp1[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp1[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'grey',lwd = .2)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp2[,1]>=loc_st) & (distances_stat_tmp2[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp2[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = 'black',lwd = 0.4)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    distances_stat_tmp22 = distances_stat_tmp1[which(distances_stat_tmp1[,5]!=0),]
    # for (ii in 1:(length(locs)-1)) {
    #   loc_st = locs[ii]
    #   loc_ed = locs[ii+1]
    #   locst = seq(loc_st,loc_ed,rsl*1000)
    #   tmp1 = which((distances_stat_tmp22[,1]>=loc_st) & (distances_stat_tmp22[,1]<loc_ed))
    #   if (length(tmp1)>0) {
    #     tmp1 = matrix(as.matrix(distances_stat_tmp22[tmp1,]),ncol=5)
    #     for (jj in (ii+1):122) {
    #       # for (jj in 121:121) {  # create the small links for 240-242 MB locus
    #       loc_st = locs[jj]
    #       loc_ed = locs[jj+1]
    #       loced = seq(loc_st,loc_ed,rsl*1000)
    #       tmp2 = which((as.numeric(tmp1[,2])>=loc_st) & ((as.numeric(tmp1[,2]))<=loc_ed))
    #       if(length(tmp2)!=0){
    #         tmp2 = matrix(tmp1[tmp2,],ncol=5)
    #         for (kk in 1:length(locst)) {
    #           tmpp = which((tmp2[,1]-1)==locst[kk])
    #           if(length(tmpp)!=0){
    #             tmpp = matrix(tmp2[tmpp,],ncol=5)
    #             for (kkk in 1:length(loced)) {
    #               tmppp = which((tmpp[,2]-1)==loced[kkk])
    #               if(length(tmppp)!=0){
    #                 lwdd = (10/tmpp[tmppp,4])
    #                 if(tmppp!=0){
    #                   circos.link(ii*2-2, kk, jj*2-2, kkk, col = five_cols[bl],lwd = 0.5)
    #                 }
    #               }
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }
    
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
    data$link_group <- "Other"
    RC_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp2[,1]),
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp2[,2])
    ))
    data$link_group[RC_indices] <- "RC"
    LRC_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp3[,1]),
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp3[,2])
    ))
    data$link_group[LRC_indices] <- "LRC"
    sLRC_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp22[,1]),
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp22[,2])
    ))
    data$link_group[LRC_indices] <- "sLRC"
    
    # Create a group column in the data
    data$locus_group <- "Other"
    lc1_group_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp2[,1])
    ))
    data$locus_group[lc1_group_indices] <- "LR-Loci-1"
    all_colors[lc1_group_indices] = 'pink'
    lc2_group_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp2[,2])
    ))
    data$locus_group[lc2_group_indices] <- "LR-Loci-2"
    all_colors[lc2_group_indices] = 'lightgreen'
    
    l1_group_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp1[,1])
    ))
    data$locus_group[l1_group_indices] <- "Loci-1"
    all_colors[l1_group_indices] = 'red'
    l2_group_indices <- unique(c(
      which(IMR_2$Coordinate[tmp] %in% distances_stat_tmp1[,2])
    ))
    data$locus_group[l2_group_indices] <- "Loci-2"
    all_colors[l2_group_indices] = 'green'
    
    ## Create the base plot for html ####
    plot <- plot_ly()%>%
      layout(
        title = paste0("chr: ",ss,"; metaLRC: ",bl))
    
    # Add all points
    plot <- add_trace(plot,
                      data = data, 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers+lines",
                      marker = list(size = 7, color = all_colors, symbol = "circle",opacity = 0.3),
                      line = list(color = "lightgrey", width = 1),
                      name = "Other",
                      visible = TRUE,
                      text = data$cord,  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0))
    # Add two ends points of LR in the regions
    plot <- add_trace(plot,
                      data = subset(data, locus_group == "LR-Loci-1"), 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 7, color = all_colors[lc1_group_indices], symbol = "circle"),
                      name = "LR-End-1",
                      legendgroup = "LR-End-1",
                      showlegend = T,
                      visible = TRUE,
                      text = data$cord[which(data$locus_group=="LR-Loci-1")],  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0))
    
    plot <- add_trace(plot,
                      data = subset(data, locus_group == "LR-Loci-2"), 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 7, color = 'lightgreen', symbol = "circle"),
                      name = "LR-End-2",
                      legendgroup = "LR-End-2",
                      showlegend = T,
                      visible = TRUE,
                      text = data$cord[which(data$locus_group=="LR-Loci-2")],  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0))
    
    # Add two ends points of metaLRC
    plot <- add_trace(plot,
                      data = subset(data, locus_group == "Loci-1"), 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 7, color = all_colors[l1_group_indices], symbol = "circle"),
                      name = "End-1",
                      legendgroup = "End-1",
                      showlegend = T,
                      visible = TRUE,
                      text = data$cord[which(data$locus_group=="Loci-1")],  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0))
    
    plot <- add_trace(plot,
                      data = subset(data, locus_group == "Loci-2"), 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 7, color = 'green', symbol = "circle"),
                      name = "End-2",
                      legendgroup = "End-2",
                      showlegend = T,
                      visible = TRUE,
                      text = data$cord[which(data$locus_group=="Loci-2")],  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0))
    
    # Add LR links
    L1 = unique(distances_stat_tmp2[,1])
    for (pp in 1:(length(L1)-1)) {
      if ((L1[pp+1]-L1[pp])==250000) {
        pp1 <- which(IMR_2$Coordinate[tmp] == L1[pp])
        pp2 <- which(IMR_2$Coordinate[tmp] == L1[pp+1])
        
        plot <- add_trace(plot,
                          x = c(data$x[pp1], data$x[pp2]),
                          y = c(data$y[pp1], data$y[pp2]),
                          z = c(data$z[pp1], data$z[pp2]),
                          type = "scatter3d",
                          mode = "lines",
                          line = list(color = "lightgrey", width = 1),
                          name = "LR",
                          legendgroup = "LR",
                          showlegend = (pp == 1),
                          visible = TRUE)
      }
    }
    L1 = unique(distances_stat_tmp2[,2])
    for (pp in 1:(length(L1)-1)) {
      if ((L1[pp+1]-L1[pp])==250000) {
        pp1 <- which(IMR_2$Coordinate[tmp] == L1[pp])
        pp2 <- which(IMR_2$Coordinate[tmp] == L1[pp+1])
        
        plot <- add_trace(plot,
                          x = c(data$x[pp1], data$x[pp2]),
                          y = c(data$y[pp1], data$y[pp2]),
                          z = c(data$z[pp1], data$z[pp2]),
                          type = "scatter3d",
                          mode = "lines",
                          dash = "dash",
                          line = list(color = "lightgrey", width = 1),
                          name = "LR",
                          legendgroup = "LR",
                          showlegend = F,
                          visible = TRUE)
      }
    }
    # Add LRC links
    for (pp in 1:nrow(distances_stat_tmp3)) {
      pp1 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp3[pp,1])
      pp2 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp3[pp,2])
      
      plot <- add_trace(plot,
                        x = c(data$x[pp1], data$x[pp2]),
                        y = c(data$y[pp1], data$y[pp2]),
                        z = c(data$z[pp1], data$z[pp2]),
                        type = "scatter3d",
                        mode = "lines",
                        line = list(color = "darkgrey", width = 3),
                        name = "LRC",
                        legendgroup = "LRC",
                        showlegend = (pp == 1),
                        visible = TRUE)
    }
    # Add sLRC links
    for (pp in 1:nrow(distances_stat_tmp22)) {
      pp1 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp22[pp,1])
      pp2 <- which(IMR_2$Coordinate[tmp] == distances_stat_tmp22[pp,2])
      
      plot <- add_trace(plot,
                        x = c(data$x[pp1], data$x[pp2]),
                        y = c(data$y[pp1], data$y[pp2]),
                        z = c(data$z[pp1], data$z[pp2]),
                        type = "scatter3d",
                        mode = "lines",
                        line = list(color = "black", width = 3),
                        name = "sLRC",
                        legendgroup = "sLRC",
                        showlegend = (pp == 1),
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
        xaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = T, title = ""),
        yaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = T, title = ""),
        zaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = T, title = "")
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
    
    plot <- plot %>% layout(
      scene = list(
        camera = list(
          up = list(x = 0, y = 0, z = 1),
          center = list(x = 0, y = 0, z = 0),
          eye = list(x = 1.25, y = 1.25, z = 1.25)
        ),
        dragmode = "turntable"
      )
    ) %>%
      onRender("
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
    htmlwidgets::saveWidget(plot,paste0("metaLRC_",ss,"-metaLRC_",bl,".html"))
    print(paste0("metaLRC_",ss,"-metaLRC_",bl,".html"))
  }
}
