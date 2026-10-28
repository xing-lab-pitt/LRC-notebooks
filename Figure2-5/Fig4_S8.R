# =====================================================================================================
# Fig. 4 and Fig. S8
# =====================================================================================================
# INPUTS:   NMuMG_rep1_pos0..9.csv  (seqFISH+ spots incl. one column per epigenetic mark, z-scored)
#           41586_2025_8838_MOESM5_ESM.csv (probe table)
#           dist_NMuMG_rep1_pos*-chr1-allele*_100Mb.csv 
#           4DNFIHM89EGL.bw
#           plus chromosome2.tsv / distances_chr2L100-3000.txt
#           mesc_pairs/pairs/*.pairs.gz                    single-cell pairs (GEO GSM9115533)
#           GSM9115533_mesc.RNAcounts.gene.b6.format.tsv.gz / .cast.format.tsv.gz   allele-specific counts (gene x cell)
#           mouse_gene_annotation.rds  
####################################################################################################
library(geometry)
library(dplyr) 
library(ggplot2)
library(circlize)
options(scipen = 999)  
library(gridExtra)
library(plotly)

## DNA seqFISH data clean up####
probes = read.csv('41586_2025_8838_MOESM5_ESM.csv')
probes_ID = strsplit(probes$name, split = '-')
probes_ID = do.call(rbind,probes_ID)
probes = cbind(probes_ID,probes$name,probes$chrom,probes$start)
colnames(probes) = c('chr','probe_ID','names','chrom','start')
probes = data.frame(probes)

## LRC and epigenetics ####
### clean ####
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
for (cc in paste0('chr',c(1:2))) {
  loci_pair_LRC = matrix(0, nrow = 1, ncol = 3)
  colnames(loci_pair_LRC) = c('loc1','loc2','cell')
  for (file_name in file_names) { 
    for (allel in c('1','2')) {
      fls = paste0("ist_", file_name, '-', cc, '-allele',allel,"_100Mb.csv")
      print(fls)
      seqFISH_dist <- read.csv(fls)
      print(dim(seqFISH_dist))
      
      cell_ID = colnames(seqFISH_dist[,3:ncol(seqFISH_dist)])
      cell_ID = strsplit(cell_ID,split = '.', fixed = TRUE)
      cell_ID = do.call(rbind,cell_ID)
      
      LRC = which(seqFISH_dist[,3:ncol(seqFISH_dist)]<299,arr.ind = T)
      loci_pair_LRC_tmp =cbind(seqFISH_dist[LRC[,1],1], seqFISH_dist[LRC[,1],2],paste0(file_name,'-allel',allel,'-',cell_ID[LRC[,2],2]))
      
      loci_pair_LRC = rbind(loci_pair_LRC, loci_pair_LRC_tmp)
    }
  }
  loci_pair_LRC = loci_pair_LRC[-1,]
  print(length(loci_pair_LRC))
  write.csv(loci_pair_LRC, file = paste0("LRC_", cc,".csv"), row.names = FALSE)
}

for (cc in paste0('chr',c(1:2))) {
  probes_tmp = which(probes$chrom==cc)
  probes_tmp = probes[probes_tmp,c(1,2,5)]
  
  fls = paste0("LRC_", cc,".csv")
  print(fls)
  loci_pair_LRC = read.csv(fls)
  split_list <- strsplit(loci_pair_LRC$cell, "-")
  cells = do.call(rbind, split_list)
  loci_pair_LRC2 = cbind(loci_pair_LRC$loc1,loci_pair_LRC$loc2,cells)
  colnames(loci_pair_LRC2) = c('loc1', 'loc2','pos','allel','cell')
  loci_pair_LRC2 = data.frame(loci_pair_LRC2)

  lp_cols = c('loc1','loc2','pos','allel','cell')
  sd_cols = c('loc', colnames(read.csv(paste0(file_names[1], '.csv'), nrows = 1)))

  loci_pair_LRC3 = matrix(rep(0,length(c(lp_cols,'cc',paste0('loc1-',sd_cols),paste0('loc2-',sd_cols)))),nrow = 1)
  colnames(loci_pair_LRC3) = c(lp_cols,'cc',paste0('loc1-',sd_cols),paste0('loc2-',sd_cols))
  for (file_name in file_names) {
    fls = paste0(file_name,'.csv')
    print(fls)
    seqFISH_DNA <- read.csv(fls)
    tmp = which(seqFISH_DNA$chrom==cc)
    seqFISH_DNA = seqFISH_DNA[tmp,]
    
    loc_tmp = strsplit(seqFISH_DNA$name, "-")
    loc_tmp = do.call(rbind,loc_tmp)
    tmp = match(loc_tmp[,2],probes_tmp$probe_ID)
    seqFISH_DNA2 = cbind(probes_tmp$start[tmp],seqFISH_DNA)
    colnames(seqFISH_DNA2)[1] = 'loc'
    
    loci_pair_LRC_tmp = loci_pair_LRC2[which(loci_pair_LRC2$pos==file_name),]
    print(dim(loci_pair_LRC_tmp))

    for (ii in 1:nrow(loci_pair_LRC_tmp)) {
      if(ii%%500==0){
        print(ii)
      }
      allel_tmp = sub(".*(\\d)$", "\\1", loci_pair_LRC_tmp$allel[ii])
      tmp1 = which((seqFISH_DNA2$dbscan_allele==allel_tmp)&(seqFISH_DNA2$cellID==loci_pair_LRC_tmp$cell[ii])&((as.numeric(seqFISH_DNA2$loc)/1000)==loci_pair_LRC_tmp$loc1[ii]))[1]
      tmp2 = which((seqFISH_DNA2$dbscan_allele==allel_tmp)&(seqFISH_DNA2$cellID==loci_pair_LRC_tmp$cell[ii])&((as.numeric(seqFISH_DNA2$loc)/1000)==loci_pair_LRC_tmp$loc2[ii]))[1]
      tmp = c(loci_pair_LRC_tmp[ii,],cc,seqFISH_DNA2[tmp1,],seqFISH_DNA2[tmp2,])
      loci_pair_LRC3 = rbind(loci_pair_LRC3,tmp)
    }
  }
  loci_pair_LRC3 = loci_pair_LRC3[-1,]
  write.csv(loci_pair_LRC3, file = paste0("LRC_", cc,"-2.csv"), row.names = FALSE)
}
loci_pair_LRC3  = data.frame(loci_pair_LRC3)

####################################################################################################
# Fig. 4c-f
####################################################################################################
cc = 'chr1'
fls = paste0("LRC_", cc,"-2.csv")
print(fls)
loci_pair_LRC = read.csv(fls)
print(dim(loci_pair_LRC))
tmp = paste(loci_pair_LRC$pos,loci_pair_LRC$allel,loci_pair_LRC$cell)
cells = unique(tmp)
print(length(cells))

### grouped to 250kb ####
loci_bin = seq(0,198000000,250000)
loci_pair_LRC$posA_bin <- loci_bin[findInterval(loci_pair_LRC$loc1*1000, loci_bin)]
loci_pair_LRC$posB_bin <- loci_bin[findInterval(loci_pair_LRC$loc2*1000, loci_bin)]
loci_pair_LRC$posA_bin_tmp = paste0(loci_pair_LRC$posA_bin,"_",loci_pair_LRC$pos,"_",loci_pair_LRC$allel,"_",loci_pair_LRC$cell,"_",loci_pair_LRC$cc)
loci_pair_LRC$posB_bin_tmp = paste0(loci_pair_LRC$posB_bin,"_",loci_pair_LRC$pos,"_",loci_pair_LRC$allel,"_",loci_pair_LRC$cell,"_",loci_pair_LRC$cc)
epi_marks = c('H3K4me2','H3K27ac','SF3A66','H3K9me3','H3K27me3','EZH2','LaminB1')  
nselct = c("loc1","loc2","posA_bin_tmp", "posB_bin_tmp",paste0("loc1.",epi_marks),paste0("loc2.",epi_marks))
loci_pair_LRC2 = loci_pair_LRC[,nselct]

result_df_dplyr <- loci_pair_LRC2 %>%
  group_by(posA_bin_tmp, posB_bin_tmp) %>%  
  summarise(
    across(
      .cols = -c(loc1, loc2),      
      .fns = ~median(.x, na.rm = TRUE),
      .names = "{.col}_md"        
    ),
    .groups = "drop"  
  ) %>%
  as.data.frame()
print(dim(result_df_dplyr))
loci_pair_LRC_bin = result_df_dplyr

tmp = match(result_df_dplyr$posA_bin_tmp,loci_pair_LRC$posA_bin_tmp)
loci_pair_LRC_bin = cbind(loci_pair_LRC[tmp,c("posA_bin","posB_bin","pos","allel", "cell", "cc")],loci_pair_LRC_bin)

LRC_pairs = paste(loci_pair_LRC_bin$posA_bin/1000, loci_pair_LRC_bin$posB_bin/1000)
LRC_pair_top_idx <-which(table(LRC_pairs)%>%as.numeric>10)
print(length(LRC_pair_top_idx))
LRC_pair_top = table(LRC_pairs)[LRC_pair_top_idx]
split_list <- strsplit(names(LRC_pair_top), " ")
result_matrix <- do.call(rbind, split_list)
LRC_pair_top2 <- matrix(as.numeric(result_matrix), ncol = 2)
print(LRC_pair_top)
LRC_pairs = unique(LRC_pairs)
print(length(LRC_pairs))

LRC_site = c(loci_pair_LRC_bin$posA_bin/1000,loci_pair_LRC_bin$posB_bin/1000)
LRC_site_top <- head(sort(table(LRC_site), decreasing = TRUE), 200)
print(LRC_site_top)
LRC_site = unique(LRC_site)%>%sort()
print(length(LRC_site))

### check the epigenetics of pairs of loci when they colocalizaed or not ####
cc = 'chr1'
file_names = paste0('NMuMG_rep1_pos',c(1:10)-1)
processed_data <- list()
for (ii in 1:length(file_names)) {
  fls = paste0(file_names[ii], '.csv')
  print(fls)
  seqFISH_DNA <- read.csv(fls)
  
  # move duplicated loci+pos+cell+allele
  seqFISH_DNA_tmp = cbind(seqFISH_DNA$fov, seqFISH_DNA$cellID, seqFISH_DNA$name, seqFISH_DNA$dbscan_ldp_nbr_allele)
  duplicated_all <- duplicated(seqFISH_DNA_tmp) | duplicated(seqFISH_DNA_tmp, fromLast = TRUE)
  
  if (length(which(duplicated_all)) > 0) {
    seqFISH_DNA = seqFISH_DNA[-which(duplicated_all), ]
  }
  
  split_list <- strsplit(seqFISH_DNA$name, "-")
  loci <- do.call(rbind, split_list)
  loci[, 2] <- as.numeric(loci[, 2])
  
  chr_indices = which(seqFISH_DNA$chr == cc)
  seqFISH_DNA_tmp = seqFISH_DNA[chr_indices, ]
  seqFISH_DNA_tmp$loci = loci[chr_indices, 2]
  
  processed_data[[ii]] = seqFISH_DNA_tmp
}
# combine
seqFISH_DNA <- do.call(rbind, processed_data)
# write.csv(seqFISH_DNA, file = paste0("seqFISH-", cc,".csv"), row.names = FALSE)

### plot epi marks of LRC_pair_top ####
cc = 'chr1'
seqFISH_chr1 = read.csv(paste0("seqFISH-", cc,".csv"))

probes = read.csv('41586_2025_8838_MOESM5_ESM.csv')
probes_ID = strsplit(probes$name, split = '-')
probes_ID = do.call(rbind,probes_ID)
probes = cbind(probes_ID,probes$name,probes$chrom,probes$start)
colnames(probes) = c('chr','probe_ID','names','chrom','start')
probes = data.frame(probes)
probes_tmp = probes[which(probes$chr==cc),]
tmp = match(seqFISH_chr1$name,probes_tmp$names)
seqFISH_chr1_loci = cbind(seqFISH_chr1$name,probes_tmp$names[tmp],as.numeric(probes_tmp$start[tmp])/1000)
seqFISH_chr1$loci = as.numeric(probes_tmp$start[tmp])/1000


loci_bin = seq(0,198000000,250000)
seqFISH_chr1$loci_bin <- loci_bin[findInterval(seqFISH_chr1$loci*1000, loci_bin)]

seqFISH_chr1$loci_bin_tmp = paste0(seqFISH_chr1$loci_bin,"_",seqFISH_chr1$fov,"_",seqFISH_chr1$dbscan_ldp_nbr_allele,"_",seqFISH_chr1$cellID,"_",seqFISH_chr1$chrom)

nselct = c("loci","loci_bin_tmp", epi_marks)
seqFISH_chr1_2 = seqFISH_chr1[,nselct]

result_df_dplyr <- seqFISH_chr1_2 %>%
  group_by(loci_bin_tmp) %>%  
  summarise(
    across(
      .cols = -c(loci),       
      .fns = ~median(.x, na.rm = TRUE),
      .names = "{.col}_md"         
    ),
    .groups = "drop"  
  ) %>%
  as.data.frame()
print(dim(result_df_dplyr))
seqFISH_chr1_bin = result_df_dplyr
tmp = match(seqFISH_chr1_bin$loci_bin_tmp,seqFISH_chr1$loci_bin_tmp)
seqFISH_chr1_bin$loci_bin =seqFISH_chr1$loci_bin[tmp]/1000
seqFISH_chr1_bin$fov = seqFISH_chr1$fov[tmp]
seqFISH_chr1_bin$cellID = seqFISH_chr1$cellID[tmp]
seqFISH_chr1_bin$dbscan_ldp_nbr_allele = seqFISH_chr1$dbscan_ldp_nbr_allele[tmp]

cell_info = paste0('allel', seqFISH_chr1_bin$dbscan_ldp_nbr_allele, ' NMuMG_rep1_pos', seqFISH_chr1_bin$fov, ' ', seqFISH_chr1_bin$cellID)
print(length(unique(cell_info)))

####################################################################################################
# Fig. 4b  ####
####################################################################################################
# epi_marks = paste0(epi_marks,"_md")
# plot_list <- list()
par(mar = c(2,2,1,1))
par(mfrow = c(2,4))
for (ee in epi_marks) {
  tmp = seqFISH_chr1_bin[[paste0(ee,"_md")]]
  abl = quantile(tmp)[4]
  hist(as.numeric(tmp),main = ee, xlim = c(-2,5), col = 'lightblue')
  # abline(v = -1,col = "darkred",lty = 2)
  # abline(v = 1,col = "darkred",lty = 2)
  abline(v = quantile(tmp)[4],col = 'blue',lty = 2)
  box(lwd = 3)
}

####################################################################################################
# Fig. 4c-f  ####
####################################################################################################
for (ll in 1:nrow(LRC_pair_top2)) {
  par(mar = c(2,2,1,1))
  par(mfrow = c(2,4))
  for (ee in epi_marks) {
    print(ee)
    tmp = seqFISH_chr1_bin[[paste0(ee,"_md")]]
    abl = quantile(tmp)[4]
    matched_row_numbers <- which(loci_pair_LRC_bin[, 1]/1000 == as.numeric(LRC_pair_top2[ll, 1]) & loci_pair_LRC_bin[, 2]/1000 == as.numeric(LRC_pair_top2[ll, 2]))
    matched_rows <- loci_pair_LRC_bin[matched_row_numbers, ]
    cell_tmp = paste(matched_rows$allel, matched_rows$pos, matched_rows$cell)
    col = rainbow(nrow(matched_rows)*2)
    
    matched_row_numbers_st <- which((as.numeric(seqFISH_chr1_bin$loci_bin) == LRC_pair_top2[ll, 1]))
    matched_row_numbers_ed <- which((as.numeric(seqFISH_chr1_bin$loci_bin) == LRC_pair_top2[ll, 2]))
    
    matched_st = which(cell_info[matched_row_numbers_st]%in%cell_info[matched_row_numbers_ed])
    seqFISH_chr1_st = seqFISH_chr1_bin[matched_row_numbers_st[matched_st],]
    matched_ed = which(cell_info[matched_row_numbers_ed]%in%cell_info[matched_row_numbers_st])
    seqFISH_chr1_ed = seqFISH_chr1_bin[matched_row_numbers_ed[matched_ed],]
    tmp = match(cell_info[matched_row_numbers_st[matched_st]],cell_info[matched_row_numbers_ed[matched_ed]])
    seqFISH_chr1_ed = seqFISH_chr1_ed[tmp,]
    
    vector1 = seqFISH_chr1_st[[paste0(ee,'_md')]]
    print(head(seqFISH_chr1_st[[paste0(ee,'_md')]]))
    vector2 = seqFISH_chr1_ed[[paste0(ee,'_md')]]
    df = data.frame(
      value = c(vector1,vector2),
      group = c(rep('vector1',length(vector1)),rep('vector2',length(vector2)))
    )
    pt1 = matched_rows[[paste0('loc1.',ee,'_md')]]
    pt2 = matched_rows[[paste0('loc2.',ee,'_md')]]
    
    plot(
      x = vector1,          
      y = vector2,          
      xlim = c(-1.5,4.5),
      ylim = c(-1.5,4.5),
      main = paste0(LRC_pair_top2[ll,1]/1e3, '-',LRC_pair_top2[ll,2]/1e3, 'Mb ',ee),  
      xlab = "epi on locus 2",                
      ylab = "epi on locus 1",                
      pch = 16,             
      col = "grey",      
      cex = 1,              
      alpha = 0.7           
    )
    box(lwd = 3)
    points(pt1, pt2, pch = 19, cex = 2, col = rainbow(nrow(matched_rows)))
    # Add diagonal line (y = x) - core code for diagonal
    abline(a = 0, b = 1, col = "darkgreen", lwd = 2, lty = "dashed")    
    abline(h = abl, col = "blue", lwd = 1, lty =3)
    abline(v = abl, col = "blue", lwd = 1, lty = 3)
    
    # seqFISH_chr1_st = seqFISH_chr1[matched_row_numbers_st,]
    # duplicated_rows <- duplicated(cbind(seqFISH_chr1_st$fov, seqFISH_chr1_st$cellID, seqFISH_chr1_st$dbscan_ldp_nbr_allele))
    # seqFISH_chr1_st <- seqFISH_chr1_st[!duplicated_rows, ]
    # seqFISH_chr1_ed = seqFISH_chr1[matched_row_numbers_ed,]
    # duplicated_rows <- duplicated(cbind(seqFISH_chr1_ed$fov, seqFISH_chr1_ed$cellID, seqFISH_chr1_ed$dbscan_ldp_nbr_allele))
    # seqFISH_chr1_ed <- seqFISH_chr1_ed[!duplicated_rows, ]
    # 
    # vector1 = seqFISH_chr1_st[[ee]]
    # print(head(seqFISH_chr1_st[[ee]]))
    # vector2 = seqFISH_chr1_ed[[ee]]
    # df = data.frame(
    #   value = c(vector1,vector2),
    #   group = c(rep('vector1',length(vector1)),rep('vector2',length(vector2)))
    # )
    # df$group = factor(df$group)
    # ee1 = paste0('loc1.',ee,'_md')
    # ee2 = paste0('loc2.',ee,'_md')
    # x = c(as.numeric(matched_rows[[ee1]]),as.numeric(matched_rows[[ee2]]))
    # points_df <- data.frame(
    #   x = x,
    #   y = rep(0,length(x))
    # )
    # # Add row identifier
    # points_df$row_id <- rep(1:(length(x)/2),each = 2)
    # points_df$row_id = factor(points_df$row_id)
    # # Create a color palette that scales with your points
    # n_points <- length(unique(points_df$row_id))
    # point_colors <- scales::hue_pal()(n_points)  # Automatically generates n distinct colors
    # p = ggplot(df, aes(x = value)) +
    #   geom_density(aes(color = group), linewidth = 0.5) +
    #   geom_point(data = points_df, 
    #              aes(x = x, y = y, color = row_id),
    #              size = 10, shape = 42, alpha = 1, show.legend = FALSE) +
    #   labs(title = paste0(LRC_pair_top2[ll,1]/1e3, '-',LRC_pair_top2[ll,2]/1e3, 'Mb ',ee),
    #        x = "Values",
    #        y = "Density",
    #        color = "Densities",
    #        fill = "Points") +
    #   theme_minimal()+
    #   theme(legend.position = "none")
    # 
    # # plot(p)
    # 
    # plot_list[[length(plot_list) + 1]] <- p
  }
  # combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 2)
  # print(combined_plot)
  # plot_list <- list()
}

####################################################################################################
# G. Fig. S8a  ####
####################################################################################################
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
tt1 = which(lner>=100000000)

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
print(dim(conformation_matrix))
print(head(distances_stat2))
print(head(IMR_2))

bin_thresh = 1000000
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
conformation_matrix22 = conformation_matrix[tmp,]
print(dim(distances_stat22))
print(dim(conformation_matrix22))

# compare with A/B compartment ####
library(rtracklayer)
region <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(start = 1, end = 242000000)
)
# Import the BigWig file
fls = c('4DNFIRL31X9M','4DNFI342UZP1','4DNFIHM89EGL')
bw_data <- "Downloads/4DNFI342UZP1.bw"
fls = paste0("4DNFIHM89EGL.bw")
print(fls)
# Import bigBed data
features <- import(
  con = fls,  # Replace with your bigBed file path
  format = "BigWig",
  which = region
)
# Convert to data frame for easier handling
features_df <- as.data.frame(features)
features_df = rbind(features_df,c("chr2",242000001,242250000,250000,"*",0.79))

distances_stat222 = distances_stat22[which(distances_stat22$end!=242000001),]

####################################################################################################
# A/B score at both ends of each sLRC (mean of the bin and the preceding bin), A-A / B-B / A-B fractions, density comparison all loci vs sLRC loci
####################################################################################################
feature_st = 0
feature_ed = 0
for (ii in 1:nrow(distances_stat222)) {
  tmp = as.numeric(distances_stat222[ii,])
  st_tmp = which((as.numeric(features_df$start)<=tmp[1])&(as.numeric(features_df$end)>tmp[1]))
  ed_tmp = which((as.numeric(features_df$start)<=tmp[2])&(as.numeric(features_df$end)>tmp[2]))
  if ((length(st_tmp)!=0)&(length(ed_tmp)!=0)) {
    feature_st = c(feature_st,mean(as.numeric(features_df$score[max((st_tmp-1),0):st_tmp])))
    feature_ed = c(feature_ed,mean(as.numeric(features_df$score[max((ed_tmp-1),0):ed_tmp])))
  }else{
    feature_st = c(feature_st,NA)
    feature_ed = c(feature_ed,NA)
  }
}
feature_st = feature_st[-1]
feature_ed = feature_ed[-1]
feature_LRC = cbind(distances_stat222,feature_st,feature_ed)
# plot(feature_LRC$feature_st,feature_LRC$feature_ed,pch = '.')
tmp = which((feature_LRC$feature_st>0.3)&(feature_LRC$feature_ed>0.3))
print(length(tmp)/nrow(feature_LRC))
tmp = which((feature_LRC$feature_st<(-0.3))&(feature_LRC$feature_ed<(-0.3)))
print(length(tmp)/nrow(feature_LRC))
tmp = which((feature_LRC$feature_st>(0.3))&(feature_LRC$feature_ed<(-0.3)))
print(length(tmp)/nrow(feature_LRC))
tmp = which((feature_LRC$feature_st<(-0.3))&(feature_LRC$feature_ed>(0.3)))
print(length(tmp)/nrow(feature_LRC))
# Create example data (replace with your actual data)
tmp = unique(feature_LRC$start)
set1 <- na.omit(as.numeric(features_df$score)) 
set2 <- unique(na.omit(c(as.numeric(feature_LRC$feature_st),as.numeric(feature_LRC$feature_ed))))

dens1 <- density(set1)
dens2 <- density(set2)

plot(dens1,
     type = "n",  # Don't plot the line yet
     xlim = c(-1.5,2),
     ylim = c(0,max(max(dens1$y),(max(dens2$y)))),
     main = "Distribution Comparison",
     xlab = "Value",
     ylab = "Density")

# Add filled areas with transparency
# For set1 (blue)
polygon(c(dens1$x, rev(dens1$x)),
        c(dens1$y, rep(0, length(dens1$y))),
        col = rgb(0, 0, 1, 0.5),  # Blue with 50% transparency
        border = NA)
# For set2 (red)
polygon(c(dens2$x, rev(dens2$x)),
        c(dens2$y, rep(0, length(dens2$y))),
        col = rgb(0, 0.8, 0, 0.5),  # Red with 50% transparency
        border = NA)
lines(dens1, col = "blue", lwd = 2)
lines(dens2, col = "green", lwd = 2)

####################################################################################################
# S8b-h ####
####################################################################################################
# 3D chromosome contact =====
pairs_files <- list.files("mesc_pairs/pairs",
                          pattern = "\\.pairs\\.gz$",
                          full.names = TRUE)
print(pairs_files)
print(length(pairs_files))
cat("Loading 3D Genome data...\n")
#bin and windws
bin_size <- 250000
window_size <- 100000
process_single_file <- function(file_path) {
  
  # 1. Rapid file reading (select only the required columns during the read process to save memory)
  genome_data <- fread(file_path, select = c(2, 3, 4, 5, 8, 9), 
                       col.names = c("chr1", "pos1", "chr2", "pos2", "phase0", "phase1"))
  
  # 2. Vectorized operations + one-pass filtering (eliminating cumbersome `which()` and `rmv` calls)
  processed_data <- genome_data %>%
    mutate(
      # calculate center
      anchor1 = round(pos1 / bin_size) * bin_size,
      anchor2 = round(pos2 / bin_size) * bin_size,
      
      # windows
      bin1 = ifelse(abs(pos1 - anchor1) <= window_size, anchor1, NA),
      bin2 = ifelse(abs(pos2 - anchor2) <= window_size, anchor2, NA),
      
      # add column cellID (file_path / cl)
      cellID = file_path
    ) %>%
    # [Core Filtering Step] Performs all filtering in a single pass—cleaner code, faster execution.
    filter(
      !is.na(bin1),                                # Filter out rows where bin1 is NA.
      !is.na(bin2),                                # Filter out rows where bin2 is NA.
      chr1 == chr2,                                # have to be intra-chromosomal
      abs(bin2 - bin1) > 100000000                 # >100000000
    ) %>%
    # delete temperol anchor
    select(-anchor1, -anchor2)
  
  # 3.Split the processed data by allele and return it packaged as a list.
  list(
    phase0 = processed_data %>% filter(phase0 == 1),
    phase1 = processed_data %>% filter(phase1 == 1)
  )
}

# Batch processing and merging (avoiding for-loops with rbind) ####
cat("Starting batch processing of all files...\n")

# lapply iterates through all files in pairs_files and stores the result (a list) for each file.
# This process is significantly faster than using a for-loop.
all_results <- lapply(pairs_files, process_single_file)

cat("File processing complete; performing final merge...\n")
# Use dplyr's bind_rows to merge all results at once (requiring only a single memory allocation)!
# Extract the "phase0" component from all results and merge them
genome_data_binned1 <- bind_rows(lapply(all_results, `[[`, "phase0"))
# Extract the "phase1" component from all results and merge them
genome_data_binned2 <- bind_rows(lapply(all_results, `[[`, "phase1"))

# change the column order if wanted
col_order <- c("cellID", "chr1", "pos1", "chr2", "pos2", "phase0", "phase1", "bin1", "bin2")
genome_data_binned1 <- genome_data_binned1 %>% select(all_of(col_order))
genome_data_binned2 <- genome_data_binned2 %>% select(all_of(col_order))

# check dimesion
cat("phase0 (C57BL/6J) dimesion:\n")
print(dim(genome_data_binned1))

cat("phase1 (CAST/EiJ) dimesion:\n")
print(dim(genome_data_binned2))

# drop the placeholder first row (all zeros) used at initialization
genome_data_binned1 <- as.data.table(genome_data_binned1)[-1, ]
genome_data_binned2 <- as.data.table(genome_data_binned2)[-1, ]

fwrite(genome_data_binned1, "mesc.genome_data_binned1.csv")
fwrite(genome_data_binned2, "mesc.genome_data_binned2.csv")

## input ####
# The genes located at or near those loci
genome_data_binned1 <- fread("mesc.genome_data_binned1.csv")
genome_data_binned2 <- fread("mesc.genome_data_binned2.csv")
print(dim(genome_data_binned1))
print(head(genome_data_binned1))
print(dim(genome_data_binned2))
print(head(genome_data_binned2))
# extract a clean cell ID from the file path that matches RNA matrix
# column names (e.g. ".../pairs/P10001.something.pairs.gz" -> "P10001")
# *** CHECK THIS REGEX against your actual file names / RNA colnames ***
extract_cell_id <- function(path) {
  sub(".*?(P[0-9]+).*", "\\1", path)
}
genome_data_binned1[, cell := extract_cell_id(cellID)]
genome_data_binned2[, cell := extract_cell_id(cellID)]
print(dim(genome_data_binned1))
print(head(genome_data_binned1))
print(dim(genome_data_binned2))
print(head(genome_data_binned2))

# sanity check: do extracted IDs overlap RNA matrix column names?
cat("Example extracted cell IDs (phase0):\n")
print(head(unique(genome_data_binned1$cell)))
cat("Example extracted cell IDs (phase1):\n")
print(head(unique(genome_data_binned2$cell)))

# File paths (ensure these are in your working directory)
rna_file <- "GSM9115533_mesc.RNAcounts.gene.total.format.tsv.gz"
allel1 <- "GSM9115533_mesc.RNAcounts.gene.b6.format.tsv.gz"
allel2 <- "GSM9115533_mesc.RNAcounts.gene.cast.format.tsv.gz"

cat("\nLoading RNA expression data...\n")
RNA_B6 <- fread(allel1)
print(dim(RNA_B6))
print(RNA_B6[1:5,1:5])
RNA_cast <- fread(allel2)
print(dim(RNA_cast))
print(RNA_cast[1:5,1:5])
cat("Overlap (B7):", length(intersect(unique(genome_data_binned1$cell),
                                      colnames(RNA_B6))), "\n")
cat("Overlap (CAST):", length(intersect(unique(genome_data_binned2$cell),
                                        colnames(RNA_cast))), "\n")
target_cell = unique(c(colnames(RNA_B6),colnames(RNA_cast)))[-1]
print(length(target_cell))
print(head(target_cell))

## Gene activity at the ends of distal (>100Mb) ==================

## Gene -> 250kb bin annotation  ---------------
# 
# # to Ensembl mice data
# ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# attributes <- listAttributes(ensembl)
# head(attributes)
# 
# # gene ID, chr, loci, gene name...
# gene_annotation <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name", 
#                  "chromosome_name", "start_position", "end_position",
#                  "gene_biotype"),
#   mart = ensembl
# )
# head(gene_annotation)
# gene_annotation <- gene_annotation %>%
#   filter(chromosome_name %in% c(as.character(1:19), "X", "Y"))
# #save to local
# saveRDS(gene_annotation, "mouse_gene_annotation.rds")

gene_annotation <- readRDS("XingD/mouse_gene_annotation.rds")

gene_tab <- as.data.table(gene_annotation)[
  , .(gene = external_gene_name,
      chr  = chromosome_name,
      start = start_position)
]

# drop genes with no symbol, or duplicated symbols (e.g. different
# loci sharing a name) -- keep the first occurrence per symbol+chr
gene_tab <- gene_tab[gene != "" & !is.na(start)]
gene_tab <- unique(gene_tab, by = c("gene", "chr"))

# same binning rule as the contacts
gene_tab[, anchor := round(start / bin_size) * bin_size]
gene_tab[, bin    := ifelse(abs(start - anchor) <= window_size, anchor, NA_real_)]
gene_tab <- gene_tab[!is.na(bin)]
gene_tab[, anchor := NULL]

# biomaRt chromosome_name is "1","2",...,"X","Y" (no "chr" prefix).
# *** CHECK: do genome_data_binned$chr1 values look like "chr1" or "1"? ***
cat("Example chr values in contacts:", head(unique(genome_data_binned1$chr1)), "\n")
cat("Example chr values in contacts:", head(unique(genome_data_binned2$chr1)), "\n")
cat("Example chr values in gene_tab:", head(unique(gene_tab$chr)), "\n")
# uncomment ONE of these if formats don't match:
gene_tab[, chr := paste0("chr", chr)]                    # "1" -> "chr1"
# genome_data_binned1[, chr1 := sub("^chr", "", chr1)]     # "chr1" -> "1"
# genome_data_binned2[, chr1 := sub("^chr", "", chr1)]

cat("Begin converting the RNA matrix format and matching bin coordinates....\n")
# 1. Wide-to-long format conversion (highly efficient melt)
# the first line of RNA_B6 is gene，then P10001, P10002...
rna_long_B6 <- melt(RNA_B6, id.vars = "gene", variable.name = "cell", value.name = "expr")
rna_long_cast <- melt(RNA_cast, id.vars = "gene", variable.name = "cell", value.name = "expr")

# 2. Perform an inner join (merge) with your gene_tab.
# This step automatically discards genes falling within the "gray zone (NA)" we have defined.
rna_binned_B6 <- merge(rna_long_B6, gene_tab, by = "gene", all.x = FALSE)
rna_binned_cast <- merge(rna_long_cast, gene_tab, by = "gene", all.x = FALSE)
setDT(rna_binned_B6)
setDT(rna_binned_cast)
rna_binned_B6[,Allele:='B6']
rna_binned_cast[,Allele:='cast']
rna_binned = rbind(rna_binned_B6, rna_binned_cast)
# 3. Aggregate
# Core logic: Group by chromosome (chr), anchor (bin), and cell, and sum the gene expression levels for all genes falling within the same 250kb range.
rna_binned_agg_B6 <- as.data.table(rna_binned_B6)[, .(bin_expr = sum(expr)), by = .(chr, bin, cell)]
rna_binned_agg_cast <- as.data.table(rna_binned_cast)[, .(bin_expr = sum(expr)), by = .(chr, bin, cell)]

# 4. create index (Bin_ID): "chr1_3500000"
rna_binned_agg_B6[, bin_id := paste(chr, as.character(bin), sep = "_")]
rna_binned_agg_cast[, bin_id := paste(chr, as.character(bin), sep = "_")]

# 5. back to the wide form, add 0 for NA (dcast)
rna_bin_matrix_B6 <- dcast(rna_binned_agg_B6, bin_id ~ cell, value.var = "bin_expr", fill = 0)
rna_bin_matrix_cast <- dcast(rna_binned_agg_cast, bin_id ~ cell, value.var = "bin_expr", fill = 0)

cat("\nRNA expression with bins \n")
print(dim(rna_bin_matrix_B6))
print(dim(rna_bin_matrix_cast))
print(rna_bin_matrix_B6[1:5, 1:5])
print(rna_bin_matrix_cast[1:5, 1:5])

# 1. remove the first column (bin id) and revert to matrix 
expr_matrix_B6 <- as.matrix(rna_bin_matrix_B6[,-1])
expr_matrix_cast <- as.matrix(rna_bin_matrix_cast[,-1])
# 2. row names = bin_id
rownames(expr_matrix_B6) <- rna_bin_matrix_B6$bin_id
rownames(expr_matrix_cast) <- rna_bin_matrix_cast$bin_id
# 3. double check
print(dim(expr_matrix_B6))
print(expr_matrix_B6[1:5, 1:5])

cat("\n=== Step 2: Extract and filter global 3D interaction pairs ===\n")

# statistic % >100Mb contact
global_pairs_phase0 <- genome_data_binned1 %>%
  filter(chr1 == chr2, abs(bin2 - bin1) > 100000000) %>%
  group_by(chr1, bin1, bin2) %>%
  summarize(
    contact_cells = list(unique(cell)),
    n_cells = length(unique(cell)),
    .groups = "drop"
  )
global_pairs_phase1 <- genome_data_binned2 %>%
  filter(chr1 == chr2, abs(bin2 - bin1) > 100000000) %>%
  group_by(chr1, bin1, bin2) %>%
  summarize(
    contact_cells = list(unique(cell)),
    n_cells = length(unique(cell)),
    .groups = "drop"
  )
#combine
global_pairs_phase0$Allele <- "B6"
global_pairs_phase1$Allele <- "CAST"
global_pairs <- rbind(global_pairs_phase0, global_pairs_phase1)

# filter out too rare contact
# cell number is 960， Contact cells >= 20，No-Contact cells >= 20
total_cells_B6 <- colnames(expr_matrix_B6)
total_cells_cast <- colnames(expr_matrix_cast)

min_cells_required <- 20

valid_pairs_phase0 <- global_pairs_phase0 %>%
  filter(n_cells >= min_cells_required & (length(total_cells_B6) - n_cells) >= min_cells_required)
valid_pairs_phase1 <- global_pairs_phase1 %>%
  filter(n_cells >= min_cells_required & (length(total_cells_cast) - n_cells) >= min_cells_required)
print(dim(valid_pairs_phase0))
print(dim(valid_pairs_phase1))
print(dim(global_pairs_phase0))
print(dim(global_pairs_phase1))

valid_pairs_phase0$bin1_id <- paste(valid_pairs_phase0$chr1, valid_pairs_phase0$bin1, sep = "_")
valid_pairs_phase0$bin2_id <- paste(valid_pairs_phase0$chr1, valid_pairs_phase0$bin2, sep = "_")
valid_pairs_phase1$bin1_id <- paste(valid_pairs_phase1$chr1, valid_pairs_phase1$bin1, sep = "_")
valid_pairs_phase1$bin2_id <- paste(valid_pairs_phase1$chr1, valid_pairs_phase1$bin2, sep = "_")
valid_pairss = rbind(valid_pairs_phase0,valid_pairs_phase1)

valid_pairs <- valid_pairss %>%
  filter(bin1_id %in% c(rownames(expr_matrix_B6),rownames(expr_matrix_cast)) & bin2_id %in%c(rownames(expr_matrix_B6),rownames(expr_matrix_cast)) )
cat(sprintf("after filtering，there are %d in the genome fit the QC for long-range contact \n", nrow(valid_pairs)))

####################################################################################################
# Fig. S8e-h  gene expression at the two anchors of each tested pair: contact vs no-contact alleles
####################################################################################################

# 1. Re-label alleles by combining cell and Allele
rna_binned[, allele_id := paste(cell, Allele, sep = "_")]
valid_pairs[, allele_id := lapply(contact_cells, function(x) paste(x, Allele, sep = "_"))]

# 2. Get gene-bin mapping
gene_bins <- unique(rna_binned[, .(gene, chr, bin, allele_id)])
gene_bins[, active := FALSE]
gene_bins$active[which(rna_binned$expr>0)]=TRUE

# 3. Create complete gene-allele matrix
all_genes <- unique(rna_binned$gene)
all_alleles <- unique(rna_binned$allele_id)
all_gene_allele <- CJ(gene = all_genes, allele_id = all_alleles)

all_gene_allele <- merge(all_gene_allele, 
                         gene_bins[, .(gene, allele_id, active)], 
                         by = c("gene", "allele_id"), 
                         all.x = TRUE)
all_gene_allele[is.na(active), active := FALSE]

# 4. Get genes for each bin
bin_to_genes <- gene_bins[, .(genes = list(unique(gene))), by = .(chr, bin)]

# 5. Expand valid_pairs to pair-gene-allele level
valid_pairs_temp <- copy(valid_pairs)
valid_pairs_temp[, pair_id := .I]

# Get genes for bin1
valid_pairs_bin1 <- merge(valid_pairs_temp, bin_to_genes, 
                          by.x = c("chr1", "bin1"), 
                          by.y = c("chr", "bin"), 
                          all.x = TRUE)
setnames(valid_pairs_bin1, "genes", "genes_bin1")

# Get genes for bin2
valid_pairs_bin2 <- merge(valid_pairs_temp, bin_to_genes, 
                          by.x = c("chr1", "bin2"), 
                          by.y = c("chr", "bin"), 
                          all.x = TRUE)
setnames(valid_pairs_bin2, "genes", "genes_bin2")

# Combine
valid_pairs_genes <- merge(valid_pairs_bin1[, .(pair_id, chr1, bin1, bin2, n_cells, Allele, contact_cells, genes_bin1)],
                           valid_pairs_bin2[, .(pair_id, genes_bin2)],
                           by = "pair_id")

valid_pairs_genes[, genes_in_pair := mapply(c, genes_bin1, genes_bin2, SIMPLIFY = FALSE)]
valid_pairs_genes[, genes_in_pair := lapply(genes_in_pair, unique)]

# Expand to gene level
# Use a loop or apply to avoid list column issues
expanded_list <- list()
for(i in 1:nrow(valid_pairs_genes)) {
  genes_vec <- valid_pairs_genes$genes_in_pair[[i]]
  contact_alleles_vec <- unlist(lapply(valid_pairs_genes$contact_cells[[i]], 
                                       function(x) paste(x, valid_pairs_genes$Allele[i], sep = "_")))
  
  for(g in genes_vec) {
    expanded_list[[length(expanded_list) + 1]] <- data.table(
      pair_id = valid_pairs_genes$pair_id[i],
      chr1 = valid_pairs_genes$chr1[i],
      bin1 = valid_pairs_genes$bin1[i],
      bin2 = valid_pairs_genes$bin2[i],
      n_cells = valid_pairs_genes$n_cells[i],
      Allele = valid_pairs_genes$Allele[i],
      gene = g,
      contact_alleles = list(contact_alleles_vec)
    )
  }
}

valid_pairs_expanded <- rbindlist(expanded_list)

# Remove duplicates (using gene, pair_id, and contact_alleles as keys)
valid_pairs_expanded <- unique(valid_pairs_expanded, by = c("pair_id", "gene"))

cat("Number of gene-pair combinations:", nrow(valid_pairs_expanded), "\n")

# 6. For each gene-pair combination, calculate statistics
gene_pair_stats_list <- list()
counter <- 1

for(i in 1:nrow(valid_pairs_expanded)) {
  print(i)
  gene_i <- valid_pairs_expanded$gene[i]
  pair_id_i <- valid_pairs_expanded$pair_id[i]
  contact_alleles_vec <- unlist(valid_pairs_expanded$contact_alleles[i])
  
  # Get all alleles for this gene
  gene_allele_data <- all_gene_allele[gene == gene_i]
  
  # Contact alleles for this gene in this specific pair
  contact_data <- gene_allele_data[allele_id %in% contact_alleles_vec]
  contact_total <- nrow(contact_data)
  contact_active <- sum(contact_data$active)
  contact_prop <- if(contact_total > 0) contact_active / contact_total else 0
  
  # No-contact alleles for this gene
  nocontact_data <- gene_allele_data[!allele_id %in% contact_alleles_vec]
  nocontact_total <- nrow(nocontact_data)
  nocontact_active <- sum(nocontact_data$active)
  nocontact_prop <- if(nocontact_total > 0) nocontact_active / nocontact_total else 0
  
  # Mean expression in active contact alleles
  contact_active_expr <- rna_binned[expr > 0 & gene == gene_i & allele_id %in% contact_alleles_vec, 
                                    mean(expr, na.rm = TRUE)]
  if(length(contact_active_expr) == 0) contact_active_expr <- 0
  
  # Mean expression in active no-contact alleles
  nocontact_active_expr <- rna_binned[expr > 0 & gene == gene_i & !allele_id %in% contact_alleles_vec, 
                                      mean(expr, na.rm = TRUE)]
  if(length(nocontact_active_expr) == 0) nocontact_active_expr <- 0
  
  gene_pair_stats_list[[counter]] <- data.table(
    gene = gene_i,
    pair_id = pair_id_i,
    chr1 = valid_pairs_expanded$chr1[i],
    bin1 = valid_pairs_expanded$bin1[i],
    bin2 = valid_pairs_expanded$bin2[i],
    contact_total = contact_total,
    contact_active = contact_active,
    contact_prop = contact_prop,
    nocontact_total = nocontact_total,
    nocontact_active = nocontact_active,
    nocontact_prop = nocontact_prop,
    contact_mean_expr = contact_active_expr,
    nocontact_mean_expr = nocontact_active_expr
  )
  
  counter <- counter + 1
}

gene_pair_stats <- rbindlist(gene_pair_stats_list)

# 7. Filter out rows with NA or insufficient data
gene_pair_stats_filtered <- gene_pair_stats[
  !is.na(contact_prop) & 
    !is.na(nocontact_prop) & 
    contact_total >= 3 & 
    nocontact_total >= 3
]

cat("Total gene-pair combinations:", nrow(gene_pair_stats), "\n")
cat("Filtered gene-pair combinations:", nrow(gene_pair_stats_filtered), "\n")
cat("Unique genes:", length(unique(gene_pair_stats_filtered$gene)), "\n")
cat("Unique pairs:", length(unique(gene_pair_stats_filtered$pair_id)), "\n")

# 8. Create plots
if(nrow(gene_pair_stats_filtered) > 0) {
  plot_data <- as.data.frame(gene_pair_stats_filtered)
  
  # Plot 1: Active proportion
  p1 <- ggplot(plot_data, aes(x = contact_prop, y = nocontact_prop)) +
    geom_point(size = 2, alpha = 0.5, color = "#2C3E50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    labs(
      title = "Gene Active Proportion: Contact vs No-Contact Alleles",
      subtitle = paste("Each point = one gene in one contact pair (n =", nrow(plot_data), ")"),
      x = "Active Proportion in Contact Alleles",
      y = "Active Proportion in No-Contact Alleles"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    ) +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1)
  
  print(p1)
  
  # Plot 2: Mean expression
  p2 <- ggplot(plot_data, aes(x = contact_mean_expr, y = nocontact_mean_expr)) +
    geom_point(size = 2, alpha = 0.5, color = "#2C3E50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    labs(
      title = "Mean Expression in Active Alleles: Contact vs No-Contact",
      subtitle = paste("Each point = one gene in one contact pair (n =", nrow(plot_data), ")"),
      x = "Mean Expression in Contact Alleles",
      y = "Mean Expression in No-Contact Alleles"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    ) +
    coord_fixed()
  
  print(p2)
  
  # Save results
  # write.csv(gene_pair_stats, "gene_pair_stats_all.csv", row.names = FALSE)
  # write.csv(gene_pair_stats_filtered, "gene_pair_stats_filtered.csv", row.names = FALSE)
  # 
  # ggsave("active_proportion_per_pair.pdf", p1, width = 8, height = 7)
  # ggsave("mean_expression_per_pair.pdf", p2, width = 8, height = 7)
  
  # Show distribution
  pairs_per_gene <- gene_pair_stats_filtered[, .N, by = gene]
  cat("\nPairs per gene summary:\n")
  print(summary(pairs_per_gene$N))
  
} else {
  cat("No gene-pair combinations passed the filtering criteria.\n")
  cat("Try lowering the threshold (contact_total >= 1 & nocontact_total >= 1)\n")
}


library(data.table)
library(ggplot2)

# Get all unique pairs
unique_pairs <- unique(plot_data$pair_id)

# Create directory for plots
dir.create("pair_plots", showWarnings = FALSE)

# Calculate average active proportion for each point
setDT(plot_data)
plot_data[, avg_prop := (contact_prop + nocontact_prop) / 2]
plot_data$contact_mean_expr[which(is.na(plot_data$contact_mean_expr))]=0
plot_data$nocontact_mean_expr[which(is.na(plot_data$nocontact_mean_expr))]=0
plot_data[, avg_epxr := (contact_mean_expr + nocontact_mean_expr) / 2]
plot_data$contact_mean_expr = (plot_data$contact_mean_expr*plot_data$contact_active)/plot_data$contact_total
plot_data$nocontact_mean_expr = (plot_data$nocontact_mean_expr*plot_data$nocontact_active)/plot_data$nocontact_total
# Create combined plot
p_combined <- ggplot(plot_data, aes(x = contact_mean_expr, y = nocontact_mean_expr)) +
  geom_point(aes(size = avg_prop, color = factor(pair_id)), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  # scale_size_continuous(
  #   name = "Avg Active Proportion",
  #   range = c(1, 8),
  #   limits = c(0, 1)
  # ) +
  labs(
    title = "Gene Expression: Contact vs No-Contact Alleles",
    subtitle = paste("Each point = one gene in one contact pair (n =", nrow(plot_data), ")"),
    x = "Mean Expression in Contact Alleles",
    y = "Mean Expression in No-Contact Alleles",
    color = "Pair ID"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  ) +
  coord_fixed()

print(p_combined)

p_combined <- ggplot(plot_data, aes(x = contact_prop, y = nocontact_prop)) +
  geom_point(aes(size = avg_epxr, color = factor(pair_id)), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  # scale_size_continuous(
  #   name = "Avg Active Proportion",
  #   range = c(1, 8),
  #   limits = c(0, 1)
  # ) +
  labs(
    title = "Active gene proprotion: Contact vs No-Contact Alleles",
    subtitle = paste("Each point = one gene in one contact pair (n =", nrow(plot_data), ")"),
    x = "Active gene proprotion in Contact Alleles",
    y = "Active gene proprotion in No-Contact Alleles",
    color = "Pair ID"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  ) +
  coord_fixed()

print(p_combined)

# # Save
# ggsave("combined_expression_bubble.pdf", p_combined, width = 12, height = 8)

# Loop through each pair and create plots
unique_pairs = c(52, 57, 19) # examples
for(pid in unique_pairs) {
  
  # Subset data for this pair
  pair_data <- plot_data[which(plot_data$pair_id == pid),]
  pair_data$contact_mean_expr[which(is.na(pair_data$contact_mean_expr))]=0
  pair_data$nocontact_mean_expr[which(is.na(pair_data$nocontact_mean_expr))]=0
  # Skip if too few genes
  if(nrow(pair_data) < 3) next
  
  # Get pair info for title
  title_text <- paste0("Pair ", pid, ": ", 
                       unique(pair_data$chr1)[1], "_", 
                       format(unique(pair_data$bin1)[1], scientific = FALSE), " & ",
                       format(unique(pair_data$bin2)[1], scientific = FALSE))
  
  # Plot 1: Active proportion
  p1 <- ggplot(pair_data, aes(x = contact_prop, y = nocontact_prop)) +
    geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    labs(
      title = title_text,
      subtitle = paste("n =", nrow(pair_data), "genes"),
      x = "Active Proportion in Contact Alleles",
      y = "Active Proportion in No-Contact Alleles"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1)
  
  # Plot 2: Mean expression
  p2 <- ggplot(pair_data, aes(x = contact_mean_expr, y = nocontact_mean_expr)) +
    geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    labs(
      title = title_text,
      subtitle = paste("n =", nrow(pair_data), "genes"),
      x = "Mean Expression in Contact Alleles",
      y = "Mean Expression in No-Contact Alleles"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    xlim(0, 2) +
    ylim(0, 2) +
    coord_fixed()
  # print(p1)
  print(p2)
  # Save
  # ggsave(paste0("pair_plots/pair_", pid, "_proportion.pdf"), p1, width = 8, height = 7)
  # ggsave(paste0("pair_plots/pair_", pid, "_expression.pdf"), p2, width = 8, height = 7)
  
  cat("Saved pair", pid, "\n")
}
#52， 57，19
cat("Done! Check 'pair_plots' folder.\n")

####################################################################################################
# Fig. S8c, S8d  cell-cycle phase (S / G2M marker z-scores), contact penetrance per phase, chi-square test (simulated p-value) ####
####################################################################################################
# Cell cycle #####
cat("\n=== Step 1: RNA Inference for Cell Cycle Scoring ===\n")

# 1. Define highly conserved mouse cell cycle marker genes (Canonical)
s_genes <- c("Mcm2", "Mcm4", "Mcm5", "Mcm6", "Pcna", "Slbp", "Gins2", "Ung", "Fen1")
g2m_genes <- c("Top2a", "Cdk1", "Ccnb1", "Ccnb2", "Aurka", "Birc5", "Mki67", "Tpx2", "Ndc80")

all_markers <- c(s_genes, g2m_genes)

# 2. Extract marker genes from both alleles (B6 and CAST)
# Cell cycle is a global cellular state, so we must sum the alleles
b6_sub <- RNA_B6[gene %in% all_markers]
cast_sub <- RNA_cast[gene %in% all_markers]
print(dim(b6_sub))
print(dim(cast_sub))
# Ensure the genes are in the exact same order for matrix addition
setkey(b6_sub, gene)
setkey(cast_sub, gene)
genes_kept <- intersect(b6_sub$gene, cast_sub$gene)

# Convert to pure numeric matrices
mat_b6 <- as.matrix(b6_sub[, -1, with=FALSE])
mat_cast <- as.matrix(cast_sub[, -1, with=FALSE])

# Total un-phased expression for marker genes
# check diff(match(colnames(mat_b6), colnames(mat_cast)))
mat_total <- mat_b6 + mat_cast
rownames(mat_total) <- genes_kept

# 3. Simple Normalization (Counts Per 10k to account for sequencing depth differences)
# Calculate total UMI across the entire transcriptome per cell
total_umi <- colSums(as.matrix(RNA_B6[, -1, with=FALSE])) + colSums(as.matrix(RNA_cast[, -1, with=FALSE]))
mat_norm <- t(t(mat_total) / total_umi) * 10000
mat_log <- log1p(mat_norm) # log(x+1) transformation

# 4. Calculate Phase Scores
s_kept <- intersect(s_genes, genes_kept)
g2m_kept <- intersect(g2m_genes, genes_kept)

s_score <- colMeans(mat_log[s_kept, , drop=FALSE])
g2m_score <- colMeans(mat_log[g2m_kept, , drop=FALSE])

# 5. Assign Cell Cycle Phase based on Z-scores
cc_df <- data.frame(
  cell = colnames(mat_log),
  S_Score = s_score,
  G2M_Score = g2m_score
)

# Standardize scores (Z-score) to find relative activation
cc_df$S_z <- scale(cc_df$S_Score)[, 1]
cc_df$G2M_z <- scale(cc_df$G2M_Score)[, 1]

# Default everything to G1
cc_df$Phase <- "G1"
# If S is highly active and dominates G2M, it's S phase
cc_df$Phase[cc_df$S_z > 0 & cc_df$S_z > cc_df$G2M_z] <- "S"
# If G2M is highly active and dominates S, it's G2/M phase
cc_df$Phase[cc_df$G2M_z > 0 & cc_df$G2M_z > cc_df$S_z] <- "G2M"

# Lock the factor order for consistent plotting
cc_df$Phase <- factor(cc_df$Phase, levels = c("G1", "S", "G2M"))

cat("Cell Cycle Assignment Summary:\n")
print(table(cc_df$Phase))


cat("\n=== Step 2 & 3: Penetrance Calculation and Statistical Testing ===\n")

# Define a function to test cell-cycle dependency for a single 3D interaction pair
test_cell_cycle <- function(i) {
  row_data <- valid_pairs[i, ]
  # Retrieve cells that have this contact
  contact_cells <- intersect(unlist(row_data$contact_cells), cc_df$cell)
  
  # Mark cell status in the metadata
  temp_df <- cc_df
  temp_df$Contact_Status <- ifelse(temp_df$cell %in% contact_cells, "Contact", "No_Contact")
  
  # Create Contingency Table (Rows: Phase, Cols: Status)
  tbl <- table(temp_df$Phase, temp_df$Contact_Status)
  
  # Ensure the table has both columns before testing
  if("Contact" %in% colnames(tbl) && "No_Contact" %in% colnames(tbl)) {
    # Perform Chi-square test with simulated p-value (robust for small counts)
    test_res <- chisq.test(tbl, simulate.p.value = TRUE)
    
    # Calculate Penetrance (Probability of Contact given the Phase)
    # Formula: Contact / (Contact + No_Contact)
    g1_pen  <- tbl["G1", "Contact"] / sum(tbl["G1", ])
    s_pen   <- tbl["S", "Contact"] / sum(tbl["S", ])
    g2m_pen <- tbl["G2M", "Contact"] / sum(tbl["G2M", ])
    
    return(data.frame(
      pair_id = paste(row_data$bin1_id, row_data$bin2_id, sep="|"),
      Allele = row_data$Allele,
      P_Value = test_res$p.value,
      G1 = g1_pen,
      S = s_pen,
      G2M = g2m_pen
    ))
  }
  return(NULL)
}

# Run the test across all valid pairs
cc_results_list <- pblapply(1:nrow(valid_pairs), test_cell_cycle)
cc_results_df <- bind_rows(cc_results_list)

# Filter for statistically significant pairs (P < 0.05)
sig_cc_pairs <- cc_results_df %>% filter(P_Value < 0.05)

cat(sprintf("\nFound %d >100Mb pairs showing SIGNIFICANT cell-cycle dependency (P < 0.05).\n", nrow(sig_cc_pairs)))


cat("\n=== Step 4: High-Dimensional Visualization (Heatmap) ===\n")

# If we found significant pairs, plot them!
if(nrow(sig_cc_pairs) > 0) {
  
  # 1. Melt the data for ggplot heatmap
  hm_data <- melt(sig_cc_pairs[, c("pair_id", "G1", "S", "G2M")],
                  id.vars = "pair_id", variable.name = "Phase", value.name = "Penetrance")
  
  # 2. Row-wise Z-score normalization (Highlights the PREFERRED phase for each pair)
  hm_data <- hm_data %>%
    group_by(pair_id) %>%
    mutate(Z_Score = scale(Penetrance)[, 1]) %>%
    ungroup()
  
  # 3. Create a clean display name for the Y-axis
  hm_data$Display_Name <- sub("chr1_", "", hm_data$pair_id)
  
  # 4. Order the pairs computationally (clustering) to make the heatmap look professional
  # We sort them by which phase has the highest Z-score
  order_df <- hm_data %>%
    group_by(pair_id) %>%
    summarize(Max_Phase = Phase[which.max(Z_Score)], .groups = "drop") %>%
    arrange(Max_Phase)
  
  hm_data$pair_id <- factor(hm_data$pair_id, levels = order_df$pair_id)
  
  # 5. Plot the Z-score Heatmap
  p_heatmap <- ggplot(hm_data, aes(x = Phase, y = pair_id, fill = Z_Score)) +
    geom_tile(color = "white", size = 0.5) +
    # Diverging color scale: Blue (Depleted) -> White (Neutral) -> Red (Enriched)
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0, 
                         name = "Relative\nEnrichment\n(Z-Score)") +
    labs(
      title = "Cell-Cycle Specific Long-Range Chromatin Contacts",
      subtitle = "Row-normalized contact probability (P < 0.05)",
      x = "Cell Cycle Phase",
      y = "Significant >100Mb Loci Pairs"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(), # Hide long messy bin IDs
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p_heatmap)
  
} else {
  cat("\nNo pairs showed significant cell cycle bias. This implies these long-range contacts are stable across the cell cycle!\n")
}

# Plot ALL pairs sorted by p-value (significant ones highlighted)

cat("\n=== Step 4: High-Dimensional Visualization (Heatmap) ===\n")

# 1. Melt the data for ALL pairs (not just significant ones)
hm_data_all <- melt(cc_results_df[, c("pair_id", "G1", "S", "G2M", "P_Value")],
                    id.vars = c("pair_id", "P_Value"), 
                    variable.name = "Phase", 
                    value.name = "Penetrance")

# 2. Row-wise Z-score normalization
hm_data_all <- hm_data_all %>%
  group_by(pair_id) %>%
  mutate(Z_Score = scale(Penetrance)[, 1]) %>%
  ungroup()

# 3. Sort pairs by P-value (ascending) for plotting order
# Get unique pair_ids sorted by P-value
pair_order <- cc_results_df %>%
  arrange(P_Value) %>%
  distinct(pair_id) %>%  # Add this to remove duplicates
  pull(pair_id)

# Now set the factor levels
hm_data_all$pair_id <- factor(hm_data_all$pair_id, levels = pair_order)

# 4. Create significance annotation column
hm_data_all$Significant <- ifelse(hm_data_all$P_Value < 0.05, "Yes", "No")

# 5. Plot ALL pairs with significance annotation
p_heatmap_all <- ggplot(hm_data_all, aes(x = Phase, y = pair_id, fill = Z_Score)) +
  geom_tile(color = "white", size = 0.3) +
  # Add significance markers on the left side
  geom_tile(data = hm_data_all %>% 
              group_by(pair_id) %>% 
              slice(1) %>% 
              mutate(Phase = "G1"),  # Place marker at left of G1 column
            aes(x = -0.5, fill = NULL, color = Significant), 
            fill = NA, size = 1, width = 0.2) +
  scale_color_manual(values = c("Yes" = "red", "No" = "gray80"), 
                     name = "P < 0.05") +
  # Diverging color scale
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0, 
                       name = "Relative\nEnrichment\n(Z-Score)") +
  labs(
    title = "Cell-Cycle Specific Long-Range Chromatin Contacts",
    subtitle = sprintf("All %d pairs sorted by P-value (Red marker = P < 0.05)", 
                       nrow(cc_results_df)),
    x = "Cell Cycle Phase",
    y = "Loci Pairs (sorted by P-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_heatmap_all)
# Option 2: Alternative - Facet by significance
p_heatmap_facet <- ggplot(hm_data_all, aes(x = Phase, y = pair_id, fill = Z_Score)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0, 
                       name = "Relative\nEnrichment\n(Z-Score)") +
  facet_grid(Significant ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Cell-Cycle Specific Long-Range Chromatin Contacts",
    subtitle = sprintf("Top: Significant (P < 0.05, n=%d) | Bottom: Non-significant (n=%d)",
                       sum(cc_results_df$P_Value < 0.05),
                       sum(cc_results_df$P_Value >= 0.05)),
    x = "Cell Cycle Phase",
    y = "Loci Pairs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_heatmap_facet)

# Option 3: If you want to highlight significant pairs in the heatmap itself
p_heatmap_highlight <- ggplot(hm_data_all, aes(x = Phase, y = pair_id, fill = Z_Score)) +
  geom_tile(aes(color = Significant), size = 0.8) +
  scale_color_manual(values = c("Yes" = "red", "No" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0, 
                       name = "Relative\nEnrichment\n(Z-Score)") +
  labs(
    title = "Cell-Cycle Specific Long-Range Chromatin Contacts",
    subtitle = sprintf("All %d pairs sorted by P-value (Red border = P < 0.05)", 
                       nrow(cc_results_df)),
    x = "Cell Cycle Phase",
    y = "Loci Pairs (sorted by P-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_heatmap_highlight)

# Save plots
ggsave("heatmap_all_pairs_sorted_by_pvalue.pdf", p_heatmap_all, width = 10, height = 20)
ggsave("heatmap_facet_by_significance.pdf", p_heatmap_facet, width = 10, height = 20)
ggsave("heatmap_highlight_significant.pdf", p_heatmap_highlight, width = 10, height = 20)

cat(sprintf("\nPlotted %d total pairs (significant: %d, non-significant: %d)\n", 
            nrow(cc_results_df),
            sum(cc_results_df$P_Value < 0.05),
            sum(cc_results_df$P_Value >= 0.05)))

####################################################################################################
# Fig. S8b  mESC chr1 pseudo-bulk A/B compartments ####
# CHR1 ONLY (TARGET_CHR = 'chr1')
####################################################################################################
# pseudo bulk Hi-C and A/B comparment ####

cat("\n=== Step 1: Generate Pseudo-bulk Matrix for a Target Chromosome ===\n")

TARGET_CHR <- "chr1"
BIN_SIZE <- 250000

pairs_files <- list.files("mesc_pairs/pairs",
                          pattern = "\\.pairs\\.gz$",
                          full.names = TRUE)
print(pairs_files)
print(length(pairs_files))
cat("Loading 3D Genome data...\n")
TARGET_CHR = 'chr1'
# 2. Define a function to process a SINGLE file safely
read_and_filter_chr <- function(file_path) {
  # Read ONLY the structural columns (chr1, pos1, chr2, pos2) to save massive amounts of memory
  # In standard .pairs format, these are usually columns 2, 3, 4, 5
  tmp <- fread(file_path, select = c(2, 3, 4, 5), col.names = c("chr1", "pos1", "chr2", "pos2"))
  
  # Instantly filter for only the target chromosome (intra-chromosomal contacts)
  tmp_filtered <- tmp[chr1 == TARGET_CHR & chr2 == TARGET_CHR]
  
  return(tmp_filtered)
}

# 3. Execute the loop across all files (This will take a few minutes, pbapply shows a progress bar)
list_of_chr_contacts <- pblapply(pairs_files, read_and_filter_chr)

# 4. Combine the list of small data.tables into one giant pseudo-bulk data.table
cat("Combining files...\n")
unfiltered_contacts <- rbindlist(list_of_chr_contacts)

# 5. Apply the binning math
cat("Binning coordinates...\n")
unfiltered_contacts[, bin1 := round(pos1 / BIN_SIZE) * BIN_SIZE]
unfiltered_contacts[, bin2 := round(pos2 / BIN_SIZE) * BIN_SIZE]

cat(sprintf("\nSuccessfully loaded %d contacts for %s.\n", nrow(unfiltered_contacts), TARGET_CHR))

# 1. Aggregate all cells into a single Pseudo-bulk dataset
pseudo_bulk <- unfiltered_contacts %>%
  filter(chr1 == TARGET_CHR, chr2 == TARGET_CHR) %>%
  group_by(bin1, bin2) %>%
  summarize(counts = n(), .groups = "drop")

# 2. Determine the full range of bins for this chromosome
all_bins <- seq(min(c(pseudo_bulk$bin1, pseudo_bulk$bin2), na.rm = TRUE), 
                max(c(pseudo_bulk$bin1, pseudo_bulk$bin2), na.rm = TRUE), 
                by = BIN_SIZE)
n_bins <- length(all_bins)

# 3. Create an empty dense symmetric matrix
contact_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)
rownames(contact_matrix) <- all_bins
colnames(contact_matrix) <- all_bins

# 4. Fill the matrix (both upper and lower triangles to make it symmetric)
idx_row <- match(pseudo_bulk$bin1, all_bins)
idx_col <- match(pseudo_bulk$bin2, all_bins)

for(i in 1:nrow(pseudo_bulk)) {
  contact_matrix[idx_row[i], idx_col[i]] <- pseudo_bulk$counts[i]
  contact_matrix[idx_col[i], idx_row[i]] <- pseudo_bulk$counts[i] # Symmetry
}

cat(sprintf("Pseudo-bulk dense matrix created: %d x %d bins.\n", n_bins, n_bins))

cat("\n=== Step 2: Observed / Expected (O/E) Normalization ===\n")

# Expected contacts depend on the genomic distance between bins
oe_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)
rownames(oe_matrix) <- all_bins
colnames(oe_matrix) <- all_bins

# Calculate expected values for each diagonal (genomic distance)
for (d in 0:(n_bins - 1)) {
  # Get indices for the current diagonal
  diag_indices <- cbind(1:(n_bins - d), (1 + d):n_bins)
  
  # Calculate mean contact frequency at this distance
  expected_val <- mean(contact_matrix[diag_indices], na.rm = TRUE)
  
  # Prevent division by zero
  if (expected_val == 0) expected_val <- 1e-6 
  
  # Divide Observed by Expected
  oe_matrix[diag_indices] <- contact_matrix[diag_indices] / expected_val
  
  # Fill the symmetric lower diagonal
  if (d > 0) {
    sym_indices <- cbind((1 + d):n_bins, 1:(n_bins - d))
    oe_matrix[sym_indices] <- contact_matrix[diag_indices] / expected_val
  }
}

# Log2 transform the O/E matrix (standard practice for compartment calling)
oe_matrix_log <- log2(oe_matrix + 0.01) # Add pseudocount to avoid log2(0)


cat("\n=== Step 3: Pearson Correlation and PCA (PC1) ===\n")

# Calculate Pearson Correlation of the O/E matrix
cor_matrix <- cor(oe_matrix_log, method = "pearson")

# Handle any NAs (which happen if a bin has absolutely zero coverage)
cor_matrix[is.na(cor_matrix)] <- 0

# Perform Principal Component Analysis (PCA)
pca_res <- prcomp(cor_matrix, center = TRUE, scale. = TRUE)

# Extract the First Principal Component (PC1)
pc1_scores <- pca_res$x[, 1]
compartments_df <- data.frame(
  bin = all_bins,
  PC1 = pc1_scores
)


cat("\n=== Step 4: Auto-Orient PC1 using RNA Expression (A/B Assignment) ===\n")

# PC1 is arbitrary in its sign (+ or -). We need to ensure that the Active (A) 
# compartment is positive. We will correlate PC1 with your previously calculated 
# RNA expression to correctly orient it.

# 1. Sum RNA expression for each bin across all cells (using your existing rna_binned_agg_B6)
# Ensure you are only looking at the target chromosome
rna_chr_summary <- rna_binned_agg_B6 %>%
  filter(chr == TARGET_CHR) %>%
  group_by(bin) %>%
  summarize(total_expr = sum(bin_expr), .groups = "drop")

# 2. Merge PC1 scores with RNA expression
compartments_oriented <- compartments_df %>%
  left_join(rna_chr_summary, by = "bin") %>%
  mutate(total_expr = ifelse(is.na(total_expr), 0, total_expr))

# 3. Check correlation. If negative, flip the PC1 sign!
pc1_rna_cor <- cor(compartments_oriented$PC1, compartments_oriented$total_expr, method = "spearman")

if (pc1_rna_cor < 0) {
  cat("Flipping PC1 sign so Positive = A compartment (Active)...\n")
  compartments_oriented$PC1 <- -compartments_oriented$PC1
}

# Assign labels
compartments_oriented <- compartments_oriented %>%
  mutate(Compartment = ifelse(PC1 > 0, "A", "B"))

cat("Compartment calling complete. Summary:\n")
print(table(compartments_oriented$Compartment))


cat("\n=== Step 5: Visualization of A/B Compartments ===\n")

# Create a classic track-style bar plot for A/B compartments
p_compartments <- ggplot(compartments_oriented, aes(x = bin / 1e6, y = PC1, fill = Compartment)) +
  geom_bar(stat = "identity", width = (BIN_SIZE / 1e6)) +
  scale_fill_manual(values = c("A" = "#E64B35", "B" = "#4DBBD5")) +
  labs(
    title = paste("A/B Compartments for", TARGET_CHR),
    subtitle = "Positive PC1 indicates the active A compartment (RNA-oriented)",
    x = "Genomic Position (Mb)",
    y = "PC1 Score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_compartments)

setDT(compartments_oriented)

pc1_min <- min(compartments_oriented$PC1, na.rm = TRUE)
pc1_max <- max(compartments_oriented$PC1, na.rm = TRUE)

# Normalize to -1.5 ~ 1.5
compartments_oriented[, PC1_norm := (PC1)]
summary(compartments_oriented$PC1_norm)

compartments_oriented[, Compartment_new := fcase(
  PC1_norm > 0.3, "A",
  PC1_norm < -0.3, "B",
  between(PC1_norm, -0.3, 0.3), "Boundary"
)]
print(table(compartments_oriented$Compartment))
valid_bins <- unique(c(valid_pairs$bin1, valid_pairs$bin2))
compartments_oriented[, in_valid_pairs := ifelse(bin %in% valid_bins, TRUE, FALSE)]
print(table(compartments_oriented$in_valid_pairs))

colors <- c("All bins" = "#1A237E", "Valid pairs bins" = "#2E7D32")
ggplot() +
  geom_density(
    data = compartments_oriented,
    aes(x = PC1_norm, fill = "All bins", color = "All bins"),
    alpha = 0.5, adjust = 1,
    linewidth = 0.8
  ) +
  geom_density(
    data = compartments_oriented[in_valid_pairs == TRUE, ],
    aes(x = PC1_norm, fill = "Valid pairs bins", color = "Valid pairs bins"),
    alpha = 0.5, adjust = 1,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    name = "Bins",
    values = c("All bins" = "#1A237E", "Valid pairs bins" = "#2ca02c")
  ) +
  scale_color_manual(
    name = "Bins",
    values = c("All bins" = "#1A237E", "Valid pairs bins" = "#2ca02c")
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Distribution of PC1 values",
    subtitle = "Deep blue: all bins | Dark green: bins in valid pairs",
    x = "Normalized PC1 (-1.5 ~ 1.5)",
    y = "Density"
  ) +
  theme(
    legend.position = "top",
    legend.box = "horizontal"
  )

compartments_oriented[, PC1_zscore := scale(PC1)]
compartments_oriented[, Compartment := fcase(
  PC1_zscore > 0.3, "A",
  PC1_zscore < -0.3, "B",
  between(PC1_zscore, -0.3, 0.3), "Boundary"
)]

comp_lookup <- compartments_oriented[, .(bin, Compartment)]
setkey(comp_lookup, bin)
valid_pairs[, comp1 := comp_lookup[.(bin1), Compartment]]
valid_pairs[, comp2 := comp_lookup[.(bin2), Compartment]]
cat("Pairs with missing compartment info:\n")
print(valid_pairs[is.na(comp1) | is.na(comp2), .N])

valid_pairs[, pair_type := fcase(
  comp1 == "A" & comp2 == "A", "A-A",
  comp1 == "B" & comp2 == "B", "B-B",
  (comp1 == "A" & comp2 == "B") | (comp1 == "B" & comp2 == "A"), "A-B",
  comp1 == "Boundary" | comp2 == "Boundary", "Contains_Boundary",
  default = "Other"
)]

result <- valid_pairs[, .N, by = pair_type]
result[, Percentage := round(N / sum(N) * 100, 2)]

result <- result[order(factor(pair_type, levels = c("A-A", "B-B", "A-B", "Contains_Boundary", "Other"))), ]

total_row <- data.table(
  pair_type = "Total",
  N = sum(result$N),
  Percentage = 100
)
result <- rbind(result, total_row)

print(result)
