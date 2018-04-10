#directory_membrane
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(xlsx)
library(readxl)

#step1: read peptides_nogly
#step2: remove the rows begin with CON and REV at "Leading razor protein" column
#step3: read homo sapines transmembrane databse
#step4: filter by keywords: "TRANSMEM" at "Transmembrane" column
#step5: calculate "TRANSMEM" keywords frequency in each row
#step6: left_join peptides_nogly_remove with homo_tran by (leading.razor.protein = Entry)
#step7: remove is.na(peptide_nogly_combine_homo_tran) and barplot
#step8: filter: tran_num == 1 and barplot

dir <- "E:\\Program Files\\R program practice\\gao_wn\\gwn_git_transmembrane\\Pancreatic_project"
if(!dir.exists(dir))
{
  dir.create(dir)
}

#step1: read table
peptides_nogly <- read.table(paste(dir,"\\peptides_nogly.txt", sep = ""),sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#step2: remove CON and REV
peptides_nogly_remove <- peptides_nogly[!grepl("^CON|^REV",peptides_nogly$Leading.razor.protein),]

#step3: read homo sapines transmembrane databse
h_tran_db <- read.delim(paste(dir,"\\Homo_sapiens(162191)_transmembrane.txt",sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#step4: filter by keywords: "TRANSMEM" at "Transmembrane" column
h_tran <- h_tran_db[grepl("^TRANSMEM", h_tran_db$Transmembrane),]

#step5: calculate "TRANSMEM" keywords frequency in each row
tran_num <- sapply(h_tran$Transmembrane,function(x){length(strsplit(x,"TRANSMEM")[[1]])-1})

h_tran$tran_num <- tran_num

#step6: left_join peptides_nogly_remove with homo_tran by (leading.razor.protein = Entry)
peptide_combine_tran <- left_join(peptides_nogly_remove,h_tran, by = c("Leading.razor.protein" = "Entry"))

name_col <- c("Leading.razor.protein", "Start.position", "End.position","Transmembrane","tran_num",
              "Proteins","Protein.names.x","Gene.names.x","Protein.names.y","Gene.names.y")

peptide_combine_tran1 <- peptide_combine_tran[,colnames(peptide_combine_tran) %in% name_col]

#step7: remove is.na(peptide_nogly_combine_homo_tran)
peptide_combine_tran2 <- peptide_combine_tran1[!is.na(peptide_combine_tran1$tran_num),]

frequency_num <-  table(peptide_combine_tran2$tran_num)

pdf(paste(dir, "\\barplot_transmembrane_number.pdf",sep = ""),width = 12)
barplot(frequency_num,main = "Transmembrane number", xlab = "transmembrane num", ylab = "Frequency")
dev.off()

#step8:Analyze only one transmembrane 
peptide_num1 <- peptide_combine_tran2[peptide_combine_tran2$tran_num==1,]

pep_tran <- str_extract_all(peptide_num1$Transmembrane,"[0-9]+[0-9]")

tran_start <- c();tran_end <- c()
for(i in 1:nrow(peptide_num1))
{
  tran_start[i] <- pep_tran[[i]][1]
  tran_end[i] <- pep_tran[[i]][2]
}

peptide_num1$tran_start <- as.numeric(tran_start)
peptide_num1$tran_end <- as.numeric(tran_end)

peptide_num1$tran_location <- "location"
pos_locate_v2 <- function(peptide_num1)
{
  for(i in 1:nrow(peptide_num1))
  {
    if(peptide_num1$End.position[i] < peptide_num1$tran_start[i])
    {peptide_num1$tran_location[i] <- "1_extracellular"}
    else if(peptide_num1$End.position[i] >= peptide_num1$tran_start[i] & peptide_num1$End.position[i] <= peptide_num1$tran_end[i] & peptide_num1$Start.position[i] < peptide_num1$tran_start[i])
    {peptide_num1$tran_location[i] <- "2_extra_trans"}
    else if(peptide_num1$Start.position[i] >= peptide_num1$tran_start[i] & peptide_num1$End.position[i] <= peptide_num1$tran_end[i])
    {peptide_num1$tran_location[i] <- "3_transmembrane"}
    else if(peptide_num1$Start.position[i] >= peptide_num1$tran_start[i] & peptide_num1$Start.position[i] <= peptide_num1$tran_end[i] & peptide_num1$End.position[i] > peptide_num1$tran_end[i])
    {peptide_num1$tran_location[i] <- "4_cyto_trans"}
    else if(peptide_num1$Start.position[i] < peptide_num1$tran_start[i] & peptide_num1$End.position[i] > peptide_num1$tran_end[i])
    {peptide_num1$tran_location[i] <- "6_cyto_trans_extra"}
    else
    {peptide_num1$tran_location[i] <- "5_cytoplasmic"}
  }
  return(peptide_num1)
}

peptide_num3 <- pos_locate_v2(peptide_num1)

write.table(peptide_num3, file = paste(dir, "\\tran_num_only1.txt",sep = ""),sep = "\t",row.names = FALSE,quote = FALSE)

peptide_num4 <- table(peptide_num3$tran_location)

pdf(paste(dir, "\\barplot_transmembrane_distributuion_number.pdf",sep = ""),width = 15)
barplot(peptide_num4,main = "distribution of peptide sequence", xlab = "zone", ylab = "Frequency")
dev.off()

#step8_EGFR and PDGFR intensity
rtk <- read.csv(paste(dir,"\\RTK.csv",sep = ""),sep = ",",header = TRUE,stringsAsFactors = FALSE)

rtk_tran1 <- right_join(peptide_num3,rtk,by = c("Leading.razor.protein" = "Protein"))

rtk_tran2 <- rtk_tran1[!is.na(rtk_tran1$Start.position),]

rtk_pdffra <- rtk_tran2[grepl("PDGFRA",rtk_tran2$Gene.names.x),]
rtk_pdffra_pos <- rtk_pdffra[,c(2,5,3,4,11:13)]
pdgfra_table <- table(rtk_pdffra_pos$tran_location)
pdgfra_table <- as.data.frame(pdgfra_table)

pdgfra_join <- left_join(rtk_pdffra_pos,peptide_combine_tran)
pdgfra_join <- pdgfra_join[-nrow(pdgfra_join),]

pdgfr_join_extra <- pdgfra_join[pdgfra_join$tran_location == "1_extracellular",]
pdgfr_extra_intensity <- select(pdgfr_join_extra,matches("intensity\\."))
pdgfr_extra_intensity[pdgfr_extra_intensity ==0] <- NA 
pdgfr_extra_mean <- apply(pdgfr_extra_intensity,2,function(x){mean(x,na.rm = TRUE)})  
  
pdgfr_join_cyto <- pdgfra_join[pdgfra_join$tran_location == "5_cytoplasmic",] 
pdgfr_cyto_intensity <- select(pdgfr_join_cyto,matches("intensity\\."))
pdgfr_cyto_intensity[pdgfr_cyto_intensity == 0] <- NA
pdgfr_cyto_mean <- apply(pdgfr_cyto_intensity,2,function(x){mean(x,na.rm = TRUE)})
pdgfra_ratio <- pdgfr_extra_mean/pdgfr_cyto_mean

rtk_pdffrb <- rtk_tran2[grepl("PDGFRB",rtk_tran2$Gene.names.x),]
rtk_pdffrb_pos <- rtk_pdffrb[,c(2,5,3,4,11:13)]
pdgfrb_table <- table(rtk_pdffrb_pos$tran_location)
pdgfrb_table <- as.data.frame(pdgfrb_table)


rtk_egfr <- rtk_tran2[grepl("EGFR",rtk_tran2$Gene.names.x),]
rtk_egfr_pos <- rtk_egfr[,c(2,5,3,4,11:13)]
egfr_table <- table(rtk_egfr_pos$tran_location)
egfr_table <- as.data.frame(egfr_table)

#step8_total intensity
#join_table
peptide_total_join <- left_join(peptide_num3,peptide_combine_tran,by = c("Leading.razor.protein","Start.position","End.position"))


peptide_total_join2 <- select(peptide_total_join,matches("Leading.razor.protein|Gene.names|position|tran_location|intensity\\.")) 
pep_pro_gene <- select(peptide_total_join2,matches("Leading.razor.protein|Gene.names\\.x\\.x"))
pep_pro_gene <- pep_pro_gene[!duplicated(pep_pro_gene$Leading.razor.protein),]

unique_protein <- unique(peptide_total_join2$Leading.razor.protein)

#i <- unique_protein[1]
pep_total_res <- data.frame()
for(i in c(unique_protein))
{
  pep_extra <- peptide_total_join2[peptide_total_join2$tran_location == "1_extracellular" & peptide_total_join2$Leading.razor.protein %in% i,]
  pep_cyto <-  peptide_total_join2[peptide_total_join2$tran_location == "5_cytoplasmic" & peptide_total_join2$Leading.razor.protein %in% i,]
  
  pep_extra_intensity <- select(pep_extra,matches("intensity\\.."))
  pep_extra_notzero <- apply(pep_extra_intensity,2,function(x){sum(x!=0)})
  pep_extra_notzero <- as.data.frame(pep_extra_notzero)
  pep_extra_notzero <- t(pep_extra_notzero)
  colnames(pep_extra_notzero) <- gsub("Intensity\\.","extra_num_",colnames(pep_extra_notzero))
  
  if(nrow(pep_extra_intensity)!=0)
  {pep_extra_intensity[pep_extra_intensity == 0] <- NA}
  
  pep_cyto_intensity <- select(pep_cyto,matches("intensity\\.."))
  pep_cyto_notzero <- apply(pep_cyto_intensity,2,function(x){sum(x!=0)})
  pep_cyto_notzero <- as.data.frame(pep_cyto_notzero)
  pep_cyto_notzero <- t(pep_cyto_notzero)
  colnames(pep_cyto_notzero) <- gsub("Intensity\\.","cyto_num_",colnames(pep_cyto_notzero)) 
  
  if(nrow(pep_cyto_intensity)!=0)
  {pep_cyto_intensity[pep_cyto_intensity == 0] <- NA}
  
  pep_extra_mean <- apply(pep_extra_intensity,2,function(x){mean(x,na.rm = TRUE)})

  pep_cyto_mean <- apply(pep_cyto_intensity,2,function(x){mean(x,na.rm = TRUE)})

  pep_e_c <- pep_extra_mean/pep_cyto_mean
  
  pep_extra_mean <- as.data.frame(pep_extra_mean)
  pep_extra_mean <- t( pep_extra_mean)
  colnames(pep_extra_mean) <- gsub("Intensity","Extra_Intensity",colnames(pep_extra_mean))
  
  pep_cyto_mean <- as.data.frame(pep_cyto_mean)
  pep_cyto_mean <- t(pep_cyto_mean)
  colnames(pep_cyto_mean) <- gsub("Intensity", "Cyto_Intensity",colnames(pep_cyto_mean))
  
  pep_e_c <- as.data.frame(pep_e_c)
  pep_e_c <- t(pep_e_c)
  colnames(pep_e_c) <- gsub("Intensity\\.","ratio_",colnames(pep_e_c))
  
  pep_res <- data.frame(Leading.razor.protein = i, extra_num <- nrow(pep_extra), cyto_num <- nrow(pep_cyto),pep_e_c,pep_extra_notzero,pep_cyto_notzero,pep_extra_mean,pep_cyto_mean,stringsAsFactors = FALSE)
  pep_total_res <- rbind(pep_total_res,pep_res)
}

pep_total_res2 <- left_join(pep_total_res,pep_pro_gene,by = "Leading.razor.protein")
pep_total_res2 <- pep_total_res2[pep_total_res2$extra_num....nrow.pep_extra.!=0 | pep_total_res2$cyto_num....nrow.pep_cyto.!=0,]

pep_total_res3 <- pep_total_res2[,c(1,109,2,3)]
for(i in 4:24)
{
  pep_total_res3 <- cbind(pep_total_res3,pep_total_res2[,c(i,i+21,i+42,i+63,i+84)])
}


for(i in 1:nrow(pep_total_res2))
{
  if(pep_total_res2$extra_num....nrow.pep_extra.[i]==0)
  {
    pep_total_res2[i,c(4:24)] <- "extra_not_identified"
  }
  else if(pep_total_res2$cyto_num....nrow.pep_cyto.[i]==0)
  {
    pep_total_res2[i,c(4:24)] <- "cyto_not_identified"
  }
  # else if(pep_total_res2$extra_num....nrow.pep_extra.[i]!=0 & pep_total_res2$cyto_num....nrow.pep_cyto.[i]!=0)
  # {
  #   pep_total_res2[i,c(4:24)][is.na(pep_total_res2[i,c(4:24)])] <- "not_identified"
  # }
}


write.table(pep_total_res2,file = paste(dir,"\\peptide_intensity_ratio.txt",sep = ""),sep = "\t",row.names = FALSE,quote = FALSE)

write.table(pep_total_res3,file = paste(dir,"\\peptide_intensity_ratio_v2.txt",sep = ""),sep = "\t",row.names = FALSE,quote = FALSE)





#step9:Analyze the peptide with transmembrane_num >=2 

peptide_combine_tran2 <- peptide_combine_tran2[peptide_combine_tran2$tran_num!=1,]

pep_tran_list <- lapply(peptide_combine_tran2$Transmembrane,function(x){gregexpr("TRANSMEM\\s{1}\\d{1,4}\\s{1}\\d{1,4}",x)})

grep_start <- list(); grep_add <- list();grep_stop <- list();tran_res <- list()
for(i in 1:nrow(peptide_combine_tran2))
{
  grep_start[[i]] <- unlist(pep_tran_list[[i]])
  grep_add[[i]] <- attr(pep_tran_list[[i]][[1]],'match.length')
  grep_stop[[i]] <- grep_start[[i]] + grep_add[[i]]
  tran_res[[i]] <- substring(peptide_combine_tran2$Transmembrane[i],grep_start[[i]],grep_stop[[i]])
}

pep_tran_list2 <- sapply(tran_res,function(x){str_extract_all(unlist(x),"\\d+")})

pep_tran_list3 <- list(); pep_tran_list4 <- list()
for(i in 1:nrow(peptide_combine_tran2))
{
  pep_tran_list3[[i]] <- unlist(pep_tran_list2[[i]])
  pep_tran_list3[[i]] <- as.numeric(pep_tran_list3[[i]])
  
}

complete_interval <- function(x)
{
  lenx <- length(x)
  seq_oddnum <- seq(1,length(x),by = 2)
  seq_evenum <- seq(2,length(x), by =2)
  odd_num <- x[seq_oddnum]-1
  even_num <- x[seq_evenum]+1
  combine_num <- 1
  
  for(i in 1:length(odd_num))
  {
    combine_num <- c(combine_num,odd_num[i],even_num[i])
  }  
  combine_num <- c(combine_num,9999)
  
  extra_1st <- seq(1,length(combine_num),4)
  extra_last <- seq(2,length(combine_num),4)
  cyto_1st <- seq(3,length(combine_num),4)
  cyto_last <- seq(4,length(combine_num),4)
  
  extra_interval <- c()
  for(i in 1:length(extra_1st))
  {
    extra_interval <- c(extra_interval,combine_num[extra_1st][i],combine_num[extra_last][i])
  }
  
  cyto_interval <- c()
  for(i in 1:length(cyto_1st))
  {
    cyto_interval <- c(cyto_interval,combine_num[cyto_1st][i],combine_num[cyto_last][i])  
  }
  return(list(x,extra_interval,cyto_interval))
}

#================================end complete_interval====================================

#temp_interval <- complete_interval(pep_tran_list3[[4]])

total_interval <- lapply(pep_tran_list3,function(x){complete_interval(unlist(x))})
#total_interval contains transmembrane; extracellular; cytoplasmic

#total_interval[[3]][[3]]
#i <- 1
#tran_num[3]

interval_distribution_transmembrane <- function(total_interval)
{
  tran_num <- c()
  for(i in 1:nrow(peptide_combine_tran2))
  {
    tran_num[i] <- length(total_interval[[i]][[1]])
  }
  
  col_num <- max(tran_num,na.rm = TRUE)
  
  #cost long time
  temp_vec <- rep(0,col_num)
  temp_mat <- t(data.frame(temp_vec))
  #temp_mat <- matrix(0,nrow = 1,ncol =col_num )
  for(i in 1:nrow(peptide_combine_tran2))
  {
    temp_mat <- rbind(temp_mat,total_interval[[i]][[1]])
  }
  temp_mat1 <- temp_mat[-1,]
  temp_mat1 <- as.data.frame(temp_mat1)
  
  odd_num <- seq(from = 1,to = col_num,by = 2) 
  even_num <- seq(from = 2,to = col_num,by = 2)
  
  temp_mat2 <- matrix(0,nrow = nrow(peptide_combine_tran2),ncol = 1)
  temp_array <- c()
  
  for(i in c(odd_num))
  {
    temp_array <- paste(temp_mat1[,i],temp_mat1[,i+1],sep = "_")
    temp_mat2 <- cbind(temp_mat2,temp_array)
  }
  temp_mat2 <- temp_mat2[,-1]
  
  peptide_combine_tran3 <- cbind(peptide_combine_tran2,temp_mat2)
  
  col_name <- c("Leading.razor.protein","Start.position","End.position","tran_num","temp_array")
  
  peptide_combine_tran4 <- peptide_combine_tran3[,names(peptide_combine_tran3) %in% col_name]
  peptide_combine_tran4[,c(5:((col_num)/2+4))] <- apply(peptide_combine_tran4[,c(5:((col_num)/2+4))],2,as.character)
  
  df_melt <- melt(peptide_combine_tran4,id = c("Leading.razor.protein","Start.position","End.position","tran_num"))
  
  df_melt$variable <- NULL
  df_melt2 <- unique(df_melt)
  
  df_melt2$tran_start <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[1]}))
  df_melt2$tran_end <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[2]}))
  
  df_melt2$tran_start <- as.numeric(df_melt2$tran_start)
  df_melt2$tran_end <- as.numeric(df_melt2$tran_end)
  return(df_melt2)
}


melt_transmembrane <- interval_distribution_transmembrane(total_interval)

interval_distribution_extracellular <- function(total_interval)
{
  tran_num <- c()
  for(i in 1:nrow(peptide_combine_tran2))
  {
    tran_num[i] <- length(total_interval[[i]][[2]])
  }
  
  col_num <- max(tran_num)
  
  #cost long time
  temp_vec <- rep(0,col_num)
  temp_mat <- t(data.frame(temp_vec))
  #temp_mat <- matrix(0,nrow = 1,ncol =col_num )
  for(i in 1:nrow(peptide_combine_tran2))
  {
    temp_mat <- rbind(temp_mat,total_interval[[i]][[2]])
  }
  temp_mat1 <- temp_mat[-1,]
  temp_mat1 <- as.data.frame(temp_mat1)
  
  odd_num <- seq(from = 1,to = col_num,by = 2) 
  even_num <- seq(from = 2,to = col_num,by = 2)
  
  temp_mat2 <- matrix(0,nrow = nrow(peptide_combine_tran2),ncol = 1)
  temp_array <- c()
  
  for(i in c(odd_num))
  {
    temp_array <- paste(temp_mat1[,i],temp_mat1[,i+1],sep = "_")
    temp_mat2 <- cbind(temp_mat2,temp_array)
  }
  temp_mat2 <- temp_mat2[,-1]
  
  peptide_combine_tran3 <- cbind(peptide_combine_tran2,temp_mat2)
  
  
  col_name <- c("Leading.razor.protein","Start.position","End.position","tran_num","temp_array")
  
  peptide_combine_tran4 <- peptide_combine_tran3[,names(peptide_combine_tran3) %in% col_name]
  peptide_combine_tran4[,c(5:((col_num)/2+4))] <- apply(peptide_combine_tran4[,c(5:((col_num)/2+4))],2,as.character)
  
  df_melt <- melt(peptide_combine_tran4,id = c("Leading.razor.protein","Start.position","End.position","tran_num"))
  
  df_melt$variable <- NULL
  df_melt2 <- unique(df_melt)
  
  df_melt2$tran_start <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[1]}))
  df_melt2$tran_end <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[2]}))
  
  df_melt2$tran_start <- as.numeric(df_melt2$tran_start)
  df_melt2$tran_end <- as.numeric(df_melt2$tran_end)
  return(df_melt2)
}
melt_extra <- interval_distribution_extracellular(total_interval)

interval_distribution_cytoplasmic <- function(total_interval)
{
  tran_num <- c()
  for(i in 1:nrow(peptide_combine_tran2))
  {
    tran_num[i] <- length(total_interval[[i]][[3]])
  }
  
  col_num <- max(tran_num)
  
  #cost long time
  temp_vec <- rep(0,col_num)
  temp_mat <- t(data.frame(temp_vec))
  for(i in 1:nrow(peptide_combine_tran2))
  {
    temp_mat <- rbind(temp_mat,total_interval[[i]][[3]])
  }
  temp_mat1 <- temp_mat[-1,]
  temp_mat1 <- as.data.frame(temp_mat1)
  
  odd_num <- seq(from = 1,to = col_num,by = 2) 
  even_num <- seq(from = 2,to = col_num,by = 2)
  
  temp_mat2 <- matrix(0,nrow = nrow(peptide_combine_tran2),ncol = 1)
  temp_array <- c()
  
  for(i in c(odd_num))
  {
    temp_array <- paste(temp_mat1[,i],temp_mat1[,i+1],sep = "_")
    temp_mat2 <- cbind(temp_mat2,temp_array)
  }
  temp_mat2 <- temp_mat2[,-1]
  
  peptide_combine_tran3 <- cbind(peptide_combine_tran2,temp_mat2)
  
  col_name <- c("Leading.razor.protein","Start.position","End.position","tran_num","temp_array")
  
  peptide_combine_tran4 <- peptide_combine_tran3[,names(peptide_combine_tran3) %in% col_name]
  peptide_combine_tran4[,c(5:((col_num)/2+4))] <- apply(peptide_combine_tran4[,c(5:((col_num)/2+4))],2,as.character)
  
  df_melt <- melt(peptide_combine_tran4,id = c("Leading.razor.protein","Start.position","End.position","tran_num"))
  
  df_melt$variable <- NULL
  df_melt2 <- unique(df_melt)
  
  df_melt2$tran_start <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[1]}))
  df_melt2$tran_end <- unlist(lapply(df_melt2$value,function(x){unlist(strsplit(x,split = "_"))[2]}))
  
  df_melt2$tran_start <- as.numeric(df_melt2$tran_start)
  df_melt2$tran_end <- as.numeric(df_melt2$tran_end)
  return(df_melt2)
}
melt_cyto <- interval_distribution_cytoplasmic(total_interval)


#transmembrane_location annotation
melt_transmembrane_cp <- melt_transmembrane
melt_transmembrane_cp$location <- "unknown"
for(i in 1:nrow(melt_transmembrane_cp))
{
  if(melt_transmembrane_cp$Start.position[i] >= melt_transmembrane_cp$tran_start[i] & melt_transmembrane_cp$End.position[i] <= melt_transmembrane_cp$tran_end[i])
  {
    melt_transmembrane_cp$location[i] <- "3_transmembrane"
  }
}
melt_transmembrane_cp$value <- NULL
melt_transmembrane_cp$tran_start <- NULL
melt_transmembrane_cp$tran_end <- NULL
melt_transmembrane_uni <- unique(melt_transmembrane_cp)

#extracellular_location annotation
melt_extra_cp <- melt_extra
melt_extra_cp$location <- "unknown"
for(i in 1:nrow(melt_extra_cp))
{
  if(melt_extra_cp$Start.position[i] >= melt_extra_cp$tran_start[i] & melt_extra_cp$End.position[i] <= melt_extra_cp$tran_end[i])
  {
    melt_extra_cp$location[i] <- "1_extracellular"
  }
}
melt_extra_cp$value <- NULL
melt_extra_cp$tran_start <- NULL
melt_extra_cp$tran_end <- NULL
melt_extra_cp_uni <- unique(melt_extra_cp)


#cytoplasmic_location annotation
melt_cyto_cp <- melt_cyto
melt_cyto_cp $location <- "unknown"
for(i in 1:nrow(melt_cyto_cp))
{
  if(melt_cyto_cp$Start.position[i] >= melt_cyto_cp$tran_start[i] & melt_cyto_cp$End.position[i] <= melt_cyto_cp$tran_end[i])
  {
    melt_cyto_cp$location[i] <- "5_cytoplasmic"
  }
}
melt_cyto_cp$value <- NULL
melt_cyto_cp$tran_start <- NULL
melt_cyto_cp$tran_end <- NULL
melt_cyto_cp_uni <- unique(melt_cyto_cp)

#combine three result
combine_res <- rbind(melt_transmembrane_uni,melt_extra_cp_uni,melt_cyto_cp_uni)
combine_res_uni <- unique(combine_res)

combine_res_anno <- combine_res_uni[combine_res_uni$location != "unknown",]

col_join_name <- colnames(combine_res_anno)[1:4]
peptide_combine_tran2_join <- peptide_combine_tran2[,names(peptide_combine_tran2) %in% col_join_name]

temp_join <- left_join(peptide_combine_tran2_join,combine_res_anno,by=col_join_name)

temp_join$location[is.na(temp_join$location)] <- "unknown"

write.table(temp_join, file = paste(dir,"\\tran_num_not1.txt",sep = ""), sep = "\t",row.names = FALSE, quote = FALSE)

table_anno_interval <- table(temp_join$location)

pdf(paste(dir, "\\barplot_interval_distribution.pdf",sep = ""),width = 12)
barplot(table_anno_interval,main = "interval_distribution", xlab = "interval_type", ylab = "Frequency")
dev.off()

#combine tran_num == 1 and tran_num !=1 result
#combine the two results of peptide_num3 and temp_join
col_join_name2 <- c(col_join_name, "tran_location")
tran_one <- peptide_num3[,names(peptide_num3) %in% col_join_name2]
colnames(tran_one)[5] <- "location"

tran_above_one <- temp_join[order(temp_join$tran_num),] 
tran_bind <- rbind(tran_one, tran_above_one)

#read topology domain
homo <- read_xlsx(paste(dir,"\\Homo_sapiens(20316)_transmembrane_topology.xlsx",sep = ""),sheet = 1)
homo_tran <- homo[grepl("TRANSMEM",homo$Transmembrane),]

#join tran_bind with topology

tran_bind_topo <- left_join(tran_bind,homo_tran,by = c("Leading.razor.protein" = "Entry"))
tran_bind_topo_noNA <- tran_bind_topo[!is.na(tran_bind_topo$`Topological domain`),]
colnames(tran_bind_topo_noNA)

col_join_name3 <- c("Leading.razor.protein", "Start.position", "End.position", "tran_num", "location", "Transmembrane", "Topological domain")
tran_bind_topo_short <- tran_bind_topo_noNA[,names(tran_bind_topo_noNA) %in% col_join_name3]

#topology_extracellular
tran_bind_topo_ex <- tran_bind_topo_short[grepl("Extracellular", tran_bind_topo_short$`Topological domain`),]


tran_bind_topo_ex_fun <- function(peptide_combine_tran2)
{ 
  pep_tran_list <- lapply(peptide_combine_tran2$Transmembrane,function(x){gregexpr("\\d{1,4}\\s{1}\\d{1,4}\\s{1}Extracellular",x)})
  
  grep_start <- list(); grep_add <- list();grep_stop <- list();tran_res <- list()
  for(i in 1:nrow(peptide_combine_tran2))
  {
    grep_start[[i]] <- unlist(pep_tran_list[[i]])
    grep_add[[i]] <- attr(pep_tran_list[[i]][[1]],'match.length')
    grep_stop[[i]] <- grep_start[[i]] + grep_add[[i]]
    tran_res[[i]] <- substring(peptide_combine_tran2$Transmembrane[i],grep_start[[i]],grep_stop[[i]])
  }
  
  pep_tran_list2 <- sapply(tran_res,function(x){str_extract_all(unlist(x),"\\d+")})
  
  pep_tran_list3 <- list(); pep_tran_list4 <- list()
  for(i in 1:nrow(peptide_combine_tran2))
  {
    pep_tran_list3[[i]] <- unlist(pep_tran_list2[[i]])
    pep_tran_list3[[i]] <- as.numeric(pep_tran_list3[[i]])
  }
  return(pep_tran_list3)
}