## This function loads libraries
load_required_libraries <- function(){
  require(dplyr)
  require(reshape2)
  require(manipulateWidget)
  require(manipulate)
  require(plotly)
  require(ggplot2)
  require(limma)
  require(bio3d)
  require(seqinr)
  require(xlsx)
  require(png)
  require(grid)
  require(sequences)
  require(seqinr)
  require(muscle)
  require(Biostrings)
}

source("hetero_dimer_info.R") ## This file contains information about your heterodimer
## write the kernel for the Gaussian convolution
write_convolution_kernel <- function(sse_x,sse_y, file){
  # needs secondary structure of the proteins and the output file
  res <- matrix(1,length(sse_y),length(sse_x))
  for (y in 1:nrow(res)) {
    for (x in 1:ncol(res)) {
      if (sse_y[y] == sse_x[x]) res[y,x] <- 2 
    }
  }
  colnames(res)  <- 1:ncol(res)
  row.names(res) <- 1:nrow(res)
  write.table(res,file,sep = "\t",row.names = F, col.names = F, quote = F)
}

## Predicts the secondary stucture of a protein from fasta file or sequence
predict_sse_from_fasta_psipred <- function(fasta, type = "file") {
  # fasta: fasta file or sequence
  # type: file or sequence
  # Predicts secondary structure of protein from fasta or simple sequence using PSIPRED (no plot)
  require(limma)
  require(bio3d)
  if (type == "sequences") {
    file     <- "/tmp/fasta_file.fasta"
    id       <- "WHATEVER"
    bio3d::write.fasta(ids = id, seqs = fasta, file = file)
    cmd_r    <- paste("fasta_file.ss fasta_file.ss2 fasta_file.horiz",sep = "")
    sse_file <- "fasta_file.ss2"
  }
  else if (type == "file") {
    file  <- fasta
    str1  <- strsplit2(file,"/")[,ncol(strsplit2(file,"/"))] 
    dom   <- strsplit2(str1,"[.]")[,1]
    cmd_r <- paste(dom,c("ss","ss2","horiz"),sep = ".")
    sse_file <- cmd_r[2]
  }
  else (stop("please provide the type of the sequences"))
  cmd <- paste("runpsipred_single",file, sep = " ")
  system(cmd)
  sse <- read.table(sse_file)
  file.remove(cmd_r)
  sse     <- cbind.data.frame(pos = sse$V1, sse = sse$V3)
  sse2    <- as.character(sse$sse)
  sse2
}

## Split the DCA output file to have the region on interest
split_dca_matrix <- function(dca,l) {
  # l is the length of the first protein in the cMSA
  # return list of interacting prot (AA,AB,BB)
  colnames(dca) <- c("X1","X2","X3","X4","X5","X6")
  dat     <- dca %>% select(X1,X3,X6) %>% mutate(X6 = abs(X6)); colnames(dat) <- c("V1","V2","V3")
  t_a_a   <- dat %>% filter(V1 <= l & V2 <= l)
  t_a_a_1 <- t_a_a[,c(2,1,3)]; colnames(t_a_a_1) <- c("V1","V2","V3")
  t_a_a   <- rbind(t_a_a,t_a_a_1)
  t_a_b   <- dat  %>% filter(V1 <= l & V2 >  l) %>% mutate(V2 = V2 - l) 
  t_b_a   <- dat  %>% filter(V1 >  l & V2 <= l) %>% mutate(V1 = V1 - l)
  t_b_b   <- dat  %>% filter(V1 > l & V2 > l)   %>% mutate(V1 = V1 - l, V2 = V2 - l)
  t_b_b_1 <- t_b_b[,c(2,1,3)]; colnames(t_b_b_1) <- c("V1","V2","V3")
  t_b_b   <- rbind(t_b_b,t_b_b_1)
  res     <- list(AA = t_a_a, AB = t_a_b, BA = t_b_a, BB = t_b_b)
  res
}

### Interactive plot the dca. It runs with an external c++ program that should te compiled and added to the PATH
#### Pispred compiler should be downloaded and added to th PATH too
plot_dca_convolution    <- function(dimer_info, min_signal = 0.001, type = "AB") {
  # plots the DCA matrix and uses the manipulate widget to shape the parameters
  # 
  out_dca_file <- "/tmp/dca_con_file.txt"
  out_sse_file <- "/tmp/sse_con_file.txt"
  res_con_file <- "/tmp/res_conv_file.txt"
  
  fasta_x   <- dimer_info$sequence_A                #  get the sequences to compute the lengths of the proteins
  fasta_y   <- dimer_info$sequence_B  
  seq_x     <- bio3d::read.fasta(fasta_x, to.upper = T)$ali
  seq_y     <- bio3d::read.fasta(fasta_y, to.upper = T)$ali
  l         <- length(seq_x)      # length protein X
  l_y       <- length(seq_y)      # length protein Y 
  names     <- dimer_info$gene_names                # read the protein symbol
  
  sse_x <- predict_sse_from_fasta_psipred(fasta_x)  
  sse_y <- predict_sse_from_fasta_psipred(fasta_y)
  
  dca_dat <- read.table(dimer_info$dca_signal)
  if (type == "AB") {dat_split   <- split_dca_matrix(dca_dat,l = l)$AB; xla   <- paste(names[1]); yla   <- paste(names[2])}
  if (type == "AA") {dat_split   <- split_dca_matrix(dca_dat,l = l)$AA; xla   <- paste(names[1]); yla   <- paste(names[1]); sse_y <- sse_x; l_y <- l}
  if (type == "BB") {dat_split   <- split_dca_matrix(dca_dat,l = l)$BB; xla   <- paste(names[2]); yla   <- paste(names[2]); sse_x <- sse_y; l   <- l_y}
  mi_EC <- min(dat_split$V3) # min and max EC signal 
  ma_EC <- max(dat_split$V3)
  
  cat("length protein x: ",l ,"length protein y: ",l_y,"\n")
  
  dca            <- dcast(dat_split,V2 ~ V1, value.var = 'V3')
  dca            <- dca[,-1]
  
  write.table(dca,out_dca_file,sep = "\t",row.names = F, col.names = F, quote = F)
  write_convolution_kernel(sse_x,sse_y, out_sse_file) # Uses external c function for the convolution
  manipulateWidget({ 
    nb_rows <- round(nrow(dat_split)*top_int/100) 
    if (type_dat == "convolved")  {
      cmd <- paste("convolve_dca ",out_dca_file, out_sse_file,l_y,
                   l, int_size, xy_gauss, xy_gauss, coef, 
                   res_con_file, sep = " "
      )
      cmd_f <- paste("rm ",res_con_file,sep = " ")
      system(cmd_f)
      system(cmd)
      res_con <- read.table(res_con_file,header = F)
      res_mel <- melt(t(res_con))
      mat     <- cbind.data.frame(dat_split,power = res_mel$value)
      mm      <- mat  %>% filter(V3 >= min_ECs) %>% arrange(desc(power)) %>% slice(1:nb_rows) %>% arrange(power)
      if (save_file) {
        file_n <- paste("Convolved_ECs_",names[1],"_",names[2],"_prop_",top_int,".xlsx", sep = "")
        dat_f  <- data.frame(seq_a = mm[,1], res_a = seq_x[mm[,1]],seq_b = mm[,2], res_b = seq_y[mm[,2]], ECs = mm[,3], Convol_ECs = mm[,4])
        write.xlsx(dat_f,file = file_n, row.names = FALSE)
      }
      mm     <- mm  %>% select(V1,V2,power)
    }
    if (type_dat == "original")   {
      mm <- dat_split  %>% filter(V3 >= min_ECs) %>% arrange(desc(V3)) %>% slice(1:nb_rows) %>% arrange(V3)
      if (save_file) {
        file_n <- paste("ECs_",names[1],"_",names[2],"_prop_",top_int,".xlsx", sep = "")
        dat_f  <- data.frame(seq_a = mm$V1, res_a = seq_x[mm$V1],seq_b = mm$V2, res_b = seq_y[mm$V2], ECs = mm$V3, Convol_ECs = mm$V4)
        write.xlsx(dat_f,file = file_n, row.names = FALSE)
      }
    }
    colnames(mm) <- c('V1','V2','V3')
    
    cat("Total size: ",nrow(dat_split),"Number of contacts: ",nrow(mm),"\n")
    
    p <- ggplot() + theme_bw() 
    p     <- p + geom_point(data = mm, aes(x = V1, y = V2, colour = V3),shape = 16, size = 1.5) 
    p     <- p + scale_colour_gradientn(colours = c("blue","purple","red"))
    p     <- p +  labs(y = yla, x = xla)
    p     <- p + scale_x_continuous(limits = c(xmin - 5,xmax + 5), expand = c(0, 0)) +
                 scale_y_continuous(limits = c(ymin - 5,ymax + 5), expand = c(0, 0))
    if (1) {
      p <- p + theme(
        axis.title.x = element_text(vjust = 10,
                                    size = 25,
                                    family = "Courier",
                                    face   = "bold",
                                    color  = "black"
        ),
        axis.title.y = element_text(size = 25,
                                    family = "Courier",
                                    face  = "bold",
                                    color  = "black"
        ),
        axis.text.x = element_text(size = 20,
                                   family = "Courier",
                                   face  = "bold",
                                   color  = "black",
                                   angle = 00,
                                   vjust = .5),
        axis.text.y  = element_text(size = 20,
                                    family = "Courier",
                                    face  = "bold",
                                    color  = "black"
        ),
        plot.title   = element_text(size   = 15,
                                    face   = "bold",
                                    family = "Courier",
                                    color  = "blue",
                                    hjust  = 0.5,
                                    lineheight = 1.2
        ),
        legend.position = "none"
      ) 
    }
    ggplotly(p)
  },
  xmin       = mwSlider(1, l,   value = 1),
  xmax       = mwSlider(1, l,   value = l),
  ymin       = mwSlider(1, l_y, value = 1),
  ymax       = mwSlider(1, l_y, value = l_y),
  top_int    = mwSlider(1, 100, value = 10),
  min_ECs    = mwSlider(mi_EC, ma_EC, value = min_signal),
  save_file  = mwSelect(c(TRUE,FALSE), value = FALSE),
  xy_gauss   = mwSlider(0.0001,2,value = 0.01),
  coef       = mwSlider(0.001,1, value = 0.001),
  int_size   = mwSlider(1,35,    value = 21),
  type_dat   = mwSelect(c("convolved","original"), value = "convolved")
  )
}

### Does de concatenation (joining the 2 MSA to create the cMSA) 

compute_msa_and_align <- function(fasta1,      # first MSA
                                      fasta2,     # second MSA
                                      type_fast1 = NULL,# Any keyword to recognize the sequence
                                      type_fast2 = NULL, 
                                      one_msa = NULL, # set the order of the output cMSA
                                      ref_species, 
                                      out_file, 
                                      cut = 5
                                      ){
  
  # read 2 fasta files, concatenate them, do the alignment and return msa file
  # one_msa: will do the msa for just one_file
  # ref_species is the reference protein sequence
  ## Please set these directories
  dir_in  <- "Blast-data"
  dir_ou  <- "MSA-data"
  
  
  return_annot <- function(x) {
    f5 <-    x
    sp <- strsplit2(strsplit2(f5,"[[]")[2],"[]]")[1]
    sp
  }
  Annot_protein <- function(p){
    # p AA sequence
    name <- sapply(names(p),return_annot) # get the species
    sp   <- as.character(sapply(as.character(name), function(x) paste(strsplit2(x," "),collapse = "_"))) # collapse by _
    sp
  }
  
  give_name_from_var <- function(x,term=NULL){
    y <- strsplit2(x,"[.]")
    y <- ifelse(length(y) > 1, y[1],y)
    y <- strsplit2(y,"-")
    y <- ifelse(length(y) > 1, paste(y,collapse = "_"),y)
    if (!is.null(term)) y <- paste(y,term,sep = ".")
    y
  }
  
  collapse_proteins <- function(p1,p2) {
    p1 <- as.data.frame(p1)
    p2 <- as.data.frame(p2)
    res <- data.frame(matrix(NA,nrow = nrow(p1),ncol = 1))
    row.names(res) <- row.names(p1)
    for (t in 1:nrow(p1)){
      t1      <- as.character(p1[t,])
      nam     <- as.character(row.names(p1)[t])
      t2      <- as.character(p2[nam,])
      res[t,] <- paste(t1,t2,sep = "")
    }
    res  
  }
  
  write_msa <- function(fa,file,ref_species){
    if (file.exists(file))    file.remove(file)
    len <- fa$nb
    pos <- which(fa$nam == ref_species)
    write.fasta(toupper(fa$seq[pos]),ref_species ,file,"a" )
    ##write.fasta(fa$seq[pos],ref_species ,file,"a" )
    fa$nam <- fa$nam[-pos]
    fa$seq <- fa$seq[-pos]
    len    <- len - 1
    for (i in 1:len) {
      #write.fasta(fa$seq[i],fa$nam[i],file,"a" ) 
      write.fasta(toupper(fa$seq[i]),fa$nam[i],file,"a" ) 
    }
  }
  
  extract_type_only_protein <- function(p,type){
    tts <- sapply(names(p), function(x) length(grep(type,as.character(x))))
    tt2 <- p[which(names(p) %in% names(p)[which(tts == 1)]),]
    tt2
  }
  
  
  
  file1 <- paste(dir_in,fasta1,sep = "/")
  file2 <- paste(dir_in,fasta2,sep = "/")
  
  out1 <- give_name_from_var(fasta1)
  out2 <- give_name_from_var(fasta2)
  out_msa             <- paste(dir_ou,"/",out1,"_",out2,"_joined.fasta",sep = "");
  out_file1_file_2    <- paste(dir_ou,"/",out1,"_",out2,"_msa.fasta",sep = "");
  
  
  if (file.exists(out_file1_file_2))    file.remove(out_file1_file_2)
  if (file.exists(out_msa))    file.remove(out_msa)

    # join the 2 fasta files by species
  
  seq1 <- readAAStringSet(file1)
  seq2 <- readAAStringSet(file2)
  
  #make sure the sequences are from the same protein lineage
  if (!is.null(type_fast1)) seq1 <- extract_type_only_protein(seq1,type_fast1)
  if (!is.null(type_fast2)) seq2 <- extract_type_only_protein(seq2,type_fast2)    
  
  names(seq1) <- Annot_protein(seq1)
  names(seq2) <- Annot_protein(seq2)
  
  c_min_1 <- round(length(seq1[[ref_species]])*(1 - cut/100))
  c_max_1 <- round(length(seq1[[ref_species]])*(1 + cut/100))
  c_min_2 <- round(length(seq2[[ref_species]])*(1 - cut/100))
  c_max_2 <- round(length(seq2[[ref_species]])*(1 + cut/100))
  
  seq1_leng <- sapply(1:length(seq1), function(x) length(seq1[[x]]))
  seq2_leng <- sapply(1:length(seq2), function(x) length(seq2[[x]]))
  
  seq1  <- seq1[which(seq1_leng >= c_min_1 & seq1_leng <= c_max_1)]
  seq1  <- seq1[!duplicated(names(seq1))]
  seq2  <- seq2[which(seq2_leng >= c_min_2 & seq2_leng <= c_max_2)]
  seq2  <- seq2[!duplicated(names(seq2))]
  
  if (is.null(one_msa)) {
    
    if (as.numeric(length(seq1)) <= as.numeric(length(seq2))) { cm <-  names(seq1)[which(names(seq1) %in% names(seq2))]
    }else{cm <- names(seq2)[which(names(seq2) %in% names(seq1))] }
    
    seq1_2 <- seq1[which(names(seq1) %in% cm)]
    seq2_1 <- seq2[which(names(seq2) %in% cm)]
    seq_col <- collapse_proteins(seq1_2,seq2_1)
    
    for (i in 1:nrow(seq_col)){
      write.fasta(seq_col[i,], row.names(seq_col)[i],out_msa,"a" ) 
    }
  }
  else { 
    if (one_msa == "first"){
      seq_col <- seq1
      one_out <- out1
      print(seq_col)
    }else{
      seq_col <- seq2
      one_out <- out2
      
    }
    seq_col <- as.data.frame(seq_col)
    for (i in 1:nrow(seq_col)){
      write.fasta(seq_col[i,], row.names(seq_col)[i],out_msa,"a" ) 
    }
    out_file1_file_2      <- paste(dir_ou,"/",one_out,"_msa.fasta",sep = "")
  }    
  seq_mix <- readAAStringSet(out_msa)
  tt <- as.data.frame(names(seq_mix))
  row.names(tt) <- seq(1,length(seq_mix))
  muscle(seq_mix,quiet = T,  fastaout = out_file1_file_2);
  aln <- read.alignment(file = out_file1_file_2,format = 'fasta',forceToLower = FALSE);
  aln$nam <- as.character(tt[aln$nam,])
  write_msa(aln, out_file1_file_2, ref_species = ref_species)
}


