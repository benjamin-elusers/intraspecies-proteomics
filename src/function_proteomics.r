library(tidyverse)
library(rio)
library(openxlsx)
# WORKING DIRECTORY is current document path
if(interactive()){
  WD = dirname( .rs.api.getActiveDocumentContext()$path )
  setwd(WD)
}

load.annotation = function(){
  # Preloaded uniprot data can be generated in 5min with:
  #   uni = load.uniprot.features(tax="559292",refdb="UNIPROTKB")
  #   sgd = load.sgd.features()
  
  uni_feat = read_rds(here('data','uniprot-features.rds')) %>%
    dplyr::select(-c(REVIEWED,COMMENTS,SUBLOC))
  sgd_desc = read_rds(here('data','uniprot-sgd-annotation.rds'))
  biofunc = load.vanleeuwen2016.data(single_orf=T)
  enog_annot = eggnog_annotations_species(4891,4932)
  
  annotation = full_join(sgd_desc,uni_feat,by=c("SGD","UNIPROT"='UNIPROTKB')) %>%
    full_join(biofunc,by='ORF')  %>%
    full_join(enog_annot, by =c('ORF'='string')) %>%
    filter(!is.na(ORF) | is.na(UNIPROT)) %>%
    mutate(GENENAME = ifelse(is.na(GENENAME),ORF,GENENAME)) %>%
    mutate(letter = fct_drop(letter)) %>%
    relocate(SGD,GENENAME,ORF,UNIPROT,PNAME,
             L,FAMILIES,FUNCTION,ROLE,BIOPROCESS_all,enog_annot,
             LOC,COMPLEX,ORTHO,OTHER,KEYWORDS,
             EXISTENCE,SCORE)
  
  
  return(annotation)
}

search_zip = function(zipfile,target){
  return( unzip(zipfile, list = TRUE)$Name %>% 
            str_subset(pattern=glue::glue(".+{target}+")) )
}

sys_unzip = function(zipfile,target,outdir){
  targets = search_zip(zipfile,target)
  print(targets)
  cmd = sprintf("unzip '%s' '*%s*' -d '%s' ",zipfile,target,outdir)
  cat(cmd,'\n')
  system2(command = cmd)
  return(file.path(outdir,targets))
}

regex.sample = function(I="^int_",S="(BTT|CQC)",R="([12])", M="\\-(LL|SL|SS)", partial=T){
  if(partial){
    return( sprintf("%s?(%s?%s?%s?)",I,S,R,M) )
  }
  return(sprintf("%s(%s%s%s)",I,S,R,M))
}

R1 = regex.sample(R="1")
R2 = regex.sample(R="2")
LL = regex.sample(M="-LL")
SL = regex.sample(M="-SL")
SS = regex.sample(M="-SS")
BTT = regex.sample(S="BTT")
CQC = regex.sample(S="CQC")
INTENSITIES = regex.sample(I="int_",S="",R="",M="",partial = F)

get_group = function(sample_names,sep="_",grp_num=1){
  ngroup = unique(str_count(sample_names,sep))
  groups = str_split_fixed(sample_names,sep,n=max(grp_num)+1)
  if(length(ngroup)>1){
    warning('groups are not consistent within sample names!')
  }else if(length(grp_num)>1){ 
    warning('Combine groups to color!')
    group = apply(groups[,unique(grp_num),drop=F],1, paste0, collapse="_")
  }else{
    if(grp_num>0){
      group = apply(groups[,(grp_num),drop=F],1, paste0, collapse="_")
    }else{
      group = sample_names 
    }
  }
  return(group)
}

transpose = function(df){ 
  t_df = as.matrix(df) %>% as_tibble() %>%
    rownames_to_column() %>%
    pivot_longer(-rowname, 'variable', 'value') %>%
    pivot_wider(names_from = variable, rowname)
  return(t_df)
}

#### 0 Read data ----------------------------------------------------------
read_proteomics_results = function(datain=ms.resfile,zero.to.na=T){
  MS = import(datain)
  cat(sprintf("Total number of proteins hits: %s\n",nrow(MS)))
  
  col.ratios = str_subset( colnames(MS), "Ratio") %>% str_replace("Ratio ","") %>% str_replace_all("/ ","/")
  col.samples = str_subset( colnames(MS), "Intensity") %>% str_replace("Intensity ","int_")
  cat("Sample names:\n",sprintf("%2s. %s\n",seq_along(col.samples),col.samples))
  colnames(MS) = c('uniprot','major_uniprot','protein_names','gene_names','fasta_header','pep','upep',
                   col.ratios, col.samples,
                   'seqcov_perc','MW_kDa','qv','score')
  if(zero.to.na){
    # Convert 0-intensities to NA
    zeros = MS[,col.samples] == 0
    MS[,col.samples][zeros] = NA
  }
  return(MS)
}

read_maxquant = function(datain,zero.to.na=T, int_type='LFQ',pep_type='Peptides', sample_pattern="(wt|BTT[12]-[SL]L|(25|58))"){
  MAXQ = import(datain,showProgress=T)
  cat(sprintf("Total number of proteins hits: %s\n",nrow(MAXQ)))
  COLS = colnames(MAXQ)
  KEYS = c("Protein IDs","Majority protein IDs","Gene names","Number of proteins","Peptides","Razor + unique peptides","Unique peptides")
  INFO = c("Sequence coverage [%]","Unique sequence coverage [%]","Sequence length","Sequence lengths","Mol. weight [kDa]","Score")
  
  # Intensities are the sums of all individual peptide intensities belonging to a particular protein group. 
  # Unique and razor peptide intensities are used as default.
  intensity_type = match.arg(int_type,choices = c('LFQ','iBAQ','Intensity'), several.ok = F)
  peptide_type = match.arg(pep_type,choices = c('Peptides','Unique peptides','Razor + unique peptides'), several.ok = F)
  
  # LFQ: 
  # iBAQ: Σ intensity/#theoretical peptides
  sample_names =  str_subset( COLS, "^Identification type") %>% str_replace('Identification type ' ,"")
  col_int = str_subset( COLS, paste0("^",int_type))
  col_pep = str_subset( COLS, paste0(peptide_type," "))

  nsel = sum( str_count(sample_names,sample_pattern) )
  selected = str_detect(sample_names,sample_pattern)
  sample_names[selected] = paste0(sample_names[selected], " (x)")
  cat("Sample names:\n",
      ifelse(nsel>0,sprintf("...selecting %s samples (x)...\n",nsel),""),
      sprintf("%2s. %s\n",seq_along(sample_names),sample_names))
  
  if(nsel>0){
    col_int = str_subset(col_int,sample_pattern)
    col_pep = str_subset(col_pep,sample_pattern)
  }
#  colnames(MS) = c('uniprot','major_uniprot','protein_names','gene_names','fasta_header','pep','upep',
#                   col.ratios, col.samples,
#                   'seqcov_perc','MW_kDa','qv','score')
  if(zero.to.na){
    # Convert 0-intensities to NA
    zeros = MAXQ[,col_int] == 0
    MAXQ[,col_int][zeros] = NA
  }
  
  maxq = MAXQ %>% 
          dplyr::select(all_of(c(KEYS,col_int,col_pep))) %>% 
          janitor::clean_names(case="snake",parsing_option=1) %>% 
          rename_with(.fn = tolower)
  
  return(maxq)
}
#ms.resfile = here::here('data','proteomics-results-INPCM-EL_14975_061221.xlsx')
#ms0=read_proteomics_results(ms.resfile)

## 1 Filter data ----------------------------------------------------------
filter_hits = function(MS0=ms0,id='majority_protein_i_ds',np=2){
  # Mark contaminants
  MS0$is_contaminated = str_detect(MS0[[id]],"CON")
  # Mark reverse sequenced
  MS0$is_reversed= str_detect(MS0[[id]],"REV")
  
  # Mark proteins with at least 2 unique peptides
  MS0$has_upep = !is.na(MS0$unique_peptides) & MS0$unique_peptides >= np
  # Mark proteins that are not distinguished
  MS0$has_many = str_detect(MS0[[id]],";")
  MS1 = MS0 %>% 
    dplyr::filter(!is_contaminated & !is_reversed & has_upep & !has_many)
  cat("Discarding problematic hits...\n")
  cat(sprintf("* %5s = less than 2 unique peptides\n",sum(!MS0$has_upep)))
  cat(sprintf("* %5s = contaminated hits\n",sum(MS0$is_contaminated)))
  cat(sprintf("* %5s = reversed sequences\n",sum(MS0$is_reversed)))
  cat(sprintf("* %5s = multi-protein hits\n",sum(MS0$has_many)))
  nelim = nrow(MS0)-nrow(MS1)
  cat("-----------------------------------------\n")
  cat(sprintf(" -> %-5s hits eliminated\n",nelim))
  cat(sprintf(" => %-5s remaining hits for analysis\n",nrow(MS1)))
  
  return(MS1)
}
#ms1=filter_hits(ms0)

get_intensities = function(MS){
  intensities = MS %>% dplyr::select(uniprot,matches(INTENSITIES)) %>% column_to_rownames("uniprot")
  return(intensities)
}
#int_all = get_intensities(ms1)

get_long_intensities = function(intensities,int.col='lfq_intensity_'){

  long_int_all = intensities %>% 
    pivot_longer(cols=-uniprot,names_to=c('strain','bio','tech'), 
                 names_prefix = int.col, names_sep = '_',values_to = 'int') %>% 
    group_by(uniprot,strain) %>% 
    # Calculate number of missing values across replicates
    mutate(log10_int = log10(int), n_na_rep=sum(is.na(int)),
           ratio_na_rep=mean(is.na(int)) ) %>%
    # Calculate number of strains with missing value for each hit
    group_by(uniprot) %>%
    mutate(n_na_strains = n_distinct(unique(strain[is.na(int)])),
           na_strains = paste0(unique(strain[is.na(int)]),collapse="|")) %>% 
    # Calculate number of hits not detected across strains
    group_by(n_na_strains) %>% 
    mutate( nprot = n_distinct(uniprot),
            has_missing_strains= n_na_rep == 0 & n_na_strains >= 0,
            nprot_miss=n_distinct(uniprot[has_missing_strains])) 
  return(long_int_all)
}

get_int_byrep = function(long_intensities, fun=mean_){
  df=pivot_wider(long_intensities,
              id_cols=uniprot,
              names_from = c('media','strain'), names_glue = "{strain}-{media}", 
              values_from = starts_with("int"), values_fn=fun )
  return(df)
}

get_int_bystrain = function(long_intensities, fun=mean_){
  df=pivot_wider(long_intensities,
                 id_cols=uniprot,
                 names_from = c('media','rep'), names_glue = "{media}-{rep}", 
                 values_from = starts_with("int"), values_fn=fun )
  return(df)
}

get_int_bymedia = function(long_intensities, fun=mean_){
  df=pivot_wider(long_intensities,
                 id_cols=uniprot,
                 names_from = c('strain','rep'), names_glue = "{strain}-{rep}", 
                 values_from = starts_with("int"), values_fn=fun )
  return(df)
}

get_int_bysample = function(long_intensities, fun=mean_){
  df=pivot_wider(long_intensities,
                 id_cols=uniprot,
                 names_from = c('sample'), 
                 values_from = starts_with("int"), values_fn=fun )
  return(df)
}

count_NA = function(INT,groups){
  cat("Counting NA in each sample/groups...\n")
  # Count NA for each subgroup (media,replicates,strains)
  # for replicate
  INT$na_R1 = INT %>% dplyr::select(matches("1")) %>% is.na() %>% rowSums()
  INT$na_R2 = INT %>% dplyr::select(matches("2")) %>% is.na() %>% rowSums()
  INT$na_rep = INT$na_R1 + INT$na_R2
  # for media
  INT$na_LL = INT %>% dplyr::select(matches("-LL")) %>% is.na() %>% rowSums()
  INT$na_SL = INT %>% dplyr::select(matches("-SL")) %>% is.na() %>% rowSums()
  INT$na_SS = INT %>% dplyr::select(matches("-SS")) %>% is.na() %>% rowSums()
  INT$na_media = INT$na_LL + INT$na_SL + INT$na_SS
  # for strain
  INT$na_BTT = INT %>% dplyr::select(matches("BTT")) %>% is.na() %>% rowSums()
  INT$na_CQC = INT %>% dplyr::select(matches("CQC")) %>% is.na() %>% rowSums()
  INT$na_strain = INT$na_BTT + INT$na_CQC
  # for all samples
  INT$na_all = INT %>% dplyr::select(matches(INTENSITIES)) %>% is.na() %>% rowSums()
  return(INT)
}

remove_NA_rows = function(INT,ns=3,nm=2){

  INT2 = count_NA(INT) %>% dplyr::filter(na_BTT<ns & na_CQC<ns & (na_LL<nm | na_SL<nm) & na_SS<nm)
  
  cat("Removing hits with NA in multiple samples/groups...\n")
  cat("Max. #NAs per strain (BTT/CQC) = ",ns,"\n")
  cat("Max. #NAs per media (LL/SL/SS) = ",nm,"\n")
  
  cat(sprintf(" => %s eliminated hits\n",nrow(INT)-nrow(INT2)))
  cat(sprintf(" -> %s remaining hits\n",nrow(INT2)))
  return(INT2)
}
#int=remove_NA_rows(int_all,ns=3,nm=2)

## 2 Normalize intensities --------------------------------------------------------------
normalize_intensities = function(int.raw){
  cat("Normalize intensities to the median...\n")
  # transform positive intensities to fold-change (log2)
  int.log2 = log2(int.raw+1) %>% as.matrix
  int.log2[!is.finite(int.log2)] = NA

  # center each sample on median
  md.int = apply(int.log2,2,hablar::median_)
  int.norm = sweep(int.log2,MARGIN=2,STATS=md.int,FUN='-')
  min.int = min_(int.norm)
  int.pos = int.norm+abs(min.int)
  
  return(int.pos)
}
#int.norm=normalize_intensities(int)

## 3 Correlation of intensities ----------------------------------------------------------
compute_samples_correlation = function(datain=int.norm,as.df=F){
  #install.packages('corrr')
  cat("Compute pairwise samples correlation (Spearman)...\n")
  library(corrr)
  COR = cor(datain,met = 'spearman',use='pairwise.complete') %>% as_tibble() %>% as.matrix()
  rownames(COR)=colnames(COR)
  df.COR = correlate(datain,met = 'spearman',use = 'pairwise.complete') %>% stretch()
  if(as.df){ return(df.COR) }
  return(COR)
}
#COR=compute_samples_correlation(int.norm)

#install.packages('ggcorrplot')
#library(ggcorrplot)
# ggcorrplot(COR,lab = T,show.diag = F, outline.color = NA,hc.order = T, hc.method = 'ward.D2') + 
#   scale_fill_gradient2(limit = c(0.9,1), low = "blue", high =  "red", mid = "white", midpoint = 0.95)

# txtsize <- par('din')[2] / 1.5
# p = ggplot(df.COR, aes(x, y, fill=r)) +
#   geom_tile() + 
#   xlab(NULL) + ylab(NULL) +
#   theme_classic(base_size = 14) +
#   scale_fill_gradient2(limit = c(0.9,1), low = "blue", high =  "red", mid = "yellow", midpoint = 0.95, na.value = 'white',name='spearman') +
#   geom_text(aes(label=round(r,2)), size=txtsize * 0.8, color="grey9") +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, size = 8),
#         axis.text.y = element_text(size = 8))  
# p

draw_heatmap_samples = function(mcor,df.group,col.group,k=2){

  #install.packages('ggplotify')
  library(ggplotify)
  library(cowplot)
  #graphics.off()
  p <- as.ggplot(
        pheatmap::pheatmap(mcor,display_numbers = T,
                           annotation_row = df.group, annotation_col=df.group, annotation_colors =  col.group,
                           cutree_rows = k,gaps_row = k)
  )
  #save_plot(p,filename=here::here('output','correlation-heatmap-norm.intensities.pdf'),base_aspect_ratio = 0.8)
  # pheatmap::pheatmap(mcor,display_numbers = T, cutree_rows = 2,gaps_row = 2, 
  #                    annotation_row = df.anno, annotation_col=df.anno, annotation_colors =  color.anno,
  #                    width = 11,height=10, border_color = NA, #cellwidth = 4, cellheight = 4,
  #                    filename =here::here('output','correlation-heatmap-norm.intensities.pdf'))
  # pheatmap::pheatmap(mcor,display_numbers = T, cutree_rows = 2,gaps_row = 2, fontsize = 14,
  #                    annotation_row = df.anno, annotation_col=df.anno, annotation_colors =  color.anno,
  #                    width = 11,height=10, border_color = NA, #cellwidth = 4, cellheight = 4,
  #                    filename =here::here('output','correlation-heatmap-norm.intensities.png'))
  return(p)
}

# table( rowSums(is.na(int.norm)) )
# apply(int.norm,2,function(x){ sum(!is.na(x)) })
# 
# pheatmap::pheatmap(int.norm,annotation_col=df.anno, annotation_colors =  color.anno)
# 

# 4 Scatterplots ----------------------------------------------------------
draw_scatterplots = function(datain=ms2){
  library(GGally)
  #graphics.off()
  df = datain %>% as_tibble() %>% as.data.frame()
  p_ <- GGally::print_if_interactive
  p2 = ggpairs(df,
               upper = list(continuous=wrap(ggally_cor, display_grid = FALSE, 
                                            title="spearman", method='spearman', exact=F, use='pairwise.complete',digits=2,size = 3)),
               lower = list(continuous = wrap("smooth",col='dodgerblue',size=0.5,alpha=0.5,shape=19))
               ) + theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())
  p_(p2)
  #save_plot(p2,filename=here::here('output','scatmat-norm.intensities.pdf'),base_aspect_ratio=1,base_height = 12,base_width = 12)
  #save_plot(p2,filename=here::here('output','scatmat-norm.intensities.png'),base_aspect_ratio=1,base_height = 12,base_width = 12)
  return(p2)
}

compare_normalization = function(rawexp, normalization = c()){
  
  library(limma)
  limma::normalizeCyclicLoess()
  limma::normalizeMedianAbsValues()

  if(all){
    n = ncol(rawexp)/2
    s1 = colnames(rawexp)
    s2 = colnames(rawexp)
  }else{
    stopifnot("requires even number of samples!"=is.even(length(paired)))
    n = length(paired)/2
    s1 = odds(paired) 
    s2 = evens(paired)
  }
  par(pch=19,cex.lab=1.4,cex.axis=1.2)
  if(!onebyone){ graphics::layout(n) }

  for( i in 1:n){
    exp_unit= "(log2)-intensity"
    plot( rawexp[,s1[i]], rawexp[,s2[i]], xlab=paste(s1[i],exp_unit) ,ylab=paste(s2[i],exp_unit), col=rgb(0,0,0,0.5), cex=0.6)
    points(normexp[,s1[i]], normexp[,s2[i]],col=rgb(1,0,0,0.5), cex=0.6)
    legend('topleft',title='protein expression',legend = c('raw','norm'), 
           col = 1:2, pch = 19, border = NA, bty='n')
    if(n>1 & onebyone){ invisible(readline(prompt="Press [enter] to continue")) }
  }
}

# 5 Principal component analysis ------------------------------------------
make_pca = function(x,with_labels=F,col_by_group=2){
  pca.samples = prcomp(t(x))
  pov <- 100*pca.samples$sdev^2/sum(pca.samples$sdev^2)
  pca.data <- as.data.frame(pca.samples$x)
  ngroup = unique(str_count(rownames(pca.data),"_"))
  pca.data$samples = rownames(pca.data) #sample names
  pca.data$group =get_group(sample_names = rownames(pca.data),sep="_",grp_num = col_by_group)

  p <- ggplot(pca.data, aes(PC1,PC2)) +
    geom_point(pch = 21, size = 4, colour = "black", alpha = 0.7, aes(fill = group),show.legend = !with_labels) +
    theme_classic() +
    ylab(sprintf("Principle component 2 (%.1f%%)",pov[2])) +
    xlab(sprintf("Principle component 1 (%.1f%%)",pov[1])) + 
    labs(title = sprintf('%% variance = %.1f',sum(pov[1:2])))
  
    if(with_labels){
      p = p + geom_label_repel(aes(label=samples,fill=group), size=2.5,max.overlaps = 20,show.legend = F)
    }
  return(p)
}

# 6 Statistical difference -----------------------------------------------------
library(limma)
pairwise_condition = function(conditions = all_strains){
  cond_pair = expand.grid(S1=conditions,S2=conditions)
  true.pairs <- t(apply(cond_pair[,1:2], 1, sort))
  dup.pairs <- duplicated(true.pairs)
  unique.pairs <- true.pairs[!dup.pairs & true.pairs[,1]!=true.pairs[,2],] %>% 
    as_tibble %>% rename(cond1=V1,cond2=V2)
  return(unique.pairs)
}

compare_conditions = function(input, id_col = "UNIPROT",
                              comparison = pairwise_condition(all_strains),
                              col_group=1){
  
  if( id_col %in% colnames(input)){
    processed_data = na.omit(input) %>% column_to_rownames('UNIPROT')
  }else{
    processed_data = na.omit(input)
  }
  
  samples=colnames(processed_data)
 
  GROUP = get_group(sample_names = colnames(processed_data),sep="_",grp_num = col_group)
  
  f.df <- factor(GROUP)
  design <- model.matrix(~0+f.df)
  colnames(design) <- levels(f.df)
  fit <- lmFit(processed_data, design)

  cont.matrix <- makeContrasts(contrasts = apply(comparison,1,paste0,collapse='-'), levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit3 <- eBayes(fit2)
  return(fit3)
}

get_volcano_data = function(input_data=int_norm, min_lfc=2, min_pval=0.01, which=c('both','up','down'), topn = 20){
  datList <- list()
  combination <- pairwise_condition()
  stat_comb = apply(combination,1,paste0,collapse='-')
  
  for (i in 1:length(stat_comb)) {
    fit2 <- compare_conditions(input_data, comparison=combination, col_group = 1)
    
    d.out <- data.frame(ID = names(fit2$coefficients[,i]),
                        pValue = fit2$p.value[,i],
                        qValue = p.adjust(fit2$p.value[,i], "fdr"),
                        EffectSize = fit2$coefficients[,i],
                        comparison = stat_comb[i])
    d.out <- mutate(d.out, 
                    sig = ifelse(d.out$EffectSize > min_lfc & round(d.out$qValue, 3) < min_pval, "Upregulated",
                                 ifelse(d.out$EffectSize < (min_lfc * -1) & round(d.out$qValue, 3) < min_pval, "Downregulated", "Non significant")))
    
    if (which == 'up') {
        d.up = d.out %>% 
          filter(sig == "Upregulated") %>%
          arrange(desc(EffectSize)) %>%
          slice(1:topn)
    } else if (which == 'down') {
        d.down = d.out %>% 
          filter(sig == "Downregulated") %>%
          arrange(EffectSize) %>%
          slice(1:topn)
    } else if (which == 'both') {
        d.both = d.out %>% 
          filter(sig != "Non significant") %>%
          arrange(desc(EffectSize)) %>%
          group_by(sig) %>%
          slice(1:topn*2)
    }

    rownames(d.out) <- d.out$ID
    datName <- stat_comb[i]
    datList[[datName]] = d.out
  } #for loop close
  return(datList)
}

get_dfe = function(INPUT=int_norm, MIN_LFC=2, MIN_PVAL=0.01, WHICH='both', TOPN = 20){
  dfe <- get_volcano_data(input_data=INPUT, min_lfc=MIN_LFC, min_pval=MIN_PVAL, WHICH, topn = TOPN) %>% 
            bind_rows %>% 
            dplyr::filter(sig!="Non significant")
  return(dfe)
}
volcPlot = function(INPUT=int_norm, MIN_LFC=2, MIN_PVAL=0.01, WHICH='both', TOPN = 20){
    
  plotList <- list()
  dlist <- get_volcano_data(input_data=INPUT, min_lfc=MIN_LFC, min_pval=MIN_PVAL, WHICH, topn = TOPN)
  for( ncomp in 1:length(dlist) ){
    comparison_current = names(dlist)[ncomp]
    d=dlist[[ncomp]]
    d.down = sum(round(d$qValue, 3) < MIN_PVAL & d$EffectSize < MIN_LFC*-1) 
    d.up = sum(round(d$qValue, 3) < MIN_PVAL & d$EffectSize > MIN_LFC)
    p <- ggplot(d, aes(x=EffectSize, y=-log10(qValue), fill = sig)) +
      xlab("log2 fold change") + ylab("-log10 q-value") + labs(fill = NULL) +
      ggtitle(label = comparison_current ) +
      theme_classic(base_size = 14) +
      geom_point(pch = 21, colour = "black", alpha = 0.5, size = 1.5) +
      scale_fill_manual(values=c("Non significant" = 'gray', "Downregulated" = 'red', "Upregulated" = 'blue')) +
      scale_color_manual(values=c("Non significant" = 'gray', "Downregulated" = 'red', "Upregulated" = 'blue')) +
      theme(legend.position = 'top',
            plot.title = element_text(face  = 1, hjust =0, size = 18),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
    
    if (WHICH == "up" || WHICH == "down") {
      p <- p + geom_text_repel(data=subset(d,sig!='Non significant'),aes(label = ID,col=sig), show.legend = FALSE)
    }
    
    if(WHICH == "both"){ 
      p <- p + geom_text_repel(data=subset(d,sig!='Non significant'),aes(label = ID,col=sig), show.legend = FALSE)# +
       # geom_text_repel(aes(label = ID), show.legend = FALSE)
    }
    
    p <- p +  
        geom_vline(aes(xintercept = (MIN_LFC*-1)), lty = 'dashed', colour = 'red',lwd=0.5) +
        geom_vline(aes(xintercept = (MIN_LFC)), lty = 'dashed', colour = 'blue',lwd=0.5) +
        geom_hline(aes(yintercept = -log10(MIN_PVAL)), lty = 'dashed', colour = 'black',lwd=0.5) 
    plot(p)
    plotList[[comparison_current]] = p
  }
  # if (input$volcFeatures == "Counts") {
  #   p <- p + geom_text(aes(x = input$volcXdown, y= input$volcYdown, label=d.down)) +
  #     geom_text(aes(x = input$volcXup, y= input$volcYup, label=d.up))
  # }
  # 
  # if (input$volcFeatures == "Both") {
  #   p <- p + geom_vline(aes(xintercept = (input$UserFCCutoff*-1)),
  #                       lty = input$volcLinesType, 
  #                       colour = input$volcDown,
  #                       lwd=input$volcLinesLWD) +
  #     geom_vline(aes(xintercept = input$UserFCCutoff), 
  #                lty = input$volcLinesType, colour =input$volcUp, 
  #                lwd=input$volcLinesLWD) + 
  #     geom_text(aes(x = input$volcXdown, 
  #                   y=input$volcYdown, 
  #                   label=d.down)) +
  #     geom_text(aes(x = input$volcXup, 
  #                   y= input$volcYup, 
  #                   label=d.up))
  # } 
  return(plotList)
}
  
