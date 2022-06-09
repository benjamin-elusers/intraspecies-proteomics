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

exp2matrix = function(df_proteomics, int_col = 'intensity_', id_col = 'uniprot'){
  
  col_ids = df_proteomics %>% dplyr::select(where(is.character))
  id = match.arg(id_col, col_ids, several.ok = F)
  
  intensities = df_proteomics %>% dplyr::select(starts_with(int_col),id) 
  
}

# 
# transpose = function(df){ 
#   t_df = as.matrix(df) %>% as_tibble() %>%
#     rownames_to_column() %>%
#     pivot_longer(-rowname, 'variable', 'value') %>%
#     pivot_wider(names_from = c(variable,value), rowname)
#   return(t_df)
# }

#### 0 Read data ----------------------------------------------------------
read_proteomics_results = function(datain=ms.resfile,zero.to.na=T){
  library(rio)
  library(tidyverse)
  MS = rio::import(datain)
  cat(sprintf("Total number of proteins hits: %s\n",nrow(MS)))
  
  col.ratios = stringr::str_subset( colnames(MS), "Ratio") %>% str_replace("Ratio ","") %>% str_replace_all("/ ","/")
  col.samples = stringr::str_subset( colnames(MS), "Intensity") %>% str_replace("Intensity ","int_")
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

read_maxquant = function(datain,zero.to.na=T, int_type='LFQ',pep_type='Peptides',
                         sample_pattern=""){ #(wt|BTT[12]-[SL]L|(25|58))
  MAXQ = rio::import(datain,showProgress=T)
  cat(sprintf("Total number of proteins hits: %s\n",nrow(MAXQ)))
  COLS = colnames(MAXQ)
  KEYS = c("Protein IDs","Majority protein IDs","Gene names","Number of proteins","Peptides","Razor + unique peptides","Unique peptides")
  INFO = c("Sequence coverage [%]","Unique sequence coverage [%]","Sequence length","Sequence lengths","Mol. weight [kDa]","Score")
  
  # Intensities are the sums of all individual peptide intensities belonging to a particular protein group. 
  # Unique and razor peptide intensities are used as default.
  intensity_type = match.arg(int_type,choices = c('LFQ','iBAQ','Intensity'), several.ok = T)
  peptide_type = match.arg(pep_type,choices = c('Peptides','Unique peptides','Razor + unique peptides'), several.ok = F)
  
  # LFQ: 
  # iBAQ: Σ intensity/#theoretical peptides
  sample_names =  stringr::str_subset( COLS, "^Identification type") %>% str_replace('Identification type ' ,"")
  col_int = paste0(int_type," ") %>% purrr::map(.f = stringr::str_subset, string=COLS, negate=F) %>% unlist
  col_pep = stringr::str_subset( COLS, paste0(peptide_type," "))

  nsel = sum( stringr::str_count(sample_names,sample_pattern) )
  
  if(nsel>0){
    col_int = stringr::str_subset(col_int, sample_pattern)
    col_pep = stringr::str_subset(col_pep, sample_pattern)
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
          dplyr::rename_with(.fn = tolower)
  
  return(maxq)
}
#ms.resfile = here::here('data','proteomics-results-INPCM-EL_14975_061221.xlsx')
#ms0=read_proteomics_results(ms.resfile)

## 1 Filter data ----------------------------------------------------------
filter_hits = function(MS0=ms0,id='majority_protein_i_ds',np=2,verbose=T){
  # Mark contaminants
  MS0$is_contaminated = stringr::str_detect(MS0[[id]],"CON")
  # Mark reverse sequenced
  MS0$is_reversed= stringr::str_detect(MS0[[id]],"REV")
  
  # Mark proteins with at least 2 unique peptides
  MS0$has_upep = !is.na(MS0$unique_peptides) & MS0$unique_peptides >= np
  # Mark proteins that are not distinguished
  MS0$has_many = stringr::str_detect(MS0[[id]],";")
  MS1 = MS0 %>% 
    dplyr::filter(!is_contaminated & !is_reversed & has_upep & !has_many)
  nelim = nrow(MS0)-nrow(MS1)
  
  if(verbose){
    cat("Discarding problematic hits...\n")
    cat(sprintf("* %5s = less than 2 unique peptides\n",sum(!MS0$has_upep)))
    cat(sprintf("* %5s = contaminated hits\n",sum(MS0$is_contaminated)))
    cat(sprintf("* %5s = reversed sequences\n",sum(MS0$is_reversed)))
    cat(sprintf("* %5s = multi-protein hits\n",sum(MS0$has_many)))
    cat("-----------------------------------------\n")
    cat(sprintf(" -> %-5s hits eliminated\n",nelim))
    cat(sprintf(" => %-5s remaining hits for analysis\n",nrow(MS1)))
  }
  return(MS1)
}
#ms1=filter_hits(ms0)

get_intensities = function(MS, regex_int, col_id, uppercase=T, remove_prefix=T){
  
  intensities = MS %>% column_to_rownames({{col_id}}) %>% dplyr::select(matches(regex_int))
  if(remove_prefix){ intensities = intensities %>% rename_with(.fn = str_replace, pattern=regex_int, replacement='') }
  if(uppercase){ intensities = intensities %>% rename_with(.fn = str_to_upper) }
  return(intensities)
}
#int_all = get_intensities(ms1)

make_intensities_long = function(MS, regex_int, col_id, uppercase=T, remove_prefix=T){
  
  int_wide = get_intensities(MS, regex_int, col_id, uppercase, remove_prefix) %>% 
             rownames_to_column('id')
   
  int_prefix=hutils::longest_prefix(colnames(int_wide))
  int_long = pivot_longer(int_wide, cols=-id, 
                          names_to=c('sample'), names_prefix = int_prefix, values_to = 'int') %>% 
    mutate(log10_int = log10(int), log2_int = log2(int)) %>%
    # Calculate number of missing values per id
    group_by(id) %>% mutate( n_na_sample=sum(is.na(int)), f_na_sample=mean(is.na(int)) )
  
  return(int_long)
}
    

get_long_intensities = function(intensities,int.col='lfq_intensity_',use_log10=T){
  
  long_int_all = intensities %>% 
    pivot_longer(cols=-uniprot,names_to=c('strain','bio','tech'), 
                 names_prefix = int.col, names_sep = '_',values_to = 'int') %>% 
    group_by(uniprot,strain) %>%
    # Calculate number of missing values across replicates
    mutate(log10_int = log10(int), log2_int = log2(int), n_na_rep=sum(is.na(int)),
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


center_intensities = function(int.raw, center='median', tolog2=T ){

  if( tolog2 ){
    cat("Normalize log2-transformed intensities by the samples ",center,"...\n")
    # Transform positive intensities to fold-change (log2)
    INT = log2(int.raw+1) %>% as.matrix
    INT[!is.finite(INT)] = NA
  }else{
    cat("Normalize raw intensities to the samples ",center,"...\n")
    INT = int.raw  
  }

  # center each sample on median
  if(center == 'median'){  
    center.int = apply(INT,2,hablar::median_)
  }else if(center == 'mean'){  
    center.int = apply(INT,2,hablar::mean_)
  }
  
  int.norm = sweep(INT,MARGIN=2,STATS=center.int,FUN='-')
  min.int = min_(int.norm)
  int.pos = int.norm+abs(min.int)
  
  return(int.pos)
}
#int.norm=center_intensities(int, tolog2=T, center='median')

normalize_intensities = function(int,design=df.group){

  library(NormalyzerDE)
  m.int = int %>% mutate(across(everything(),as.numeric)) %>% as.matrix() #%>% remove_rownames() 
  experiment <- SummarizedExperiment::SummarizedExperiment(
    assays=list(raw = m.int ),
    rowData = rownames(int),
    colData = DataFrame( design, row.names = NULL),
    metadata = list(sample='sample',group='strain'),
    checkDimnames=F
  )  
  #source("https://raw.githubusercontent.com/ByrmLab/proteiNorm/master/normFunctions.R")
  
  norm <- getVerifiedNormalyzerObject('normalized_data', experiment)
  norm_res <- normMethods(norm)
  norm_res_eval <- analyzeNormalizations(norm_res)
  
  #generatePlots(normperf)
  #exp_normloess =SummarizedExperiment::SummarizedExperiment(assays = list(norm=NORM@normalizations$CycLoess), colData=df.group, metadata=list(sample='sample',group='strain'))
  #input <- getVerifiedNormalyzerObject('loess_normalized', summarizedExp = exp_normloess)
  # nst <- NormalyzerStatistics(experiment, logTrans=FALSE)
  # combination <- pairwise_condition(conditions = df.group$strain) %>% mutate( comparison=paste0(cond1,"-",cond2))
  # nst <- calculateContrasts(nst, combination$comparison, condCol="strain", leastRepCount=2)
  # annotDf <- generateAnnotatedMatrix(nst)
  # generateStatsReport(nst,jobDir = '.',jobName = 'loess_strain')
  return(norm_res_eval)
}

sum_intensities = function(exp_mat,int_prefix=''){
  tot_int = pivot_longer(exp_mat, where(is.numeric),
                         names_to = c('sample','strain','bio','tech','day'),
                         names_pattern = "(([^_]+)_([^_]+)_([^_]+)_([^_]+))",
                         names_prefix = int_prefix,
                         values_to='int') %>% 
    group_by(sample,strain,bio,tech,day) %>% 
    summarize(total_intensity=sum_(as.numeric(int))) 
  return(tot_int)
}

draw_barplot_sumint = function(tot_int, y_label = 'Total intensity (raw)' ){
  
  bp= ggplot( tot_int, aes(x=sample,y=total_intensity,fill=strain,label=sample)) + 
    geom_col() + 
    geom_text(angle=90,check_overlap = T,hjust=-0.1,y=0,size=2) + 
    ylab('Total intensity (lfq)') +
    theme(axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position='none')
  return(bp)
}

# 3 Variations of intensities --------------------------------------------------
calculate_cv = function(intensity_long, by_bio=T){
  
  CV = intensity_long %>%
    drop_na() %>%
    group_by(uniprot) %>% mutate(cv_strain = 100*sd_(int2use)/mean_(int2use))
  
  if(by_bio){
    CV = CV %>% 
      mutate(biological_rep=paste0(strain,'-',bio)) %>%
      group_by(uniprot,biological_rep) %>% mutate(cv_biorep = 100*sd_(int2use)/mean_(int2use)) %>%
      dplyr::select(uniprot,strain,cv_strain,biological_rep,cv_biorep) %>% distinct()
  }else{
    CV = CV %>% dplyr::select(uniprot,strain,cv_strain) %>% distinct()
  }
  
  return(CV)
}

draw_boxcv = function(CV, by_bio=T,plot=F){
  
  if( missing(CV) ){
    stop("run calculate_cv() on intensity dataframe in long format...")
  }
  
  cv1 = ggplot(CV) + 
    geom_boxplot(aes(y=cv_strain,x=strain,fill=strain),na.rm = T,show.legend = F) +  # trim=T,draw_quantiles = c(0.25,0.75), 
    theme(axis.text.x = element_text(angle=90,vjust = 0.6)) + 
    ylab('%CV within strains')
  if(by_bio){
    cv2 = ggplot(CV) + 
      geom_boxplot(aes(y=cv_biorep,x=biological_rep,fill=biological_rep), na.rm = T, show.legend = F) + #  trim=T,draw_quantiles = c(0.25,0.75), 
      theme(axis.text.x = element_text(angle=90,vjust = 0.6)) + 
      ylab('%CV within biological replicates')
    cv = (cv1|cv2)
  }else{
    cv = (cv1)
  }
  return(cv)
}

show_table_cv = function(CV,by_bio=T,caption=''){
  if( missing(CV) ){
    stop("run calculate_cv() on intensity dataframe in long format...")
  }

  if(by_bio){
    tab_cv = CV %>%
              group_by(biological_rep) %>%
              summarize(min=min_(cv_biorep),q25=q25(cv_biorep), md=median_(cv_biorep), q75=q75(cv_biorep), max=max_(cv_biorep)) %>% 
              kbl(x=., digits = 2,
                  caption = 'Coefficient of variation - intensities per biological replicate') %>%
              kable_paper("striped", full_width = F) %>% kable_styling()
  }else{ 
    tab_cv = CV %>% 
              group_by(strain) %>% 
              summarize(min=min_(cv_strain),q25=q25(cv_strain), md=median_(cv_strain), q75=q75(cv_strain), max=max_(cv_strain)) %>%
              kbl(.,digits = 2,
                  caption = 'Coefficient of variation - intensities per strain') %>% 
              kable_paper("striped", full_width = F) %>% kable_styling()
  }
  return(tab_cv)
}


## 4 Correlation of intensities ----------------------------------------------------------
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

# 5 Scatterplots of sample intensities -----------------------------------------
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
  return(p2)
}



norm_barplot = function(nr=NORM){
  
  library(NormalyzerDE)
  if(missing(nr)){
    stop("run first normalize_intensities(INT,DESIGN)")
  }
  nds <- nds(nr)

  methodlist <- normalizations(nr)
  filterED <- sampleReplicateGroups(nds)
  filterrawdata <- filterrawdata(nds)
  
}


draw_normalization_density = function(RAW = int_raw, INT=int_lfq, DESIGN = df.group){
  
  IDS = rownames(INT)
  NORM = normalize_intensities(INT,DESIGN)
  #boxplot(NORM@ner@avgcvmem)
  #boxplot(NORM@ner@avgmadmem)
  #boxplot(NORM@ner@avgvarmem)
  
  norm_methods = names(NORM@normalizations)
  
  #NORM@normalizations[[n]]$normalization = n
  
  df_raw = pivot_longer(data = RAW %>% rownames_to_column('uniprot'), cols = -uniprot, 
                        names_to = 'sample', values_to ='norm_int', values_transform= log2) %>% 
            mutate(type='raw intensity') %>%  left_join(DESIGN)

  norm_list = list()
  for(n in norm_methods ){
    norm_list[[n]] = pivot_longer(data = NORM@normalizations[[n]] %>% as.data.frame(row.names=IDS) %>% rownames_to_column('uniprot'),  
                                  cols = -uniprot,
                                  names_to = 'sample', values_to ='norm_int') %>% 
                     mutate(norm = n )
  }
  #norm_list$raw=df_raw
  df_norm = bind_rows(norm_list) %>% left_join(DESIGN) #%>% 
  #          left_join(df_raw,by=c('uniprot','sample'))
  
  library(geomtextpath)
  dens_plot = ggplot(df_norm,aes(x=norm_int)) + 
    geom_density(data=subset(df_raw, strain!='AMH'), mapping=aes(x=norm_int,group=strain),
                 col='black',show.legend = F,linewidth=1) +
    geom_textdensity(data=subset(df_raw, strain=='AMH'), mapping=aes(x=norm_int,group=strain,label=type),
                     col='black',show.legend = F, size=3,
                     fontface = 2,linewidth=1, hjust = 0.2, vjust = 1.5) +
    geom_density(aes(col=strain),show.legend = T) + 
    facet_wrap(~norm,nrow=2,ncol=4)  
  
  return(dens_plot)
    
  
  
  #plot(log2(INT[,1]),log2(INT[,2]))
  # #points(NORM@normalizations$Quantile[,1],NORM@normalizations$Quantile[,2],col='red')
  # 
  # library(limma)
  # 
  # n = ncol(rawexp)/2
  # s = colnames(rawexp)
  # 
  # par(pch=19,cex.lab=1.4,cex.axis=1.2)
  # cols = c( raw=rgb(0,0,0,0.5), loess=rgb(1,0,0,0.5), median=rgb(0,1,0,0.5), quantile=rgb(0,0,1,0.5))
  # for( i in 1:n){
  #   for( j in 1:n){
  #     exp_unit= "(log2)-intensity"
  #     plot( rawexp[,s[i]], rawexp[,s[j]], xlab='', ylab='',col=cols['raw'], cex=0.6)
  #     points(norm_loess[,s[i]], norm_loess[,s[j]],col=cols['loess'], cex=0.6)
  #     points(norm_median[,s[i]], norm_median[,s[j]],col=cols['median'], cex=0.6)
  #     points(norm_quantile[,s[i]], norm_quantile[,s[j]],col=cols['quantile'], cex=0.6)
  #     
  #     legend('topleft',title=NULL,legend = names(cols), col = cols, pch = 19, border = NA, bty='n')
  #     title(main="Raw vs Normalized protein expression:",
  #           xlab=paste(s[i],exp_unit) ,ylab=paste(s[j],exp_unit))
  #   }
  # }
  #  if(n>1 & onebyone){ invisible(readline(prompt="Press [enter] to continue")) }
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

# 6 Differential expression ----------------------------------------------------
library(limma)
pairwise_condition = function(conditions = all_strains,to_pair=T){
  cond_pair = expand.grid(S1=conditions,S2=conditions)
  true.pairs <- t(apply(cond_pair[,1:2], 1, sort))
  dup.pairs <- duplicated(true.pairs)
  unique.pairs <- true.pairs[!dup.pairs & true.pairs[,1]!=true.pairs[,2],] %>% 
    as_tibble() %>% rename(cond1=V1,cond2=V2)
  
  if(to_pair){
    return( apply(unique.pairs,1,paste0,collapse='-') )
  }
  return(unique.pairs)
}

compare_conditions = function(input, id_col = "uniprot",
                              comparison = pairwise_condition(all_strains,to_pair=T),
                              col_group=1){
  
  if( id_col %in% colnames(input)){
    processed_data = na.omit(input) %>% remove_rownames %>% column_to_rownames(var = id_col)
  }else{
    processed_data = na.omit(input)
  }
  
  samples=colnames(processed_data)
 
  GROUP = get_group(sample_names = colnames(processed_data),sep="_",grp_num = col_group)
  
  f.df <- factor(GROUP)
  design <- model.matrix(~0+f.df)
  colnames(design) <- levels(f.df)
  fit <- lmFit(processed_data, design)

  cont.matrix <- makeContrasts(contrasts = comparison, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit3 <- eBayes(fit2,trend=T,robust=T)
  return(fit3)
}

get_volcano_data = function(input_data=int_norm,
                            min_lfc=2, min_pval=0.01, which=c('both','up','down'), topn = 20){
  datList <- list()
  combination <- pairwise_condition(to_pair=T)
  if(missing(topn)){
    cat('selecting top 50 differentially expressed genes...\n')
    topn = 50
  }
    
  if(missing(which)){ 
    cat('selecting both up and down regulated genes...\n')
    which = 'both'
  }
  
  for (i in 1:length(combination)) {
    
    fit2 <- compare_conditions(input=input_data, id_col = 'uniprot', comparison=combination, col_group = 1)
    
    d.out <- data.frame(ID = names(fit2$coefficients[,i]),
                        pValue = fit2$p.value[,i],
                        qValue = p.adjust(fit2$p.value[,i], "fdr"),
                        EffectSize = fit2$coefficients[,i],
                        comparison = combination[i])
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
    datName <- combination[i]
    datList[[datName]] = d.out
  }

  return(datList)
}

get_dfe = function(INPUT=int_norm, MIN_LFC=2, MIN_PVAL=0.01, WHICH='both', TOPN = 20,
                   count_up=T, count_down=T){
  dfe <- get_volcano_data(input_data=INPUT, min_lfc=MIN_LFC, min_pval=MIN_PVAL, WHICH, topn = TOPN) %>% 
            bind_rows %>% 
            dplyr::filter(sig!="Non significant")
  
  if(count_down){
    down = dfe %>% group_by(ID) %>% dplyr::filter(sig=='Downregulated') %>%
      summarize( strains_down = paste0(comparison,collapse=' '),
                 down_AMH = str_count(strains_down,'AMH-'),
                 down_BAN = str_count(strains_down,'BAN-'),
                 down_BED = str_count(strains_down,'BED-'),
                 down_BPL = str_count(strains_down,'BPL-'),
                 down_BTT = str_count(strains_down,'BTT-'),
                 down_CMP = str_count(strains_down,'CMP-'),
                 down_CPI = str_count(strains_down,'CPI-'),
                 down_CQC = str_count(strains_down,'CQC-'))
    dfe = left_join(dfe,down)
    dfe[ grep(x=colnames(dfe),'^down_')] = 0
  }
  
  if(count_up){
    up = dfe %>% group_by(ID) %>% dplyr::filter(sig=='Upregulated') %>%  
      summarize( strains_up = paste0(comparison,collapse=' '),
                 up_AMH = str_count(strains_up,'AMH-'),
                 up_BAN = str_count(strains_up,'BAN-'),
                 up_BED = str_count(strains_up,'BED-'),
                 up_BPL = str_count(strains_up,'BPL-'),
                 up_BTT = str_count(strains_up,'BTT-'),
                 up_CMP = str_count(strains_up,'CMP-'),
                 up_CPI = str_count(strains_up,'CPI-'),
                 up_CQC = str_count(strains_up,'CQC-'))
    dfe = left_join(dfe,up) 
    dfe[ grep(x=colnames(dfe),'^up_')] = 0
  }
  return(dfe)
}

draw_volcano = function(data2plot,plot_title){
  sig_count = table(data2plot$sig)
  down = sig_count[1]
  up = last(sig_count)
  
  data_sig = subset(data2plot,sig!='Non significant')
  p <- ggplot(data2plot, 
              aes(x=EffectSize, y=-log10(qValue), fill = sig)) +
       geom_point(pch = 21, colour = "black", alpha = 0.5, size = 1.5)
  p <- p + scale_fill_manual(aesthetics = c('colour','fill'),
                             values=c("Non significant" = 'gray', 
                                      "Downregulated" = 'red', 
                                      "Upregulated" = 'blue',
                                      'is_imputed'='purple')) +
          xlab("log2 fold change") +
          ylab("-log10 q-value") + 
          labs(fill = NULL) +
          ggtitle(label = plot_title )
  
  p <- p +
       theme_classic(base_size = 14) +
       theme(legend.position = 'top',
             plot.title = element_text(face  = 1, hjust =0, size = 18),
             axis.title = element_text(size = 10)
             )

  p <- p + geom_text_repel(data=data_sig,aes(label = ID,col=sig), show.legend = FALSE)
  return(p)
}

volcPlot = function(INPUT=int_norm, IMPUTED, MIN_LFC=2, MIN_PVAL=0.01, WHICH='both', TOPN = 20, 
                    plot=F, use_plotly=T, use_label='genename'){
    
  plotList <- list()
  
  if(missing(IMPUTED)){ 
    IMPUTED = tibble( uniprot = INPUT$uniprot, 
                      is_imputed = (rowSums( is.na(INPUT))>0 )* 1, 
                      imputed = factor(is_imputed,levels = c(0,1), labels = c('not','is_imputed')))
  }
  #dlist <- get_volcano_data(input_data=INPUT, min_lfc=MIN_LFC, min_pval=MIN_PVAL, WHICH, topn = TOPN, id_col=use_label)
  all_data = get_volcano_data(INPUT, MIN_LFC, MIN_PVAL,  WHICH,TOPN) %>% bind_rows %>%
    mutate(log10_qvalue=-log10(qValue)) %>%
    dplyr::left_join(sc_identifiers, by=c('ID'='UNIPROT'), keep=T ) %>% 
    dplyr::left_join(IMPUTED,by=c('ID'='uniprot'))
  
  id2show = names(sc_identifiers)
  label2use=match.arg(tolower(use_label),tolower(id2show)) 
  col_label = grep(label2use,id2show,ignore.case = T,v=T)
  all_data$ID = all_data[[col_label]]
  
  comps = unique(all_data$comparison)
  
  for( comp in comps ){
    data_comp = subset(all_data, comparison == comp )
    p = draw_volcano(data_comp,comp)
    p = p + geom_point(data=subset(data_comp,is_imputed==1), aes(col=imputed), fill=NA, stroke=1,shape=21,show.legend = F)
    p = p + geom_vline(aes(xintercept = (MIN_LFC*-1)), lty = 'dashed', colour = 'red',lwd=0.5) +
            geom_vline(aes(xintercept = (MIN_LFC)), lty = 'dashed', colour = 'blue',lwd=0.5) +
            geom_hline(aes(yintercept = -log10(MIN_PVAL)), lty = 'dashed', colour = 'black',lwd=0.5) 
    
    if(plot && !use_plotly)
      plot(p)
    plotList[[comp]] = p
  }
  
  if(use_plotly){
    library(plotly)
    
    button_comparisons = lapply(comps, FUN = function(comp) {
     button <- list(method = 'restyle', args = list('transforms[0].value', comp), label = comp)
    })
    button_imputed <- list(
      list( method = 'restyle', args = list('transforms[1].value', 1), label = 'with_impute'),
      list( method = 'restyle', args = list('transforms[1].value', 0), label = 'no_impute'))
    button_ids = lapply(id2show, FUN = function(id) {
      button <- list(method = 'restyle', args = list("text",list(as.formula(paste0("~",id)))),label = id  )
    })
    
    
    count_data = all_data %>%
                 group_by(comparison) %>% 
                 mutate(Y = max(log10_qvalue), Y_imp =Y - 1*is_imputed) %>%
                 ungroup() %>% mutate(X = fct_recode(factor(sig),"-2"="Downregulated","0"="Non significant","2"="Upregulated") %>% unfactor) %>%
                 group_by(comparison,sig) %>% mutate(n=n()) %>% 
                 add_count(is_imputed,name = 'n_imp')
    
    #colnames(df) <- c("x", "y")
    fig <- plot_ly(data=all_data %>% arrange(sig) %>% dplyr::left_join(count_data),
                   type='scatter', mode='markers',
                   x=~EffectSize, y=~log10_qvalue, color=~sig, 
                   text = ~ID, hoverinfo = 'text', alpha = 0.3, sizes = 0.8,
                   colors=c("Non significant" = 'gray', "Downregulated" = 'red', "Upregulated" = 'blue'),
                   transforms = list(
                        list(type = 'filter', target = ~comparison, operation = '=', value = sort(comps)[1]),
                        list(type = 'filter', target = ~is_imputed, operation = '<=', value = 1)
                   )
                 ) %>%
      add_text(x=~X, y=~Y_imp, customdata =  ~n_imp,  showlegend=F, texttemplate = '%{customdata:.s}', textposition = 'outside', textfont = list(size=20), hovertemplate = ~sig) %>%
      plotly::layout( updatemenus = list(
                          list( y=1,x=-0.1,type='dropdown', active = 0, buttons = button_comparisons, name='comparison'),
                          list( y=0.2,x=-0.1,type='buttons', buttons = button_imputed, name='with_impute' ),
                          list( y=0.7,x=-0.1,type='buttons', active = 2, buttons = button_ids,name='id text')),#legend = list(x = 0, y = 100,orientation='h'),
                       uniformtext=list(minsize=16, mode='hide'))#,
                      #annotations = list(list(text = "Compare:", font=list(size=16), x=-0.2, y=1.05, xref='paper', yref='paper',showarrow=F),
                       #                  list(text="Show:", font=list(size=16), x=-0.2, y=0.9, xref='paper', yref='paper',showarrow=F))) %>%
      
    
    return(fig)
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

draw_umap_DE = function(EXP,DE){
  
  library(umap)
  set.seed(142)
  umap_fit <- EXP %>%
    select(where(is.numeric)) %>%
    column_to_rownames("ID") %>%
    scale() %>% 
    umap()
  
  umap_df <- umap_fit$layout %>%
    as.data.frame()%>%
    rename(UMAP1="V1", UMAP2="V2") %>%
    mutate(ID=row_number()) #%>%
    #inner_join(penguins_meta, by="ID")
  
  U = umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2,
               color = species)) +
    geom_point(size=3, alpha=0.5)+
    #facet_wrap(~island)+
    labs(x = "UMAP1", y = "UMAP2", subtitle="UMAP plot")+
    theme(legend.position="bottom")
    #ggsave("UMAP_plot_example2.png")
  #plot(U)
  return(U)
}
