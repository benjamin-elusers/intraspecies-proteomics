reprex::reprex({
  library(datapasta)

  find_duplicate_pairs = function(vec1, vec2){
    # Find duplicates of paired values between two vectors
      purrr::map2(vec1, vec2, ~ sort(c(.x, .y))) |> duplicated()
  }

  #### 1. get expression matrix of all samples ####
  prot_exp = df_int_bpca |> relocate(uniprot)
  dim(prot_exp)
  colnames(prot_exp)

  #### 2. define groups based on sample names ####
  df_strains = tibble(
                  # strain names
                  strain = c('AMH','BAN','BED','BPL','BTT','CMP','CPI','CQC'),
                  # strain origins
                  origin = c('wild','wild','wild', 'domesticated', 'domesticated', 'wild',  'domesticated',  'domesticated')
  )

  #  sample names are constructed on four identifiers separated by an underscore ("_"):
  #  1. strain name
  #  2. biological replicate number
  #  3. technical replicate number
  #  4. day of run
  df.group = tibble(sample = colnames(prot_exp)[-1]) |>
             separate(col  = sample, sep = '_', remove = F,
                      into = c('strain','biorep','techrep','dayrun')) |>
             left_join(df_strains,by='strain') |>
             mutate( strain  = factor(strain,levels = strains),
                     biorep  = paste0("R",biorep),
                     techrep = paste0("r",techrep),
                     dayrun  = paste0(dayrun,"_04_22")) # Experiment from April 2022

  #### 3. design matrix for experiment ####
  design <- model.matrix(~0+df.group$strain)
  colnames(design) <- strains

  #### 4. fit gene linear model given design matrix ####
  fit <- lmFit(prot_exp, design)

  #### 5. contrasts for differential expression ####
  CUTOFF_LFC <- 2 # (log2(lfc)>2 => above 4-fold)
  CUTOFF_PV  <- 1e-5 # False discovery rate < 1/100000 (less than one false-positive for all comparisons)

  ####__5.a) "all-vs-all" pairwise strategy ####
  all_pairs <- pair_strains %>% filter(!is_identical) %>% pull(var=pair,name=S1)
  strain_by_pair <- split(all_pairs,f = names(all_pairs))

  pairwise_contrast <- makeContrasts(contrasts = all_pairs, levels = design)
  pairwise_fit <- contrasts.fit(fit, pairwise_contrast)
  pairwise_DE <- eBayes(pairwise_fit,trend=T,robust=T)

  ####__5.b) "one-vs-all" strategy ####
  one_contrast <- makeContrasts(contrasts = strains,levels = design)
  one_contrast[one_contrast<1] <- -1 / (nlevels(df.group$strain)-1) # weighted mean expression from all strains
  one_fit <- contrasts.fit(fit, one_contrast)
  one_DE <- eBayes(one_fit,trend=T,robust=T)
  topTable(one_DE,coef = 1,lfc=CUTOFF_LFC,p.value=CUTOFF_PV,adjust.method='fdr',number = Inf)

  ####__5.c) check log-foldchange for the first gene ####
  gene_1 = prot_exp[1,-1] |> unlist()

  ####______ using strain mean expression ####
  strain_avg_1 = tapply(gene_1, df.group$strain, mean)
  all.equal(mean(strain_avg_1),one_fit$Amean[1]) # Average expression of all strains

  lfc_strain_1 = strain_avg_1[1] - mean(strain_avg_1[-1]) # (mean AMH) vs (mean not-AMH)
  all.equal(lfc_strain_1,one_fit$coefficients[1]) # log-foldchange with strain mean expression

  ####______ using sample average expression (AMH vs. all other strains) ####
  AMH_1 = gene_1[1:4] |> mean()
  notAMH_1 = gene_1[5:32] |> mean()
  lfc_sample_1 = AMH_1 - notAMH_1
  all.equal(lfc_sample_1,one_fit$coefficients[1])

  session_info()
})

#### 5. differentially expressed genes ####


pairwise_DE_prot = map(strain_by_pair,
                       ~topTable(pairwise_DE,  coef = .x , lfc = CUTOFF_LFC, p.value = CUTOFF_PV, adjust.method = 'fdr', number = Inf) %>% pull(uniprot)
)


one_DE_prot = map(strains,
                  ~topTable(one_DE,  coef = .x , lfc = CUTOFF_LFC, p.value = CUTOFF_PV, adjust.method = 'fdr', number = Inf) %>% pull(uniprot)
)

#setdiff(sort(pairwise_DE_prot$AMH) , sort(pairwise_DE_prot[-1] %>% purrr::flatten_chr() %>% janitor::tabyl() %>% filter(n>4) %>% pull(uniprot) ))
#setdiff(sort(DE_strains$AMH) , sort(DE_strains[-1] %>% purrr::flatten_chr() %>% janitor::tabyl() %>% filter(n>4) %>% pull(1) ))
