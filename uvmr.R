library(TwoSampleMR)
library('data.table')

#读取暴露数据
exposure0 <- fread("finngen_R9_U22_COVID19_SUSPECTED.gz",header = T)
 # #筛选强相关的变量：若5E-8筛选出来的变量较少，可适当调大P值（须有文献根据）
 exposure1<-subset(exposure0,pval<5e-8)
 exposure<-format_data(exposure1,type="exposure",
                       snp_col = "rsids",
                       chr_col ="#chrom",
                       phenotype_col = "COVID19",
                       beta_col = "beta",
                       se_col = "sebeta",
                       eaf_col="af_alt",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       pval_col = "pval")

 # #去除连锁不平衡（linkage disequilibrium）
 clumped <- ld_clump(
  dplyr::tibble(rsid=exposure$SNP, pval=exposure$pval.exposure, id=exposure$id.exposure),
  #get_plink_exe()
  plink_bin = "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/plinkbinr/bin/plink_Darwin",
  #欧洲人群参考基因组位置
  bfile = "./EUR_1kgref/EUR")
 exposure_data <- subset(exposure, SNP %in% clumped$rsid)

 
#读取结局数据
outcome <- fread("ADHD2022_iPSYCH_deCODE_PGC.meta.gz",header = T)
outcome_data<-format_data(outcome,
                          type = "outcome",
                          snp_col = "SNP",
                          beta_col = "OR",
                          se_col = "SE",
                          eaf_col = "FRQ_U_186843",
                          effect_allele_col = "A1",
                          other_allele_col = "A2",
                          pval_col = "P",
                          chr_col = "CHR")
                     
  dat <- harmonise_data(
    exposure_dat = exposure_data, 
    outcome_dat = outcome_data
  )
  #perform MR
  res <- mr(dat)
  res
  
  md <- mr_method_list()
  res <-mr(dat, method_list = c("mr_ivw_mre","mr_egger_regression","mr_weighted_median"))
  
  or_res <- generate_odds_ratios(res)
  
  #Heterogeneity P>0.05是有意义还是<0.05有意义？
  heter<- mr_heterogeneity(dat)
  
  #Horizontal pleiotropy P>0.05是有意义还是<0.05有意义？
  pleio<- mr_pleiotropy_test(dat) 
  
  
  # Append or_res to the summary_or_res dataframe
  summary_or_res <- rbind(summary_or_res, or_res)
  summary_heter <- rbind(summary_heter, heter)
  summary_pleio <- rbind(summary_pleio, pleio)
  
  
  write.csv(summary_or_res, 'OR_mainresults.csv')
  write.csv(summary_heter, 'Heterogeneity.csv')
  write.csv(summary_pleio, 'Horizontal_pleiotropy.csv')