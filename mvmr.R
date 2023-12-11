rm(list=ls())
library(data.table)
library(TwoSampleMR) 

#读取暴露数据
ADHD<-fread('ADHD2022_iPSYCH_deCODE_PGC.meta.gz')
head(ADHD)
ADHD$BETA<-log(ADHD$OR)
exposure1<-format_data(ADHD,type = 'exposure',snp_col = "SNP", beta_col = "BETA", se_col = "SE", 
                       effect_allele_col = "A1", other_allele_col = "A2", 
                       pval_col = "P", chr_col = "CHR",  
                       eaf_col = "FRQ_U_186843",
                       pos_col = "BP")

finngen<-fread('finngen_R9_F5_DEPRESSIO.gz')
head(finngen)
exposure2<-format_data(finngen,type = 'exposure',snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
                       eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", 
                       pval_col = "pval", chr_col = "#chrom", 
                       pos_col = "pos")

#合并
library(dplyr)
#exposure2<-select(exposure2,-c(gene.outcome))
new_order <- c('chr.exposure','SNP','pos.exposure','effect_allele.exposure',
               'other_allele.exposure','se.exposure','pval.exposure',
               'beta.exposure','exposure','mr_keep.exposure','pval_origin.exposure','id.exposure','eaf.exposure')  
a <- select(exposure2, new_order)
b<-exposure1
a$Phenotype <- "depression"
b$Phenotype <- "ADHD"
exp<-rbind(a,b,fill=T)

temp<-exp
temp$id.exposure <- 1
temp <- temp[order(temp$pval.exposure, decreasing = FALSE),]
temp <- subset(temp, !duplicated(SNP))

#本地clump
library(ieugwasr)
library(plinkbinr)
get_plink_exe()
head(temp)
clumped <- ld_clump(
  dplyr::tibble(rsid=temp$SNP, pval=temp$pval.exposure, id=temp$id.exposure),
  plink_bin = "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/plinkbinr/bin/plink_Darwin",
  bfile = "./EUR_1kgref/EUR")

final_exposure <- subset(exp, SNP %in% clumped$rsid)
head(final_exposure)

#取并集
m1<-subset(a,SNP %in% final_exposure$SNP)
m2<-subset(b,SNP %in% final_exposure$SNP)
m11 <- subset(m1,m1$SNP %in% m2$SNP)
m22 <- subset(m2,m2$SNP %in% m1$SNP)
exposure_dat <- rbind(m11,m22)

exposure_dat$pval <- as.numeric(exposure_dat$pval.exposure)
head(exposure_dat)
exposure_dat<- format_data(exposure_dat,type = 'exposure',snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", 
                           effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", 
                           pval_col = "pval.exposure", chr_col = "chr.exposure",  
                           eaf_col = "eaf.exposure",
                           pos_col = "pos.exposure")

#读取结局数据
finngen2<-fread('finngen_R9_F5_PTSD.gz')
head(finngen2)
outcome<-format_data(finngen2,type = 'outcome',snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
                       eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref", 
                       pval_col = "pval", chr_col = "#chrom", 
                       pos_col = "pos")
outcome_dat <- subset(outcome,SNP %in% exposure_dat$SNP)

#MVMR
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat) 
res <- mv_multiple(mvdat)
res

res_OR<-generate_odds_ratios(res$result)
res_OR

write.csv(res_OR, 'OR_mainresults.csv')
