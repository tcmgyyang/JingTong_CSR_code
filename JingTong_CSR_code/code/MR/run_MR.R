# =============================================================================
# Two-sample Mendelian randomization: hub-gene expression (eQTLGen) -> cervical
# spondylosis / disc disorders (FinnGen). Tests whether Jingtong-Granules hub genes
# are CAUSALLY associated with the disease (robust to the small transcriptome n).
# Requires: R + TwoSampleMR + ieugwasr, and an OpenGWAS token.
#   install.packages("remotes"); remotes::install_github("MRCIEU/TwoSampleMR")
#   Sys.setenv(OPENGWAS_JWT="<your token from https://api.opengwas.io>")
# Run: Rscript run_MR.R
# =============================================================================
suppressMessages({library(magrittr); library(TwoSampleMR); library(ggplot2)})
tf <- file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1])),"opengwas_token.txt")
if (nchar(Sys.getenv("OPENGWAS_JWT"))<10 && file.exists(tf)) Sys.setenv(OPENGWAS_JWT=trimws(readLines(tf,warn=FALSE)[1]))
setwd(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1])))

# --- hub genes -> eQTLGen exposure IDs in OpenGWAS (cis-eQTL of blood expression) ---
# 18 multi-source hub genes (>=2-source) UNION the original 11 hubs = 22 genes; the script
# skips any without eQTLGen instruments (handled by tryCatch/next below).
exposures <- c(
  AKT1 ="eqtl-a-ENSG00000142208", IL6  ="eqtl-a-ENSG00000136244", TNF  ="eqtl-a-ENSG00000232810",
  IL1B ="eqtl-a-ENSG00000125538", ESR1 ="eqtl-a-ENSG00000091831", MYC  ="eqtl-a-ENSG00000136997",
  EGFR ="eqtl-a-ENSG00000146648", JUN  ="eqtl-a-ENSG00000177606", HIF1A="eqtl-a-ENSG00000100644",
  MMP9 ="eqtl-a-ENSG00000100985", CCL2 ="eqtl-a-ENSG00000108691", CXCL8="eqtl-a-ENSG00000169429",
  BDNF ="eqtl-a-ENSG00000176697", IL10 ="eqtl-a-ENSG00000136634", IL2  ="eqtl-a-ENSG00000109471",
  IL4  ="eqtl-a-ENSG00000113520", BCL2L1="eqtl-a-ENSG00000171552", NFKBIA="eqtl-a-ENSG00000100906",
  TP53 ="eqtl-a-ENSG00000141510", MMP3 ="eqtl-a-ENSG00000149968", VCAM1="eqtl-a-ENSG00000162692",
  IL1A ="eqtl-a-ENSG00000115008")

# --- FinnGen outcomes (cervical spondylosis / disc disorders) ---
outcomes <- c(cervical_disc="finn-b-M13_CERVICDISC",
              spondylosis  ="finn-b-M13_SPONDYLOSIS",
              interverteb  ="finn-b-M13_INTERVERTEB")

dir.create("MR_results", showWarnings=FALSE)
all_res <- list()
for (oc_name in names(outcomes)) {
  oc <- outcomes[[oc_name]]
  for (g in names(exposures)) {
    message(sprintf(">>> %s -> %s", g, oc_name))
    res <- tryCatch({
      # 1. instruments: genome-wide significant, LD-clumped cis-eQTLs
      inst <- extract_instruments(exposures[[g]])
      if (is.null(inst) || nrow(inst) < 1) next
      # 2. outcome SNP effects
      out  <- extract_outcome_data(snps = inst$SNP, outcomes = oc)
      # 3. harmonise alleles
      dat  <- harmonise_data(inst, out)
      if (nrow(dat) < 1) next
      # 4. MR (IVW for >=2 SNPs; Wald ratio for single cis-SNP) + sensitivity
      mr_res <- mr(dat)
      mr_res$gene <- g; mr_res$outcome_name <- oc_name
      het  <- tryCatch(mr_heterogeneity(dat), error=function(e) NULL)
      pleio<- tryCatch(mr_pleiotropy_test(dat), error=function(e) NULL)
      list(mr=generate_odds_ratios(mr_res), het=het, pleio=pleio, dat=dat)
    }, error=function(e){ message("  failed: ", conditionMessage(e)); NULL })
    if (!is.null(res)) all_res[[paste(g,oc_name,sep="_")]] <- res
  }
}

# --- collate + forest plot ---
mr_tab <- do.call(rbind, lapply(all_res, function(x) x$mr))
write.csv(mr_tab, "MR_results/MR_estimates.csv", row.names=FALSE)
ivw <- subset(mr_tab, method %in% c("Inverse variance weighted","Wald ratio"))
if (nrow(ivw) > 0) {
  ivw$label <- paste(ivw$gene, ivw$outcome_name, sep=" -> ")
  p <- ggplot(ivw, aes(x=or, y=label)) +
    geom_point() + geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.2) +
    geom_vline(xintercept=1, linetype=2) + theme_bw() +
    labs(x="OR (95% CI) per SD increase in gene expression", y="",
         title="Causal effect of hub-gene expression on cervical spondylosis (MR)")
  ggsave("MR_results/Fig_MR_forest.pdf", p, width=8, height=6)
  ggsave("MR_results/Fig_MR_forest.png", p, width=8, height=6, dpi=300)
}
message("DONE -> MR_results/MR_estimates.csv + Fig_MR_forest.*")
