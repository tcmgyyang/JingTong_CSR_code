suppressMessages({library(magrittr); library(ieugwasr)})
tf <- file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE)[1])),"opengwas_token.txt")
if (nchar(Sys.getenv("OPENGWAS_JWT"))<10 && file.exists(tf)) Sys.setenv(OPENGWAS_JWT=trimws(readLines(tf,warn=FALSE)[1]))
cat("token length:", nchar(Sys.getenv("OPENGWAS_JWT")), "\n")
u <- tryCatch(user(), error=function(e) NULL)
cat("auth:", if(!is.null(u)) paste("OK ->", u$user$uid) else "FAILED", "\n\n")

outcomes <- c("finn-b-M13_CERVICDISC","finn-b-M13_SPONDYLOSIS","finn-b-M13_INTERVERTEB")
expo     <- c(IL6="eqtl-a-ENSG00000136244", IL1B="eqtl-a-ENSG00000125538",
              MMP9="eqtl-a-ENSG00000100985", TP53="eqtl-a-ENSG00000141510",
              AKT1="eqtl-a-ENSG00000142208", JUN="eqtl-a-ENSG00000177606",
              MMP3="eqtl-a-ENSG00000149968", VCAM1="eqtl-a-ENSG00000162692",
              EGFR="eqtl-a-ENSG00000146648", IL1A="eqtl-a-ENSG00000115008",
              IL4="eqtl-a-ENSG00000113520")
cat("== outcomes ==\n")
for (o in outcomes){
  info <- tryCatch(gwasinfo(o), error=function(e) NULL)
  cat(sprintf("  %-26s %s\n", o, if(!is.null(info)&&nrow(info)>0) paste0("OK ncase=",info$ncase," ncontrol=",info$ncontrol) else "NOT FOUND"))
}
cat("== exposures (eQTLGen) ==\n")
for (i in seq_along(expo)){
  info <- tryCatch(gwasinfo(expo[[i]]), error=function(e) NULL)
  cat(sprintf("  %-6s %-26s %s\n", names(expo)[i], expo[[i]], if(!is.null(info)&&nrow(info)>0) "OK" else "NOT FOUND"))
}
