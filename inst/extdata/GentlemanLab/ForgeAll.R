#######################################################################
### Use this R script to forge all the BSgenome data packages
###

library(BSgenome)

pkgs <- c(
  "BSgenome.Amellifera.BeeBase.assembly4",
  "BSgenome.Amellifera.UCSC.apiMel2",
  "BSgenome.Athaliana.TAIR.01222004",
  "BSgenome.Athaliana.TAIR.04232008",
  "BSgenome.Btaurus.UCSC.bosTau3",
  "BSgenome.Btaurus.UCSC.bosTau4",
  "BSgenome.Celegans.UCSC.ce2",
  "BSgenome.Cfamiliaris.UCSC.canFam2",
  "BSgenome.Drerio.UCSC.danRer5",
  "BSgenome.Dmelanogaster.UCSC.dm2",
  "BSgenome.Dmelanogaster.UCSC.dm3",
  "BSgenome.Ecoli.NCBI.20080805",
  "BSgenome.Ggallus.UCSC.galGal3",
  "BSgenome.Hsapiens.UCSC.hg17",
  "BSgenome.Hsapiens.UCSC.hg18",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Mmusculus.UCSC.mm8",
  "BSgenome.Mmusculus.UCSC.mm9",
  "BSgenome.Ptroglodytes.UCSC.panTro2",
  "BSgenome.Rnorvegicus.UCSC.rn4",
  "BSgenome.Scerevisiae.UCSC.sacCer1"
)

for (pkg in pkgs) {
  cat("\n")
  cat("============================================================\n")
  cat("START FORGING ", pkg, "\n", sep="")
  cat("\n")
  forgeBSgenomeDataPkg(pkg)
  cat("\n")
  cat("END FORGING ", pkg, "\n", sep="")
}

