#######################################################################
### Use this R script to forge all the BSgenome data packages
###

library(BSgenome)

pkgs <- c(
  "BSgenome.Amellifera.BeeBase.assembly4",
  "BSgenome.Amellifera.UCSC.apiMel2",
  "BSgenome.Athaliana.TAIR.04232008",
  "BSgenome.Athaliana.TAIR.TAIR9",
  "BSgenome.Btaurus.UCSC.bosTau3",
  "BSgenome.Btaurus.UCSC.bosTau4",
  "BSgenome.Celegans.UCSC.ce2",
  "BSgenome.Celegans.UCSC.ce6",
  "BSgenome.Celegans.UCSC.ce10",
  "BSgenome.Cfamiliaris.UCSC.canFam2",
  "BSgenome.Cfamiliaris.UCSC.canFam3",
  "BSgenome.Dmelanogaster.UCSC.dm2",
  "BSgenome.Dmelanogaster.UCSC.dm3",
  "BSgenome.Drerio.UCSC.danRer5",
  "BSgenome.Drerio.UCSC.danRer6",
  "BSgenome.Drerio.UCSC.danRer7",
  "BSgenome.Ecoli.NCBI.20080805",
  "BSgenome.Gaculeatus.UCSC.gasAcu1",
  "BSgenome.Ggallus.UCSC.galGal3",
  "BSgenome.Ggallus.UCSC.galGal4",
  "BSgenome.Hsapiens.UCSC.hg17",
  "BSgenome.Hsapiens.UCSC.hg18",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Mmulatta.UCSC.rheMac2",
  "BSgenome.Mmusculus.UCSC.mm8",
  "BSgenome.Mmusculus.UCSC.mm9",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Ptroglodytes.UCSC.panTro2",
  "BSgenome.Ptroglodytes.UCSC.panTro3",
  "BSgenome.Rnorvegicus.UCSC.rn4",
  "BSgenome.Rnorvegicus.UCSC.rn5",
  "BSgenome.Scerevisiae.UCSC.sacCer1",
  "BSgenome.Scerevisiae.UCSC.sacCer2",
  "BSgenome.Scerevisiae.UCSC.sacCer3"
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

