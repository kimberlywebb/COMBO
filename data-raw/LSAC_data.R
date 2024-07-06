## code to prepare `LSAC_data` dataset goes here

save_wd <- "C:/Users/kimho/Dropbox/Misclassification/Code/RPackages/GitHub/COMBO/data-raw/"
LSAC <- read.csv(paste0(save_wd, "LSAC_National_Bar_Passage_Study.csv"))

LSAC_data <- LSAC

usethis::use_data(LSAC_data, overwrite = TRUE)
