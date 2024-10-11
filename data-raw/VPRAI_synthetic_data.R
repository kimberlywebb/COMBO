## code to prepare `LSAC_data` dataset goes here

save_wd <- "C:/Users/KIW58/OneDrive - University of Pittsburgh/Documents/Misclassification Research/R Packages/GitHub/COMBO/data-raw/"
VPRAI_synthetic_data <- read.csv(paste0(save_wd, "VPRAI_Synthetic_Data.csv"))

usethis::use_data(VPRAI_synthetic_data, overwrite = TRUE)
