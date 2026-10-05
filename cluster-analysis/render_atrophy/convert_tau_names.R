
convert_oasis_tau <- function(x) {
  name_sheet <- read.csv("~/R/EDAP-data/OASIS/OASIS3_data_files/scans/dictionaries-Imaging_and_UDS_data_dictionaries/resources/pdf/files/oasis_pup_variable_crosswalk.csv") %>%
    mutate(NAME = toupper(Structure_Name)) %>%
    rename(SAS_label = SAS_Compatible_Variable_Label)
  
  # Remove PET_fSUVR_rsf_ prefix
  z <- str_remove(x, "^PET_fSUVR_rsf_")
  
  name <- filter(name_sheet, SAS_label == z)$NAME
  
  if (length(name) == 0) {
    return(x)
  } else {
    return(paste0(name, "_SUVR"))
  }
}



