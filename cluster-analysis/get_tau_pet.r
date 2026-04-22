get_tau_pet <- function() {
  # TAU PET
  ## Partial Volume Correction (PVC)
  # Tau data corrected for partial volume effects using the Geometric Transfer Matrix (GTM) approach. 
  # Using the MRI closest in time to tau scan, the GTM approach models all FreeSurfer-defined ROIs as well as regions in which offtarget
  # binding is common (e.g., choroid plexus in FTP; meninges in MK6240) to reduce
  # contamination from these regions into neighboring regions of interest. 
  require(ADNIMERGE2)
  library(epitools)
  library(lubridate)
  
  tau.pet.pvc <- ADNIMERGE2::UCBERKELEY_TAUPVC_6MM %>% 
    filter(TRACER == "FTP")
  
  bl.df <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% distinct(RID, .keep_all = TRUE)
  
  tau.pet.pvc <- left_join(tau.pet.pvc, bl.df, by = "RID") %>% 
    mutate(Months = interval(EXAMDATE.bl, SCANDATE) %/% months(1),
           Years = Months/12)
  
  return(tau.pet.pvc)
}