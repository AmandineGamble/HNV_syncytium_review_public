f_data_format <- function(data) {
  
  data = subset(data, (is.na(treatment) | treatment == "control" | treatment == "pre-treatment") & 
                  # !is.element(virus, c("HeV_NiV_chimera", "NiV_mutant")) &
                  !is.na(virus_strain))
  
  data = subset(data, select = -c(treatment_start, treatment_start_unit, 
                                  time_unit))
  
  
  data[grepl("ig", data$PCR_value) | grepl("ab", data$PCR_value), "PCR_value"] = NA
  data[grepl("ig", data$PCR_pos_prop) | grepl("ab", data$PCR_pos_prop), "PCR_pos_prop"] = NA
  data[grepl("ig", data$viral_titer_value) | grepl("ab", data$viral_titer_value), "viral_titer_value"] = NA
  data[grepl("ig", data$viral_titer_pos_prop) | grepl("ab", data$viral_titer_pos_prop), "viral_titer_pos_prop"] = NA
  data[grepl("ig", data$neutralizing_ab_value) | grepl("ab", data$neutralizing_ab_value), "neutralizing_ab_value"] = NA
  data[grepl("ig", data$neutralizing_ab_pos_prop) | grepl("ab", data$neutralizing_ab_pos_prop), "neutralizing_ab_pos_prop"] = NA
  
  # Id
  # Take individual id if available, group id otherwise
  data$id = as.character(data$individual_id)
  data[is.na(data$id), "id"] = data[is.na(data$id), "study_group"]
  data$id = paste0(data$label, "_", data$id)
  
  # Order by id
  data = data[order(data$id), ]
  
  # Save the info whether it is a precise group size or not 
  # Here anything with >2 characters is a range or group of values
  data$group_size = as.character(data$group_size)
  data$group_size_type = NA
  data[which(nchar(data$group_size) > 2), "group_size_type"] = "minimum"
  data[which(nchar(data$group_size) <= 2), "group_size_type"] = "exact"
  # Need to CHECK: whether there is a count per sample size for ranges and group of values
  # to calculate a proper arithmetic mean
  
  # QUESTION: HOW TO DEAL WITH INEXACT SAMPLE SIZES?
  
  # CHECK: take the smallest option (to be conservative on sample size) but there might be better options
  data[which(data$group_size == "1-2"), "group_size"] = 1
  data[which(data$group_size == "7-8"), "group_size"] = 7
  data[which(data$group_size == "13-19"), "group_size"] = 13
  data[which(data$group_size == "49-55"), "group_size"] = 49
  
  # Another option could be to take the mean, depends on the question and the info available in the publication
  #data[which(data$group_size == "1-2"), "group_size"] = mean(c(1, 2))
  #data[which(data$group_size == "7-8"), "group_size"] = mean(c(7, 8))
  #data[which(data$group_size == "13-19"), "group_size"] = mean(c(13, 19))
  #data[which(data$group_size == "49-55"), "group_size"] = mean(c(49, 55))
  
  # Transform as a numeric value
  data$group_size = as.numeric(data$group_size)
  
  # List of individuals or groups positive at least once (i.e., succesfully infected)
  
  # Information available
  data$infection_screened = 0
  data[!is.na(data$PCR_pos_prop) | !is.na(data$viable_virus_detect_prop) |
         !is.na(data$viral_titer_pos_prop) | !is.na(data$IHC_pos_prop) |
         !is.na(data$PCR_value) | !is.na(data$viable_virus_detect) |
         !is.na(data$viral_titer_value) | !is.na(data$IHC_pos)
       , "infection_screened"] = 1
  
  # One positive test
  data$infection_pos = NA
  data[data$infection_screened == 1, "infection_pos"] = 0
  data[which(as.numeric(data$PCR_pos_prop) > 0 | as.numeric(data$viable_virus_detect_prop) > 0 |
               as.numeric(data$viral_titer_pos_prop) > 0 | as.numeric(data$IHC_pos_prop) > 0 |
               as.numeric(data$PCR_value) > 0 | as.numeric(data$viable_virus_detect) > 0 |
               as.numeric(data$viral_titer_value) > 0 | as.numeric(data$IHC_pos) > 0)
       , "infection_pos"] = 1

  #data$tissue_group = as.character(data$tissue_group)
  #data$tissue_group = as.factor(data$tissue_group)
  data$tissue_group = factor(data$tissue_group,
                             levels = c("respiratory tract", "nasal","digestive system",
                                        "lymphatic system", "nervous system", 
                                        "cardiovascular", "urinary tract",
                                        "reproductive system",  
                                        "other", "organ_mix", "organ/fluid/swab_mix", "fluid or swab"))
  # What are NAs?
  data$tissue = factor(as.character(data$tissue))
  data$tissue = factor(data$tissue, levels = unique(data$tissue[order(data$tissue_group,as.character(data$tissue))]))
    

  
  data$species = factor(data$species, 
                        levels =  c("human", "african_green_monkey", "cynomolgus_monkey", "squirrel_monkey", "pig", "dog", "cat", "ferret", "guinea_pig", "rabbit",
                                    "hamster", "rat", "house_mouse", "mouse_NSG_xeno_human_lung", "mouse_IFNAR-KO", "horse", "black_flying_fox", "grey_headed_flying_fox", "large_flying_fox", "egyptian_fruit_bat", "chicken"))

  data$virus = factor(data$virus, 
                        levels =  c("HeV", "NiV_M", "NiV_B", "NiV_M/HeV", "CeV", "NiV_M_mutant"))
  
  
  #####
  # Keep only individuals with syncytium data
  
  # Generate column informing the number of individuals for which syncytium data are reported
  data[!is.na(data$syncytium_prop), "syncytium_screen"] = 1 
  data[is.na(data$syncytium_prop), "syncytium_screen"] = 0 
  
  data$syncytium_screen = as.numeric(as.character(data$syncytium_screen))
  data$group_size = as.numeric(as.character(data$group_size))
  data$syncytium_prop_raw = data$syncytium_prop
  data$syncytium_prop = as.numeric(as.character(data$syncytium_prop))
  
  data$syncytium_screen = data$syncytium_screen * data$group_size
  
  # Number of individuals with syncytia
  data$syncytium_n = round(data$syncytium_prop * data$syncytium_screen, 0)
  
  # Proportion of individuals with syncytia among the screened ones
  data$syncytium_p = data$syncytium_n / data$syncytium_screen
  
  
  
  #####
  # Keep only individuals with IHC data
  
  # Generate column informing the number of individuals for which IHC data are reported
  data[!is.na(data$IHC_pos_prop), "IHC_screen"] = 1 
  data[is.na(data$IHC_pos_prop), "IHC_screen"] = 0 
  
  data$IHC_screen = as.numeric(as.character(data$IHC_screen))
  data$group_size = as.numeric(as.character(data$group_size))
  data$IHC_pos_prop_raw = data$IHC_pos_prop
  data$IHC_pos_prop = as.numeric(as.character(data$IHC_pos_prop))
  
  data$IHC_screen = data$IHC_screen * data$group_size
  
  # Number of individuals with syncytia
  data$IHC_n = round(data$IHC_pos_prop * data$IHC_screen, 0)
  
  # Proportion of individuals with syncytia among the screened ones
  data$IHC_p = data$IHC_n / data$IHC_screen
  
  data$time_post_infection_raw = data$time_post_infection
  data$time_post_infection = as.numeric(as.character(data$time_post_infection))
  data$time_post_infection_approximate =  data$time_post_infection
  data[which(data$time_post_infection_raw == "11-12"), "time_post_infection_approximate"] = 12
  data[which(data$time_post_infection_raw == "16-17"), "time_post_infection_approximate"] = 17
  data[which(data$time_post_infection_raw == "23-24"), "time_post_infection_approximate"] = 24
  data[which(data$time_post_infection_raw == "5-6"), "time_post_infection_approximate"] = 6
  data[which(data$time_post_infection_raw == "3-4"), "time_post_infection_approximate"] = 4
  data[which(data$time_post_infection_raw == "5,7"), "time_post_infection_approximate"] = 6
  data[which(data$time_post_infection_raw == "5-7"), "time_post_infection_approximate"] = 6
  data[which(data$time_post_infection_raw == "9, 11"), "time_post_infection_approximate"] = 10
  data[which(data$time_post_infection_raw == "9-10"), "time_post_infection_approximate"] = 10
  data[which(data$time_post_infection_raw == "8-10"), "time_post_infection_approximate"] = 9
  data[which(data$time_post_infection_raw == "7-8"), "time_post_infection_approximate"] = 8
  data[which(data$time_post_infection_raw == "6-8"), "time_post_infection_approximate"] = 7
  
  # Generate a bunch of unique identifiant to easily remove duplicates
  data$sample_id_tissue_time = paste0(data$id, "_", data$tissue, "_", data$time_post_infection)
  data$sample_id_tissuegroup_time = paste0(data$id, "_", data$tissue_group, "_", data$time_post_infection)
  data$sample_id_tissue_time_approximate = paste0(data$id, "_", data$tissue, "_", data$time_post_infection_approximate)
  data$sample_id_tissuegroup_time_approximate = paste0(data$id, "_", data$tissue_group, "_", data$time_post_infection_approximate)
  data$sample_id_tissue = paste0(data$id, "_", data$tissue)
  data$sample_id_tissuegroup = paste0(data$id, "_", data$tissue_group)
  
  data$histo_score = as.numeric(as.character(data$histo_score))
  
  return(data)
  
}




