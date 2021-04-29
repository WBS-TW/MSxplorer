
# # NOT WORKING YET
# library(readxl)
# library(tidyverse)
# library(jsonlite)
# library(ChemmineR)
# library(classyfireR)
# library(writexl)
# 
# pos_neg_comb <- read_xlsx("D:/Projects/Naturvardsverket/SSA-NTA Damm/Experiment/Data_analysis/NORMAN suspect Dust LC 28 on 39.xlsm",
#                           sheet = "pos_neg_remove_dup")
# 
# pos_neg_comb_dupl <- pos_neg_comb %>%
#   group_by(Formula) %>%
#   tally() %>%
#   mutate(weight = 1/n)
# 
# pos_neg_comb <- pos_neg_comb %>%
#   left_join(pos_neg_comb_dupl, by = "Formula")
# 
# 
# norman_LC <- read_xlsx("D:/Projects/Naturvardsverket/SSA-NTA Damm/Experiment/Data_analysis/NORMAN suspect list dust LC.xlsx",
#                        sheet = "INDOORCT16_merged")
# 
# # Using classyfireR to group hit list
# 
# inchi_keys <- pos_neg_comb$InChIKey
# 
# res <- tibble(id = NULL, Level = NULL, Classification = NULL, CHEMONT = NULL)
# 
# classify <- get_classification(inchi_keys[1])
# class_df <- classify@classification
# class_df <- class_df %>%
#   mutate(id = paste(pos_neg_comb[1,1])) %>%
#   select(id, everything())
# res <- rbind(class_df)
# 
# for (i in 2:length(inchi_keys)) {
#   tryCatch({
#     classify <- get_classification(inchi_keys[i])
#     class_df <- classify@classification
#     class_df <- class_df %>%
#       mutate(id = paste(pos_neg_comb[i,1])) %>%
#       select(id, everything())
#     res <- res %>%
#       full_join(class_df)
#     Sys.sleep(6) # Classyfire API restricts requests to 12 per minute
#   },error = function(e){})
# }
# 
# #  grouping data
# 
# res %>%
#   group_by(Level) %>%
#   tally()
# 
# res2 <- res %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_class <- res %>%
#   filter(Level == "class") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_subclass <- res %>%
#   filter(Level == "subclass") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_lvl5 <- res %>%
#   filter(Level == "level 5") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_lvl6 <- res %>%
#   filter(Level == "level 6") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_lvl7 <- res %>%
#   filter(Level == "level 7") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# res_lvl8 <- res %>%
#   filter(Level == "level 8") %>%
#   left_join(pos_neg_comb, by = "id")
# 
# # write to excel
# write_xlsx(list(all_results = res2,
#                 class = res_class,
#                 subclass = res_subclass,
#                 level5 = res_lvl5,
#                 level6 = res_lvl6,
#                 level7 = res_lvl7,
#                 level8 = res_lvl8),
#            path = "D:/Projects/Naturvardsverket/SSA-NTA Damm/Experiment/Data_analysis/NORMAN suspect Dust LC 28 on 39_classyfire_grouped.xlsx")
# 
