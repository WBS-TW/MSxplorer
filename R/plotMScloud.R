library(tidyverse)
library(viridis)

### 2. source geom_flat_violin and make raincloud theme

# source the geom flat violin function for ggplot2
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

####
# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )
####

raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  #legend.position = "none",
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
)


### 3. Read in data

# from all feature groups
clearCache("all")
mz_data_all <- patRoon::filter(fGroupsXCMS,
                               #preAbsMinIntensity = 5000,
                               #mzRange = c(70, 1200),
                               #retentionRange = c(35, 940),
                               relMinReplicateAbundance = NULL,
                               absMinIntensity = 10000,  # 
                               maxReplicateIntRSD = NULL,
                               blankThreshold = 5, # Soil and Rwref is filtered out since PFOS was in 2KBland sample
                               removeBlanks = TRUE)
mz_data_all <- as.data.table(mz_data_all)
# from suspect feature groups
# Filter the fgroups but keep blankThreshold = NULL otherwise the filter will remove these from the samples

mz_data_susp <- patRoon::filter(fGroupsScrXCMS,
                                #preAbsMinIntensity = 500,
                                #mzRange = c(70,1200),
                                #retentionRange = c(35, 940),
                                relMinReplicateAbundance = NULL,
                                absMinIntensity = 10000,
                                maxReplicateIntRSD = NULL,
                                blankThreshold = NULL,
                                removeBlanks = FALSE)


scr <- screenSuspects(mz_data_susp, suspects, mzWindow = 0.01)
mz_data_susp <- groupFeaturesScreening(mz_data_susp, scr)
mz_data_susp <- patRoon::filter(mz_data_susp,
                                blankThreshold = 5,
                                removeBlanks = TRUE)
mz_data_susp <- as.data.frame(mz_data_susp)

#colnames(my_data)
#glimpse(my_data)


### 4. Tidy data to long format
# function for intensity intensities 
normalit<-function(m){
  (m - min(m))/(max(m)-min(m))
}


groups <- mz_data_all %>%
  select(-group, -ret, -mz) %>%
  colnames() %>%
  enframe() %>%
  rename(Sample = value) %>%
  select(-name) %>%
  mutate(Group = case_when(
    str_detect(Sample, "water") ~ "Water",
    str_detect(Sample, "Insect") ~ "Insect",
    str_detect(Sample, "larvae") ~ "Larvae",
    str_detect(Sample, "Blank") ~ "Blank",
    str_detect(Sample, "terr") ~ "Terr_inv",
    str_detect(Sample, "wref") ~ "Water",
    str_detect(Sample, "Soil") ~ "Soil",
    str_detect(Sample, "Sed") ~ "Sediment",
    str_detect(Sample, "STD") ~ "Standard",
    str_detect(Sample, "plant") ~ "Plant"
  ))



mz_data_all_l <- mz_data_all %>%
  # mutate(mz = round(mz, 3)) %>%
  select(-c(group, ret)) %>%
  pivot_longer(-mz, names_to = "Sample", values_to = "Intensity") %>%
  dplyr::filter(Intensity != 0) 
# %>%
#   group_by(Sample) %>%
#   mutate(Intensity = normalit(Intensity)) %>%  # use to normalize intensities per sample
#   ungroup()

mz_data_susp_l <- mz_data_susp %>%
  # mutate(mz = round(mz, 3)) %>%
  select(-c(group, ret)) %>%
  pivot_longer(-mz, names_to = "Sample", values_to = "Intensity") %>%
  dplyr::filter(Intensity != 0) 
# %>%
#   group_by(Sample) %>%
#   mutate(Intensity = normalit(Intensity)) %>%  # use to normalize intensities per sample
#   ungroup()

# joining both peak groups into one data frame, add groups and order by factor, blanks (if not previously filtered out) 
# were removed here since they were filtered out (zero values) in the mz_data_susp_l
mz_all_susp <- mz_data_all_l %>%
  full_join(mz_data_susp_l, by = c("mz", "Sample")) %>%
  left_join(groups, by = "Sample") 


mz_all_susp$Group <- factor(mz_all_susp$Group, levels = c("Blank", "Standard",  "Soil", "Sediment", "Water", "Plant", "Insect", "Larvae", "Terr_inv"))


mz_all_susp_group <- mz_all_susp %>%
  group_by(Group, mz) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(Intensity.x = round(Intensity.x, 0),
         Intensity.y = round(Intensity.y, 0))



### 5. Cloud plot by samples 

# all + susp
ggplot(data = mz_all, 
       aes(x =  fct_reorder(Sample, desc(Group)), y = mz, fill = Group)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  #geom_point(aes(y = mz, color = Group), position = position_jitter(width = .1), alpha = 0.4) +
  #geom_point(aes(y = mz, size = Intensity.x), color = "black",  position = position_jitter(width = .1), alpha = 0.9) +
  geom_point(aes(y = mz, color = Group, size = Intensity.x), position = position_jitter(width = .1), alpha = 0.3) +
  geom_point(aes(y = mz, size = Intensity.y, alpha = 0.8), color = "red") +
  #geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "D") +
  # scale_fill_viridis(discrete = TRUE) +
  scale_color_brewer(palette = "Pastel2") +
  scale_fill_brewer(palette = "Pastel2") +
  scale_y_continuous(breaks = seq(0,2000, 100)) +
  coord_flip() + # flip or not
  theme_bw() +
  raincloud_theme



### 5-2. Raincloud plot by group

# all + susp + group
ggplot(data = mz_all_susp_group, # use "mz_all_susp_group %>% filter(mz < 1300)" to filter to crop axis, need to change the breaks
       aes(x =  Group, y = mz, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  #geom_point(aes(y = mz, color = Group), position = position_jitter(width = .1), alpha = 0.4) +
  geom_point(aes(y = mz, color = Group, size = Intensity.x), position = position_jitter(width = .1), alpha = 0.4) +
  geom_point(aes(y = mz, size = Intensity.y, alpha = 0.8), color = "red") +
  
  #geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "D")+
  # scale_fill_viridis(discrete = TRUE) +
  scale_color_brewer(palette = "Pastel2") +
  scale_fill_brewer(palette = "Pastel2") +
  scale_y_continuous(breaks = seq(0,1200, 100)) +
  coord_flip() + # flip or not
  theme_bw() +
  raincloud_theme


# Filter mz <1300
ggplot(data = mz_all_susp_group %>% dplyr::filter(mz < 1300),
       aes(x =  Group, y = mz, fill = Group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = mz, color = Group), position = position_jitter(width = .1), alpha = 0.4) +
  geom_point(aes(y = mz, size = Intensity.y, alpha = 0.8), color = "red") +
  
  #geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "D")+
  # scale_fill_viridis(discrete = TRUE) +
  scale_color_brewer(palette = "Pastel2") +
  scale_fill_brewer(palette = "Pastel2") +
  scale_y_continuous(breaks = seq(0,1200, 100)) +
  coord_flip() + # flip or not
  theme_bw() +
  raincloud_theme


### 7. Calculate summary stats 

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

mz_sum_all_ld <- mz_data_all_l %>%
  select(-mz) %>%
  group_by(Sample) %>%
  summarise_all(funs(mean, median, lower = lb, upper = ub))

mz_sum_susp_ld <- mz_data_susp_l %>%
  select(-mz) %>%
  group_by(Sample) %>%
  summarise_all(funs(mean, median, lower = lb, upper = ub))



### 8. Same plot - replacing the boxplot with a mean and confidence interval 

# note: this gives a Warning: Ignoring unknown aesthetics: y
# It is from the geom_errorbar, but when removed gives and error
# Error in FUN(X[[i]], ...) : object 'Sensitivity' not found


mz_g <- 
  ggplot(data = mz_data_all_l, 
         aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Intensity, color = Sample), 
             position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_point(data = mz_sum_all_ld, aes(x = Sample, y = mean), 
             position = position_nudge(x = 0.3), size = 2.5) +
  geom_errorbar(data = mz_sumld, aes(ymin = lower, ymax = upper, y = mean), 
                position = position_nudge(x = 0.3), width = 0) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  coord_flip() + # flip or not?
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  raincloud_theme
mz_g