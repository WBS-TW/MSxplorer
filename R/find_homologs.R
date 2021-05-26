

# Using homol.search in "nontarget" package


find_homologs <- function(file, 
                          plotdefect = FALSE,
                          p_isotopes = p_isotopes,
                          p_elements=c("C","H", "F", "O", "S"),
                          p_use_C=FALSE,
                          p_minmz=49.96,
                          p_maxmz=49.99,
                          p_minrt=1,
                          p_maxrt=4,
                          p_ppm=TRUE,
                          p_mztol=10,
                          p_rttol=0.5,
                          p_minlength = 3,
                          p_mzfilter=FALSE,
                          p_vec_size=3E6,
                          p_mat_size=3,
                          p_R2=.98,
                          p_spar=.45,
                          p_plotit=FALSE,
                          p_deb=0) {

data("isotopes")
p_isotopes <- isotopes


df <- vroom::vroom(file)

df <- df %>%
  dplyr::select(mz, intensity, rt) %>% #Choose the correct intensity sample!
  dplyr::rename(mass = mz) %>%
  dplyr::mutate(rt = round(rt/60, 2)) %>%
  dplyr::mutate(intensity = as.integer(intensity), .keep = "unused") %>% #Need to rename each time, put in function!
  dplyr::select(mass, intensity, rt) %>%
  dplyr::filter(intensity > 0) %>%
  dplyr::arrange(rt) %>%
  as.data.frame() #need to convert to dataframe, otherwise nontarget cannot recognize

res_homologs <-  nontarget::homol.search(peaklist = df,
                             isotopes = p_isotopes,
                             elements= p_elements,
                             use_C=p_use_C,
                             minmz=p_minmz,
                             maxmz=p_maxmz,
                             minrt=p_minrt,
                             maxrt=p_maxrt,
                             ppm=p_ppm,
                             mztol=p_mztol,
                             rttol=p_rttol,
                             minlength = p_minlength,
                             mzfilter=p_mzfilter,
                             vec_size=p_vec_size,
                             mat_size=p_mat_size,
                             R2=p_R2,
                             spar=p_spar,
                             plotit=p_plotit,
                             deb=p_deb)

if (plotdefect == FALSE) {
# Using plotly
res_homologs[[1]] %>%
    dplyr::filter(`HS IDs` > 0) %>%
    dplyr::mutate(`m/z increment` = round(as.numeric(`m/z increment`), 3)) %>%
    dplyr::mutate(`m/z increment` = as.factor(`m/z increment`)) %>%
  plotly::plot_ly() %>%
    plotly::add_markers(x = df$rt, 
              y = df$mass, 
              type = "scatter", 
              name = "All peaks",
              opacity = 0.9,
              colors = "BrBG",
              marker = list(color = "lightgrey")) %>%
    plotly::add_trace(x = ~RT, 
            y = ~mz, 
            color = ~`HS IDs`,
            type = "scatter",
            mode = "lines+markers",
            name = ~`m/z increment`)

} else {
# Plot mass defect
res_homologs[[1]] %>%
    dplyr::filter(`HS IDs` > 0) %>%
    dplyr::mutate(`m/z increment` = round(as.numeric(`m/z increment`), 3)) %>%
    dplyr::mutate(`m/z increment` = as.factor(`m/z increment`)) %>%
    plotly::plot_ly() %>%
    plotly::add_markers(x = df$mass, 
              y = (df$mass - round(df$mass)), 
              type = "scatter", 
              name = "All peaks",
              opacity = 0.9,
              colors = "BrBG",
              marker = list(color = "lightgrey")) %>%
    plotly::add_trace(x = ~mz, 
            y = ~(mz-round(mz)), 
            color = ~`HS IDs`,
            type = "scatter",
            mode = "lines+markers",
            name = ~`m/z increment`)

}
}

