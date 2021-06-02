#' find_homologs_MDPlot
#' A wrapper for homol.search and plothomol from the nontarget package with explorative output using plotly. Default values are for PFAS search.
#' 
#' @param file a csv file containing 'mz', 'intensity' and 'rt' as variable names. 'rt' should be in minutes.
#' @param plotdefect logical. Whether or not to plot mass defect (mz-round(mz)) instead of rt. Defaults to FALSE.
#' @param intensity variable name. If multiple sample intensities are in the dataframe, then you can choose the specific sample intensity. Defaults to 'intensity'. 
#' @param p_elements character vector. FALSE or chemical elements in the changing units of the homologue series, e.g. c("C","H") for alkane chains or c("C", "F") for perfluorinated compounds. Used to restrict search.Elements to include in the homolog search. Defaults to: c("C", "H", "O")
#' @param p_use_C logical. For elements: take element ratio to C-atoms into account? Used to restrict search
#' @param p_minmz Defines the lower limit of the m/z window to search homologue series peaks, relative to the m/z of the one peak to search from. Absolute m/z value [u].
#' @param p_maxmz Defines the upper limit of the m/z window to search homologue series peaks, relative to the m/z of the one peak to search from. Absolute m/z value [u].
#' @param p_minrt Defines the lower limit of the retention time (RT) window to look for other homologue peaks, relative to the RT of the one peak to search from, i.e., RT+minrt. For decreasing RT with increasing HS mass, use negative values of minrt.
#' @param p_maxrt Defines the upper limit of the RT window to look for other homologue peaks, relative to the RT of the one peak to search from, i.e., RT+maxrt. See minrt.
#' @param p_ppm Should mztol be set in ppm (TRUE) or in absolute m/z [u] (FALSE)?
#' @param p_mztol m/z tolerance setting: +/- value by which the m/z of a peak may vary from its expected value. If parameter ppm=TRUE (see below) given in ppm, otherwise, if ppm=FALSE, in absolute m/z [u].
#' @param p_rttol Retention time (RT) tolerance by which the RT between two adjacent pairs of a homologue series is allowed to differ. Units as given in column 3 of peaklist argument, e.g. [min].
#' @param p_minlength Minimum number of peaks in a homologue series.
#' @param p_mzfilter Vector of numerics to filter for homologue series with specific m/z differences of their repeating units, given the tolerances in mztol. Mind charge z!
#' @param p_vec_size Vector size. Ignore unless a relevant error message is printed (then try to increase size).
#' @param p_mat_size Matrix size for recombining, multiple of input tuples. Ignore unless a relevant error message is printed (then try to increase size).
#' @param p_R2 FALSE or 0<numeric<=1. Coefficient of determination for cubic smoothing spline fits of m/z versus retention time; homologue series with lower R2 are rejected. See smooth.spline.
#' @param p_spar Smoothing parameter, typically (but not necessarily) in (0,1]. See smooth.spline.
#' @param p_plotit Logical FALSE or 0<integer<5. Intermediate plots of nearest neighbor paths, spline fits of individual homologues series >=minlength, clustered HS pairs, etc .
#' @param p_deb Debug returns, ignore.
#'
#' @return
#' @export
#'
#' @examples
#' find_homologs("./data/LCneg_ComponentsRC.csv", 
#' plotdefect = FALSE,
#' intensity = N_rep1,
#' p_elements=c("C","H", "O", "S", "F"),
#' p_use_C=FALSE,
#' p_minmz=49.9,
#' p_maxmz=50,
#' p_minrt=0.5,
#' p_maxrt=2,
#' p_ppm=TRUE,
#' p_mztol=10,
#' p_rttol=0.5,
#' p_minlength = 3,
#' p_mzfilter=FALSE,
#' p_vec_size=3E6,
#' p_mat_size=3,
#' p_R2=.98,
#' p_spar=.45,
#' p_plotit=FALSE,
#' p_deb=0)




find_homologs_MDPlot <- function(file, 
                          intensity = intensity,
                          p_elements=c("C","H", "O", "S", "F", "N"),
                          p_use_C=FALSE,
                          p_minmz=49.9,
                          p_maxmz=50,
                          p_minrt=0.5,
                          p_maxrt=2,
                          p_ppm=TRUE,
                          p_mztol=12,
                          p_rttol=0.5,
                          p_minlength = 3,
                          p_mzfilter=FALSE,
                          p_vec_size=3E6,
                          p_mat_size=3,
                          p_R2=.98,
                          p_spar=.45,
                          p_plotit=FALSE,
                          p_deb=0) {
  
  
  
  p_isotopes <- read.csv("./data/isotopes.csv")
  
  df <- file
  
  df <- df %>%
    dplyr::select(mz, {{intensity}}, rt) %>% #uses embrace to specify sample intensity
    dplyr::rename(mass = mz) %>%
    dplyr::mutate(rt = round(rt/60, 2)) %>% # convert rt in file from sec to min
    dplyr::mutate(intensity = as.integer({{intensity}}), .keep = "unused") %>% 
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
                        name = ~`m/z increment`) %>%
      plotly::layout(
        yaxis = list(
          title = "m/z"),
        xaxis = list(
          title = "Retention time")
      )
    
  } 


