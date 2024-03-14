
#' Calculates the non-isothermal Kovats retention indices (from temperature-programming, using definition of Van den Dool and Kratz: Ix = 100n + 100(tx-tn) / (tn+1 − tn)
#'
#' @param data data.frame with name and rt (in minutes)
#' @param alkaneSeries data.frame with Num (number of carbons in of the alkane) and rt (in minutes) 
#'
#' @return
#' @export
#'
#' @examples
calculateRI <- function(data, alkaneSeries) {
  
######Test data
  data <- data.frame(name = c("UnknownID=0|RT=4.55|RI=1076.5", "UnknownID=192|RT=22.944|RI=2329.8"),  rt = c(4.55, 22.944))
  alkaneSeries <- data.frame(Num = c(11, 13, 15, 17, 19, 21, 23, 25),
                             rt = c(4.90, 8.00, 11.13, 14.05, 16.70, 19.37, 22.45, 25.77)
                             )
###########
    

  m <- dplyr::case_when(
    data$rt[i] >= alkaneSeries$rt[1] & data$rt[i] < alkaneSeries$rt[2] ~ alkaneSeries$Num[1],
    data$rt[i] >= alkaneSeries$rt[2] & data$rt[i] < alkaneSeries$rt[3] ~ alkaneSeries$Num[2],
    data$rt[i] >= alkaneSeries$rt[3] & data$rt[i] < alkaneSeries$rt[4] ~ alkaneSeries$Num[3],
    data$rt[i] >= alkaneSeries$rt[4] & data$rt[i] < alkaneSeries$rt[5] ~ alkaneSeries$Num[4],
    data$rt[i] >= alkaneSeries$rt[5] & data$rt[i] < alkaneSeries$rt[6] ~ alkaneSeries$Num[5],
    data$rt[i] >= alkaneSeries$rt[6] & data$rt[i] < alkaneSeries$rt[7] ~ alkaneSeries$Num[6],
    data$rt[i] >= alkaneSeries$rt[7] & data$rt[i] <= alkaneSeries$rt[8] ~ alkaneSeries$Num[7]
    )
  
  n <- dplyr::case_when(
    data$rt[i] >= alkaneSeries$rt[1] & data$rt[i] < alkaneSeries$rt[2] ~ alkaneSeries$Num[2],
    data$rt[i] >= alkaneSeries$rt[2] & data$rt[i] < alkaneSeries$rt[3] ~ alkaneSeries$Num[3],
    data$rt[i] >= alkaneSeries$rt[3] & data$rt[i] < alkaneSeries$rt[4] ~ alkaneSeries$Num[4],
    data$rt[i] >= alkaneSeries$rt[4] & data$rt[i] < alkaneSeries$rt[5] ~ alkaneSeries$Num[5],
    data$rt[i] >= alkaneSeries$rt[5] & data$rt[i] < alkaneSeries$rt[6] ~ alkaneSeries$Num[6],
    data$rt[i] >= alkaneSeries$rt[6] & data$rt[i] < alkaneSeries$rt[7] ~ alkaneSeries$Num[7],
    data$rt[i] >= alkaneSeries$rt[7] & data$rt[i] <= alkaneSeries$rt[8] ~ alkaneSeries$Num[8]
  )
  
  
  RI <- 100*n + 100*(m-n) * (data$rt[i] - alkaneSeries[alkaneSeries$Num == n,]$rt)/(alkaneSeries[alkaneSeries$Num == m,]$rt - alkaneSeries[alkaneSeries$Num == n,]$rt)

      
     

  
}