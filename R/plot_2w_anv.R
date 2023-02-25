
#' Plotting the results of two way anova
#'
#' This function allows you to plot the results of two-way anova.
#' @param two_ways_res The results of the two way anova.
#' @param first_con the first condition for the comparisson
#' @param second_con the second condition for the comparisson
#'
#' @return a plot depicting the results of the two way anova
#'
#' @examples
#' \dontrun{
#' print(1)
#' }
#'
#' @export


plot_2w_anv <- function(two_ways_res, first_con = 'drosophila',
                        second_con = 'mouse') 
{

library(ggplot2)
  
  
  ##log10 transformation of pvalues-two way
  new_log10<-data.frame()
  for (i in 1:dim(two_ways_res)[1]){
    
    val<-as.numeric(two_ways_res[i,2])
    
    new_log10<-rbind(new_log10,-log10(val))
  }
  
  y_val<--log10( as.numeric(two_ways_res$P_value))

  #plotting  
two_way_p<-ggplot(two_ways_res, aes(x = two_ways_res$Metabolites, y = y_val, color = two_ways_res$Significance))+geom_point()+
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())



two_way_p + labs(color = "Significance")+geom_hline(yintercept=1.3,color = "black")+xlab("Metabolite") + ylab("-log10 P-value")



img_n <- paste0('ano_', first_con, '_vs_', second_con, ".png")
ggsave(img_n)
return(two_way_p)

}
