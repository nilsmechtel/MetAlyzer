##new two ways

#' Two way anova
#'
#' This function allows you to perform two-way anova.
#' @param aggregate_obj The processed data.
#' @param categorical The categories, such as Species or Tissue_type that are providing groups for anova test
#' @param log2_transformation Type of transformation, no transformation: no, log2: l2, log10: l10
#' @param impute If TRUE, zero imputation will be performed, otherwise, no imputation
#' @param first_con the first condition for the comparisson
#' @param second_con the second condition for the comparisson
#'
#' @return A data frame containing the results of the two way anova
#'
#' @examples
#' \dontrun{
#' print(1)
#' }


do_anova_2w <- function(aggregate_obj,             # The processed data
                        categorical,
                        log2_transformation = 'TRUE',   # Type of transformation, no transformation: no, log2: l2, log10: l10
                        impute = 'FALSE',            # If TRUE, zero imputation will be performed, otherwise, no imputation
                        first_con = 'drosophila',
                        second_con = 'mouse') 
{



  column_name <- rlang::sym(categorical)

two_sets <- aggregate_obj %>%
  filter(UQ(column_name) %in% c(first_con, second_con))




#take the unique list of metabolites
un_mets <- unique(two_sets$Metabolite)


res=list()
res_all=list()
two_ways_res<-data.frame(matrix(NA,    # Create empty data frame
                                nrow = length(un_mets),
                                ncol = 3))


#calculating two way anova for each metabolite
l=1
for (met in un_mets){
  
  new_list <- two_sets %>%
    filter(Metabolite %in% c(met))
  
  if (log2_transformation == TRUE){
    res.aov2 <- aov(new_list$log2_Conc ~ new_list[[categorical]] , data = new_list)
  }
  else{ if (impute == TRUE){
          res.aov2 <- aov(new_list$imputed_Conc ~ new_list[[categorical]] , data = new_list)
  }
        else{
          res.aov2 <- aov(new_list$Concentration ~ new_list[[categorical]] , data = new_list)
        }
  }
  sum_res <- summary(res.aov2)
  p_val <- sum_res[[1]][["Pr(>F)"]]
  if (p_val < 0.05){
    res<-c(met, p_val[1], 'TRUE')
  }
  else{
    res<-c(met, p_val[1], 'FALSE')
  }
  res_all<-t(as.data.frame(res))
  two_ways_res[l,]<-as.data.frame(res_all)
  l = l+1
}

two_ways_res <- two_ways_res %>%
  rename(Metabolites = X1,
         P_value = X2,
         Significance = X3)


#exporting the results

dataList_xlsx <- paste0('anova_', first_con, '_vs_', second_con, ".xlsx")
write.xlsx(two_ways_res,dataList_xlsx)

return(two_ways_res)
}