require(tidyverse)
source("assembly-map.r")

sim_number_to_directory <- function(sim_num){
    sim_dir = paste0("../data/sims/",sim_num, "/")
    return(sim_dir)
}

sum_over_reactors <- function(timeseries_data){
    average_data <- timeseries_data %>% 
                        group_by(variable,time) %>%
                        summarise(total_value = sum(value))
    return(average_data)
}

r2glm <- function(model) {

  summaryLog <- summary(model)
  1 - summaryLog$deviance / summaryLog$null.deviance

}
