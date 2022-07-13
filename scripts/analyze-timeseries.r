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

get_exp_fit_from_sim <- function(sim_number){
    # Get the slope, confidence interval and simulated R2 from an exponential fit like 
    # total_count ~ exp(AI) (or log(total_count) ~ AI)
    # sum over all reactors before fitting 
    ts_data <- read.csv(paste0(sim_number_to_directory(sim_number),"timeseries.csv"))
    average_data <- sum_over_reactors(ts_data)
    average_data <- average_data %>% filter(total_value > 1)
    # Get only the last time point TODO: come back and decide if averaging here makes sense
    max_time <- max(average_data$time)
    end_time <- average_data %>% filter(time == max_time)
    end_time$AI <- assembly_index(end_time$variable)
    # Group all compounds by assembly index
    total_AI_counts <- end_time %>% group_by(AI) %>% summarise(total_count = sum(total_value))
    if (nrow(total_AI_counts) > 3) {
        # If there's more than 3 assembly numbers present
        fit.model <- glm(log(total_count) ~ AI, data = total_AI_counts)
        # Extract fit coefficient and R^2
        alpha <- as.numeric(fit.model$coefficients["AI"])
        r2 <- r2glm(fit.model)
        # Calculate confidence intervals
        confidence_intervals <- confint(fit.model)
        alpha.minus <- confidence_intervals[2,1]
        alpha.plus <- confidence_intervals[2,2]
    }
    # If there's not enough data return NA
    else{
        alpha <- NA
        alpha.minus <- NA
        alpha.plus <- NA
        r2 <- NA
        
    }

    results = list(alpha = alpha, 
                   alpha.minus = alpha.minus,
                   alpha.plus = alpha.plus,
                   r2 = r2)
    return(results)

}