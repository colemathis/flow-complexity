require(tidyverse)
require(MASS)
source("../scripts/assembly-map.r")

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


get_AI_stats_from_sim <- function(sim_number){
    # Get the mean, median and max AI from a simulation 
    # sum over all reactors before fitting 
    ts_data <- read.csv(paste0(sim_number_to_directory(sim_number),"timeseries.csv"))
    average_data <- sum_over_reactors(ts_data)
    # Get only the last time point TODO: come back and decide if averaging here makes sense
    max_time <- max(average_data$time)
    average_data <- average_data %>% 
                        filter(total_value > 1) %>% 
                        group_by(variable) %>% 
                        summarise(average_value = sum(total_value)/max_time)
    average_data$AI <- assembly_index(average_data$variable)
    # Group all compounds by assembly index
    total_AI_average <- average_data %>% group_by(AI) %>% summarise(total_average = sum(average_value))
    
    max.AI <- max(total_AI_average$AI)
    median.AI <- median(total_AI_average$AI)
    total_AI_average$normed_value = total_AI_average$total_average/sum(total_AI_average$total_average)
    weighted.mean.AI <- sum(total_AI_average$normed_value * total_AI_average$AI)
    
    results = list(max.AI = max.AI, 
                   median.AI = median.AI,
                   weighted.mean.AI = weighted.mean.AI)

    return(results)

}


get_poisson_fit_from_sim <- function(sim_number, by_reactor=TRUE){
    # Get the mean of a poisson fit and estimate of confidence from MLE, and normalization
    # sum over all reactors before fitting 
    ts_data <- read.csv(paste0(sim_number_to_directory(sim_number),"timeseries.csv"))
    ts_data <- ts_data %>% filter(variable > 1)
    mean_data <- ts_data %>% group_by(reactor, variable) %>% summarise(total.v = sum(value))
    mean_data$AI <- assembly_index(mean_data$variable)
    coarse_grainned <- mean_data %>% 
                            group_by(reactor, AI) %>%
                            summarise(total = sum(total.v))
    
    
    if (by_reactor){
        reactors = unique(coarse_grainned$reactor)
        fit.by.reactor <- data.frame(reactor = c(), lambda = c(), lambda.minus = c(), lambda.plus = c())
        
        for (r in reactors){
            samples <- c()
            reactor.data <- coarse_grainned %>% filter(reactor == r)
            AI.range <- unique(reactor.data$AI)
            for (ai in AI.range){
                ai.copy.number <- reactor.data[reactor.data$AI == ai,]$total
                samples <- c(samples, rep(ai, ai.copy.number))
            }
            p.mle <- fitdistr(samples, densfun="poisson")
            lambda <- p.mle[[1]]
            confidence_intervals <- confint(p.mle)
            lambda.minus <- confidence_intervals[1,1]
            lambda.plus <- confidence_intervals[1,2]
            # print(lambda, lambda.minus, lambda.plus)
            this.df <- data.frame(reactor = r, 
                                  lambda = lambda,
                                  lambda.minus = lambda.minus,
                                  lambda.plus = lambda.plus)
            fit.by.reactor <- rbind(fit.by.reactor, this.df)
        }
        rownames(fit.by.reactor) <- 1:nrow(fit.by.reactor)
    }
    
    # if (nrow(total_AI_counts) > 3) {
    #     # If there's more than 3 assembly numbers present
    #     fit.model <- glm(log(total_count) ~ AI, data = total_AI_counts)
    #     # Extract fit coefficient and R^2
    #     alpha <- as.numeric(fit.model$coefficients["AI"])
    #     r2 <- r2glm(fit.model)
    #     # Calculate confidence intervals
    #     confidence_intervals <- confint(fit.model)
    #     alpha.minus <- confidence_intervals[2,1]
    #     alpha.plus <- confidence_intervals[2,2]
    # }
    # # If there's not enough data return NA
    # else{
    #     alpha <- NA
    #     alpha.minus <- NA
    #     alpha.plus <- NA
    #     r2 <- NA
        
    # }

    # results = list(alpha = alpha, 
    #                alpha.minus = alpha.minus,
    #                alpha.plus = alpha.plus,
    #                r2 = r2)
    return(fit.by.reactor)

}