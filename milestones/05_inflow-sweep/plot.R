library(ggplot2)
library(dplyr)
library(gridExtra)

################################################################################

plot_timeseries <- function() {

    timeseries <- read.csv("data/timeseries.csv")

    ts <- timeseries
    ts <- ts[ts$sim_number == 1, ]
    ts <- ts[ts$integer < 10, ]
    
    library(dplyr)
    gts <- ts %>%
        group_by(time, integer) %>%
        summarize(freq_mean = mean(frequency), .groups = 'drop')
    
    library(ggplot2)
    p <- ggplot(gts, aes(x = time,
                         y = freq_mean,
                         color = factor(integer))) +
        geom_line(size = 2) +
        labs(
             ## title = "Time series for sim 1, avg over all chemostats",
             ## x = "Time",
            ## y = "Count",
             color = "Integers") +
        theme(
            ## legend.position = "right",
              panel.background = element_rect(fill = "black"),
              plot.background = element_rect(fill = "black"),
              ## axis.text = element_text(color = "white"),
              ## axis.title = element_text(color = "white"),
              ## plot.title = element_text(color = "white"),
              ## legend.background = element_rect(fill = "black"),
              ## legend.text = element_text(color = "white"),
              ## legend.title = element_text(color = "white"),
              legend.position = "none",
              axis.text = element_blank()
              )
    
    ggsave("figs/time-series_R.pdf", plot = p)
    
}

################################################################################

plot_multi <- function() {

    timeseries <- read.csv("data/timeseries.csv")
    ts_multi <- timeseries %>% filter(sim_number == 1) %>% filter(integer < 10)
    nspecies <- 10
    nchem <- 25
    
    plots <- list()
    
    for (i in 1:nchem) {
        f <- ts_multi %>% filter(chemostat_id == i)
        
        plots[[i]] <- ggplot(f,
                             aes(x = time,
                                 y = frequency,
                                 color = factor(integer))
                             ) + 
            geom_line() +
            ylim(0, 1000) + 
            theme(legend.position = "none") +
            theme(axis.title = element_blank())  # Remove axis titles
    }

    pdf("figs/multi_R.pdf", width = 10, height = 10)
    grid_plot <- grid.arrange(grobs = plots, ncol = 5, nrow = 5)
    dev.off()
    ## ggsave("figs/multi_R.pdf", plot = grid_plot, width = 10, height = 10)

}

################################################################################

plot_timeseries()
plot_multi()
