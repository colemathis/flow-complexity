#require(tidyverse)

map_df_name <- "Assembly-10000.csv"

check_data <- function(fname){
    data_exists = FALSE 
    if (file.exists(fname)){

        df = read.csv(fname)
        n <- names(df)
        nr <- nrow(df)
        if ("integer" %in% n && "assembly.index" %in% n && nr > 0){
            data_exists = TRUE
        }
    }
    return(data_exists)
}



get_assembly_index <- function(my_int, map_df) {
    AI = NaN
    sub_df <- map_df %>% filter(integer == my_int)
    if (nrow(sub_df) == 1){
        AI = sub_df$assembly.index
    }

    return(AI)
}

vec_assembly_index <- Vectorize(get_assembly_index, "my_int")

assembly_index <- function(my_integers){

    map_df <- read.csv(map_df_name)
    AIs = vec_assembly_index(my_integers, map_df)

    return(AIs)

}

check_data(map_df_name)