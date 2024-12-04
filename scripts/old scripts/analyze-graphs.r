require(igraph)
require(tidyverse)



graph_from_sim_num <- function(sim_num){

    sim_dir = paste0("../data/sims/",sim_num, "/")
    graph_path = paste0(sim_dir, "graph.csv")
    graph_df = read.csv(graph_path)
    source_verts <- graph_df %>% filter(source_inflow == "true") # %>% select("sources") # nolint
    source_verts <- unique(source_verts$sources)
    all_verts <- unique(c(graph_df$sources, graph_df$destinations))
    vert_df = data.frame(reactor_id = all_verts)
    vert_df$source_node = vert_df$reactor_id %in% source_verts
    
    reduced_graph_df <- graph_df[c("sources", "destinations")]
    graph <- graph_from_data_frame(reduced_graph_df, directed=TRUE, vertices = vert_df) # nolint
    return(graph)

}

mean_dist_from_source <- function(graph, reactor_id){

    source_nodes <- V(graph)[vertex.attributes(graph)$source_node]
    path.lengths <- c()
    for (s in source_nodes){
        sp <- shortest_paths(graph, s, reactor_id, mode= "out")
        path.lengths <- c(path.lengths, length(sp$vpath[[1]]) - 1)
    }
    d <- mean(path.lengths)
    return(d)

}