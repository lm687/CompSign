## give a graph of the joint presence

graph_binary_joint_presence <- function(exposures, thres=0, title='', heatmap=FALSE, graph=FALSE){
  require(qgraph)

    bin_presence <- apply(exposures, 1, function(x) x>thres)
    bin_presence_outer <- outer(1:nrow(bin_presence), 1:nrow(bin_presence),
                                Vectorize(function(x,y) sum(apply(bin_presence[c(x,y),], 2, sum) == 2)/sum(bin_presence[x,])))
    rownames(bin_presence_outer) <- colnames(bin_presence_outer) <- colnames(exposures)

    order_decreasing <- 1:ncol(bin_presence_outer)
    if(heatmap){
      print(ComplexHeatmap::Heatmap(bin_presence_outer[order_decreasing,order_decreasing],
                                    cluster_rows = FALSE, cluster_columns = FALSE,
                                    column_title = paste0('Fraction of overlap in presence ', title)))
    }
    if(graph){
      qgraph(bin_presence_outer, layout='circular', minimum=0.3, theme='colorblind') ## a bit difficult to see
    }
}

