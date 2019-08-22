
library(CompSign)
data("Breast560")
data("Breast560")
exp <- count_matrix(Breast560_count)
exp <- count_matrix(Breast560)
graph_binary_joint_presence(exp, graph = TRUE)
graph_binary_joint_presence(exp, heatmap = TRUE)
graph_binary_joint_presence(exp, graph = TRUE, thres = 0.1)
