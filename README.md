# Localized-Homology
依赖库：

gudhi  
eagerpy  
POT 

两种topoloss，在topoloss.py中

1，计算点云与ground truth的 wasserstein distance，可以选一段无分叉的血管作为ground truth来维持一个聚类内的连通性，也可以用来约束整个血管的拓扑性质
topoloss_wasserstein_based(pt_gen, pt_gt=[], H0=True, H1=True, max_length_percent=0.3): 
2, 维持连通性，慎用，比较暴力
topoloss_connectivity(pt_gen, max_length_percent=0.3,topk = 3):


