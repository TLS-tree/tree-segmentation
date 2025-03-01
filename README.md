## Segmentation of Individual Trees in TLS Point Clouds via Graph Optimization

1. Run *trunk_localization.m* to achieve trunk position. （Due to GitHub's file size limitations, only a small portion of the plot has been extracted as sample data for demonstration.）
2. Run *tree_segmentation.m* to implement individual tree segmentation.

The supervoxel segmentation process relies on the l0 cut-pursuit algorithm, please citing:

- Landrieu L, Obozinski G. Cut pursuit: Fast algorithms to learn piecewise constant functions on general weighted graphs[J]. SIAM Journal on Imaging Sciences, 2017, 10(4): 1724-1766.
- Xi Z, Hopkinson C. 3D graph-based individual-tree isolation (Treeiso) from terrestrial laser scanning point clouds[J]. Remote Sensing, 2022, 14(23): 6116.


## Acknowledgement
This repository contains code from Xi's repository on artemis_treeiso (https://github.com/truebelief/artemis_treeiso), and Wang's repository on SSSC (https://github.com/dwang520/SSSC).


