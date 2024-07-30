# MHDR: Improving Drug Repositioning Through Multi-Source Disease Similarity Networks
We improved computational drug repositioning by integrating multiple disease similarity networks into multiplex and multiplex-heterogeneous networks

![Network construction](https://github.com/hauldhut/MHDR/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: To store simulation results
  - **Evidence**: To store collected evidence
- **Figures**: To store generated figures from the simulation results

## How to run
- Install R packages
  - RandomWalkRestartMH, igraph, foreach, doParallel, ROCR, ggplot2, Metrics, hash
- Download the repo
- Follow instructions in the folder **Code** to run
- 
- *Note*: For large and complex networks (i.e., multiplex disease networks, and multiplex-hetergeneous networks of drugs and diseases), it is recommended to run on a multi-core and at least 16 GB RAM computer
