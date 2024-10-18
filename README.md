# Interactive Graph Search Made Simple
## Dataset
- real datasets: data/amazon.txt, data/imagenet.txt, data/wordnet.txt
- synthetic datasets: data/syn/...
## Expriments:
### script
./run.sh will compile and run the expriments.
### result files
The experiment results will be saved in the  folder named result.
#### IGS Interaction:
| experiment    | result file |
| -------- | ------- |
| (#queries & #clicks) vs k |  amazon_poms_k; imagenet_poms_k; wordnet_poms_k; syn_poms_vary_k (synthetic dataset). |
| (#queries & #clicks) vs level | amazon_poms_lvl; imagenet_poms_lvl; wordnet_poms_lvl.|
| (#queries & #clicks) vs n |  syn_poms_vary_n   |
| (#queries & #clicks) vs r |  syn_poms_vary_r   |
| (#queries & #clicks) vs d |  syn_poms_vary_d   |

#### HPDFS-tree Computation Time:
| experiment    | result file |
| -------- | ------- |
| on real datasets: Amazon, Imagenet, Wordnet | real_hpdfs_time |
| synthetic data, time vs d|syn_hpdfs_n100w_r01_varyd
| synthetic data, time vs r |  syn_hpdfs_n100w_d30_varyr   |
| synthetic data, time vs n |  syn_hpdfs_d30_r01_varyn   |

#### IGS Computation Time:
| experiment    | result file |
| -------- | ------- |
| on real datasets: Amazon, Imagenet, Wordnet | real_poms_time |
| synthetic data, time vs d|syn_poms_vary_d_time
| synthetic data, time vs r |  syn_poms_vary_r_time   |
| synthetic data, time vs n |  syn_poms_vary_n_time   |