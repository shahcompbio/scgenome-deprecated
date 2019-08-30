from scgenome.simulation.simulation import many_poisson_bicluster

#OUTPUT_FP = "/Users/massoudmaher/data/pois_ran_walk.json"
OUTPUT_FP = "/work/shah/maherm/100t_pois_ran_walk.json"
TRIALS_PER_SET = 100
SAMPLES_PER_CLUSTER = [8]
NUM_BIN = [100, 500]
MAX_CN = [4]
ALPHA = [0.01, 0.05, 0.3, 0.6, 0.9]
INIT_LAMBDAS = [(3, 1)]
JUMP_LAMBDAS = [(1, 0.1)]

print("Starting simulations")
sims = many_poisson_bicluster(TRIALS_PER_SET, SAMPLES_PER_CLUSTER, NUM_BIN,
                              MAX_CN, ALPHA, INIT_LAMBDAS, JUMP_LAMBDAS,
                              num_cores=40)
sims.to_json(OUTPUT_FP)


