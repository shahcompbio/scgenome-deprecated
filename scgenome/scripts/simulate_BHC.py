from scgenome.simulation.simulation import do_simulations

OUTPUT_FP = "/Users/massoudmaher/data/tpois_ran_walk.json"
#OUTPUT_FP = "/work/shah/maherm/100t_pois_ran_walk.json"
PARAMS = {
    "samples_per_cluster": [8],
    "num_bin": [100, 500],
    "max_cn": [4],
    "alpha": [0.01, 0.05, 0.3, 0.6, 0.9],
    "distribution": ["poisson", "gaussian"],
    "init_lambdas": [(1, 3)],
    "jump_lambdas": [(1, 0.1)],
}
METRIC = "cityblock"

print("Starting simulations")
sims = do_simulations(2, PARAMS, naive_metric=METRIC, num_cores=None)
sims.to_json(OUTPUT_FP)
print(sims.head())


