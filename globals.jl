global models = ["ACCESS-CM2", "CMCC-ESM2", "GISS-E2-1-G", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0"]
global SSPs = ["ssp126", "ssp370"]
global time_periods = ["2041-2060", "2081-2100"]

global datadir = "." ### change to define a separate path to save data
global mech_dir = "$datadir/model_run_dec16/"
global cor_dir = "$datadir/correlative_model/outputs"

global quantiles = [0.025, 0.5, 0.975]