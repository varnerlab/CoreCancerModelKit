using CoreCancerModelKit

# setup calculation -
path_to_measurements_file = "./H72.json"
organism_id = :LINE_MDA_MB_231_ATCC

# estimate the flux -
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file, organism_id; number_of_samples=250)
