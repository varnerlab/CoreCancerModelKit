using CoreCancerModelKit

# setup calculation -
path_to_measurements_file = "./HL60_T.json"
organism_id = :LINE_HL_60_TB

# estimate the flux -
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file, organism_id; number_of_samples=10)
