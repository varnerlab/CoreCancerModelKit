using CoreCancerModelKit

# setup calculation -
path_to_measurements_file = "./H72.json"
organism_id = :LINE_MDA_MB_231_ATCC
ko_index = 290 # PGI ko

# estimate the flux -
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file, organism_id; number_of_samples=10)

# grab a soln sample -
wild_type_soln = results_array[5].flux_array
wild_type_soln_fba = results_array[5].flux_bounds_array

# run moma -
moma_soln = moma_calculation(organism_id, wild_type_soln, wild_type_soln_fba, ko_index)
