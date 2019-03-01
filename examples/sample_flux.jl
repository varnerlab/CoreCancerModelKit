using CoreCancerModelKit

# setup calculation -
path_to_measurements_file = "./H72.json"
organism_id = :LINE_MDA_MB_231_ATCC
number_of_biased_flux_samples = 10

# create a biased flux distribution -> maximize growth
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file, organism_id; number_of_samples=number_of_biased_flux_samples)

# create a default dd -
dd = generate_default_data_dictionary(organism_id)

# setup the sampling -
number_of_unbiased_flux_samples = 10
sample_archive = zeros(length(flux_distribution[:,1]))
sampling_bounds_array = []
for sample_index = 1:number_of_unbiased_flux_samples

    # draw a random number -
    r_value = rand()
    mean_flux_values = flux_distribution[:,1]
    std_error = flux_distribution[:,2]

    # create a sampling bounds array -
    sampling_bounds_array = results_array[sample_index].flux_bounds_array

    # sample -
    local_sample_archive = sample_flux_space(sampling_bounds_array,dd; number_of_samples = number_of_unbiased_flux_samples)

    # cache -
    sample_archive = [sample_archive local_sample_archive]
end

# cut off the leading zero -
sample_flux_ensemble = sample_archive[:,2:end]
