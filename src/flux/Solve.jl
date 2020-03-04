"""
TODO: Fill me in with some stuff ...
"""
function maximize_flux_at_index_array(flux_index_array::Array{Int64,1}, path_to_measurements_file::String, organism_id::Symbol; number_of_samples::Int64 = 100)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary(organism_id);

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_flux_at_index(flux_index_array, default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);
        updated_data_dictionary = constrain_measured_metabolites(updated_data_dictionary, path_to_measurements_file);

        # estimate the optimal flux distrubution -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # check the solution quality -
        if (exit_flag == 0 && status_flag == 5)

            # build a return type -
            fluxResult = VLOptimalFluxResult();
            fluxResult.objective_value = objective_value;
            fluxResult.flux_array = calculated_flux_array;
            fluxResult.dual_array = dual_value_array;
            fluxResult.uptake_array = uptake_array;
            fluxResult.exit_flag = exit_flag;
            fluxResult.status_flag = status_flag;

            # get some problem setup information -
            fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

            # cache -
            push!(results_array, fluxResult);

            # user message -
            msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p; showvalues = [(:status,msg)]);
        end
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end

"""
TODO: Fill me in with some stuff ...
"""
function maximize_specific_ha_production(path_to_measurements_file::String, organism_id::Symbol; number_of_samples::Int64 = 100)

    # error checking ...
    # TODO: check, are the path_to_measurements_file and path_to_growth_rate_file legit?

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # HA production rate has the tag -
    target_reaction_tag = "hyaluronan_synthase"

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary(organism_id);

        # pass the default dictionary to a customization method, in this case looking to max the target_reaction_tag -
        updated_data_dictionary = optimize_flux_at_index(default_data_dictionary, target_reaction_tag);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);
        updated_data_dictionary = constrain_measured_metabolites(updated_data_dictionary, path_to_measurements_file);

        # in this calculation, we need to also constrain the growth rate -
        updated_data_dictionary = constrain_specific_growth_rate(updated_data_dictionary, path_to_measurements_file)

        # estimate the optimal flux distrubution -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # check the solution quality -
        if (exit_flag == 0 && status_flag == 5)

            # build a return type -
            fluxResult = VLOptimalFluxResult();
            fluxResult.objective_value = objective_value;
            fluxResult.flux_array = calculated_flux_array;
            fluxResult.dual_array = dual_value_array;
            fluxResult.uptake_array = uptake_array;
            fluxResult.exit_flag = exit_flag;
            fluxResult.status_flag = status_flag;

            # get some problem setup information -
            fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

            # cache -
            push!(results_array, fluxResult);

            # user message -
            msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p; showvalues = [(:status,msg)]);
        end
    end

    # check - did we fail *all* the trials?
    if (isempty(results_array) == true)
        throw(ErrorException("all trials failed. Please check your problem setup"))
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end


"""
TODO: Fill me in with some stuff ...
"""
function maximize_specific_growth_rate(path_to_measurements_file::String, organism_id::Symbol; number_of_samples::Int64 = 100, crowding_parameter::Float64 = 0.0044)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary(organism_id);

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file; crowding_parameter = crowding_parameter);
        updated_data_dictionary = constrain_measured_metabolites(updated_data_dictionary, path_to_measurements_file);

        # estimate the optimal flux distrubution -
        (objective_value_tmp, calculated_flux_array_tmp, dual_value_array_tmp, uptake_array_tmp, exit_flag_tmp, status_flag_tmp) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # re-run w/-1's in the crowding cosntaint -
        idx_negative = findall(x->x<0.0,calculated_flux_array_tmp)
        add_constraint = updated_data_dictionary["additional_constraint_array"]
        add_constraint[end,idx_negative] .= -1.0
        updated_data_dictionary["additional_constraint_array"] = add_constraint
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # check the solution quality -
        if (exit_flag == 0 && status_flag == 5)

            # build a return type -
            fluxResult = VLOptimalFluxResult();
            fluxResult.objective_value = objective_value;
            fluxResult.flux_array = calculated_flux_array;
            fluxResult.dual_array = dual_value_array;
            fluxResult.uptake_array = uptake_array;
            fluxResult.exit_flag = exit_flag;
            fluxResult.status_flag = status_flag;

            # get some problem setup information -
            fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

            # cache -
            push!(results_array, fluxResult);

            # user message -
            msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p; showvalues = [(:status,msg)]);
        else
            # user message -
            msg = "Failed: $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p; showvalues = [(:status,msg)]);
        end
    end

    # check - did we fail *all* the trials?
    if (isempty(results_array) == true)
        throw(ErrorException("all trials failed. Please check your problem setup"))
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end

function sample_flux_space_with_experimental_constraints(solution_bounds_array::Array{Float64,2}, number_of_samples::Int64)

    # load the default data_dictionary -
    default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id);

    # update dictionary with experimental data?
    updated_data_dictionary = constrain_measured_fluxes(default_data_dictionary, path_to_measurements_file)

    # call the sample method -
    return sample_flux_space(solution_bounds_array, updated_data_dictionary,number_of_samples)
end
