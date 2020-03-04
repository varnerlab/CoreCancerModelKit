function convex_flux_estimation(path_to_measurements_file::String, organism_id::Symbol; number_of_samples::Int64 = 100, crowding_parameter::Float64 = 0.0044)
    
    # initalize -
    results_array = Array{VLOptimalConvexFluxResult,1}()

    # Declare a progress meter for user feedback -
    p_meter = Progress(number_of_samples,color=:yellow)

    # what is my flux budget?
    flux_budget = (1/crowding_parameter)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary(organism_id);

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file; crowding_parameter = crowding_parameter);
        updated_data_dictionary = constrain_measured_metabolites(updated_data_dictionary, path_to_measurements_file);

        # in this calculation, we need to also constrain the growth rate -
        updated_data_dictionary = constrain_specific_growth_rate(updated_data_dictionary, path_to_measurements_file)

        # ok, so now we need to seyup the Convex.jl problem -
        # Build the STM -
        S = updated_data_dictionary["stoichiometric_matrix"]
        (number_of_constraints,number_of_fluxes) = size(S);
        
        # upper and lower flux bounds -
        fba = updated_data_dictionary["flux_bounds_array"]
        LB = fba[:,1]
        UB = fba[:,2]

        # How many non-zero lower or upper species bounds do we have? (we'll need to turn these into slack variables) -
        b = zeros(number_of_constraints)
        
        # variabes and constraints -
        v = Variable(number_of_fluxes)
        p = minimize(sumsquares(S*v-b))
        p.constraints += [sum(abs(v))<=flux_budget]
        p.constraints += [v>=LB]
        p.constraints += [v<=UB]

        # solve -
        Convex.solve!(p, SCSSolver(verbose=true, max_iters=50000,eps=1e-3), warmstart=true)
        
        # get status flags, keep good solutions -
        status_flag = p.status
        calculated_flux_array = v.value 
        objective_value = p.optval

        # do we have an optimal soln?
        if (status_flag == :Optimal)
            
            # build a return type -
            fluxResult = VLOptimalConvexFluxResult()
            fluxResult.objective_value = objective_value;
            fluxResult.flux_array = calculated_flux_array;
            fluxResult.status_flag = status_flag;

            # get some problem setup information -
            fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

            # cache -
            push!(results_array, fluxResult);

            # user message -
            msg = "Success: optimal solution found. Completed $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p_meter; showvalues = [(:status,msg)]);
        else
            # user message -
            msg = "Failed: optimal was not solution found. Completed $(sample_index) of $(number_of_samples) trials ...";

            # update the progress bar -
            ProgressMeter.next!(p_meter; showvalues = [(:status,msg)]);
        end
    end

    # compute the ensemble -
    @info "Completed ...\r";

     # check - did we fail *all* the trials?
     if (isempty(results_array) == true)
        throw(ErrorException("Ooops! All trials failed. Please check your problem setup ..."))
    end

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