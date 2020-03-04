"""
TODO: Fill me in with some stuff ...
"""
function optimize_specific_growth_rate(data_dictionary::Dict{String,Any})

    # get the list of reaction names -
    list_of_reaction_name_strings = data_dictionary["list_of_reaction_name_strings"]

    # ok, we need to find the index of the EX_ha(e) rate -
    index_target_reaction = find_index_of_reaction(list_of_reaction_name_strings,"growth")
    if (index_target_reaction == nothing)
        throw(ErrorException("missing reaction tag $(target_reaction_tag)"))
    end

    # update the c-vector -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]
    objective_coefficient_array[index_target_reaction] = -1.0
    data_dictionary["objective_coefficient_array"] = objective_coefficient_array

    # return -
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function optimize_flux_at_index(data_dictionary::Dict{String,Any}, target_reaction_tag::String)

    # update the c-vector -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]

    # get the list of reaction names -
    list_of_reaction_name_strings = data_dictionary["list_of_reaction_name_strings"]

    # ok, we need to find the index of the EX_ha(e) rate -
    index_target_reaction = find_index_of_reaction(list_of_reaction_name_strings,target_reaction_tag)
    if (index_target_reaction == nothing)
        throw(ErrorException("missing reaction tag $(target_reaction_tag)"))
    end

    objective_coefficient_array[index_target_reaction] = -1.0   # cellmass is always at the end -
    data_dictionary["objective_coefficient_array"] = objective_coefficient_array

    # return -
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function optimize_flux_at_index(index_array::Array{Int64,1}, data_dictionary::Dict{String,Any})

    # update the c-vector -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]

    # iterate through the index array, and flip the flag on the coefficient array -
    for index_value in index_array
        objective_coefficient_array[index_value] = -1
    end

    # return -
    return data_dictionary
end

function bound_specific_growth_rate(data_dictionary::Dict{String,Any}, path_to_growth_rate_file::String)

    # TODO: check the growth rate file -
    growth_rate_dictionary = JSON.parsefile(path_to_growth_rate_file)

    # get the flux bounds -
    copy_flux_bounds_array = deepcopy(data_dictionary["flux_bounds_array"])

    # get the lower_bound, and upper_bound -
    lower_bound = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["lower_bound"])
    upper_bound = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["upper_bound"])
    mean_value = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["mean_value"])

    # growth rate is *always* the last value -
    copy_flux_bounds_array[end,1] = lower_bound
    copy_flux_bounds_array[end,2] = upper_bound

    # update -
    data_dictionary["flux_bounds_array"] = copy_flux_bounds_array

    # return -
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function constrain_specific_growth_rate(data_dictionary::Dict{String,Any}, path_to_growth_rate_file::String)

    # TODO: check the growth rate file -
    growth_rate_dictionary = JSON.parsefile(path_to_growth_rate_file)

    # get the flux bounds -
    copy_flux_bounds_array = deepcopy(data_dictionary["flux_bounds_array"])

    # get the lower_bound, and upper_bound -
    lower_bound = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["lower_bound"])
    upper_bound = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["upper_bound"])
    mean_value = parse(Float64,growth_rate_dictionary["specific_growth_rate_measurement"]["mean_value"])

    # compute a random factor -
    random_factor = rand()
    growth_rate = (1-random_factor)*lower_bound+random_factor*upper_bound

    # growth rate is *always* the last value -
    copy_flux_bounds_array[end,1] = growth_rate
    copy_flux_bounds_array[end,2] = growth_rate

    # update -
    data_dictionary["flux_bounds_array"] = copy_flux_bounds_array

    # return -
    return data_dictionary
end

function constrain_measured_fluxes(data_dictionary::Dict{String,Any}, path_to_measurements_file::String; crowding_parameter::Float64 = 0.0044)

    # TODO: is the path_to_measurements_file legit?
    measurements_dictionary = JSON.parsefile(path_to_measurements_file)

    # get some stuff from the dd -
    rxn_name_list = data_dictionary["list_of_reaction_name_strings"]
    flux_bounds_array = data_dictionary["flux_bounds_array"]
    (number_of_reactions,nc) = size(flux_bounds_array)

    # do we have any flux ratios?
    flux_ratio_dictionary_array = measurements_dictionary["flux_ratio_measurements"]

    # initialize -
    number_of_additional_constraints = length(flux_ratio_dictionary_array)
    constraint_counter = 1
    additional_constraint_array = zeros((number_of_additional_constraints + 1),number_of_reactions)
    for (reaction_key, local_measurement_dict) in flux_ratio_dictionary_array

        # ok, so we have a reaction key - find the index of this reaction in the reaction list -
        idx_reaction_match = find_index_of_reaction(rxn_name_list, reaction_key)

        # if we have this reaction, then update the bounds array -
        if (idx_reaction != nothing)

            # get ratio value -
            ratio_value = parse(Float64,local_measurement_dict["mean_value"])

            # get the actual index -
            idx_reaction = (getindex(idx_reaction_match))[1]

            # update the constraint -
            additional_constraint_array[constraint_counter,idx_reaction] = 1.0
            additional_constraint_array[constraint_counter,glucose_uptake_index] = -1*ratio_value

            # update -
            constraint_counter = constraint_counter + 1
        end
    end

    # crowding - always the last additional constraint -
    for reaction_index = 1:number_of_reactions
        additional_constraint_array[end,reaction_index] = 1.0
    end

    sba = data_dictionary["species_bounds_array"]
    sba[end,2] = (1.0/crowding_parameter);
    data_dictionary["species_bounds_array"] = sba

    # cache the additional constraints -
    data_dictionary["additional_constraint_array"] = additional_constraint_array

    # iterate through the measurements_dictionary and set values -
    individual_measuerment_dictionaries = measurements_dictionary["exchange_flux_measurements"]
    for (reaction_key, local_measurement_dict) in individual_measuerment_dictionaries

        # ok, so we have a reaction key - find the index of this reaction in the reaction list -
        idx_reaction_match = find_index_of_reaction(rxn_name_list, reaction_key)

        # if we have this reaction, then update the bounds array -
        if (idx_reaction_match != nothing)

            # get the measured value -
            mean_measured_value = parse(Float64,local_measurement_dict["mean_value"])

            # get the CV value -
            coefficient_of_variation = parse(Float64, local_measurement_dict["coefficient_of_variation"])

            # what is the measured value?
            distrubution = Distributions.Normal(mean_measured_value, coefficient_of_variation*mean_measured_value)
            measured_value = abs(rand(distrubution))

            # get the actual index -
            idx_reaction = idx_reaction_match

            # what is the directionality of this measurement?
            directionality = Symbol(local_measurement_dict["directionality"])
            if directionality == :input

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = 0.0

            elseif directionality == :output

                flux_bounds_array[idx_reaction,1] = 0
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :free

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :input_fixed

                flux_bounds_array[idx_reaction,1] = -1*measured_value
                flux_bounds_array[idx_reaction,2] = -1*measured_value

            elseif directionality == :output_fixed

                flux_bounds_array[idx_reaction,1] = measured_value
                flux_bounds_array[idx_reaction,2] = measured_value

            elseif directionality == :nil

                flux_bounds_array[idx_reaction,1] = -1e-6
                flux_bounds_array[idx_reaction,2] = 1e-6

            else
                # TODO: issue a warning ...
                println("Missing $(idx_reaction_match)")
            end
        else
            #println("missing $(reaction_key)??")
        end
    end

    # update -
    data_dictionary["flux_bounds_array"] = flux_bounds_array

    (number_of_reactions, nc) = size(flux_bounds_array)
    for reaction_index = 1:number_of_reactions

        msg = "$(reaction_index),$(flux_bounds_array[reaction_index,1]),$(flux_bounds_array[reaction_index,2])"
        # println(msg)

    end

    # return modified dictionary -
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function generate_exchange_flux_index_array(reaction_name_array::Array{String,1},pattern::String)

    # initialize -
    exchange_flux_index_array = Int64[]

    # we are looking for a name starting w/EX_
    number_of_reactions = length(reaction_name_array)
    for reaction_index = 1:number_of_reactions

        # grab the name -
        reaction_name = reaction_name_array[reaction_index]

        # check -
        if (startswith(reaction_name,pattern) == true)
            push!(exchange_flux_index_array,reaction_index)
        end
    end

    # return -
    return exchange_flux_index_array
end

function generate_path_to_cellmass_file(organism_id)

    # load the map -
    path_cellmass_map_file = "$(path_to_package)/cobra/config/data/cellmass/Map.json"

    # map -
    cellmass_file_dictionary = JSON.parsefile(path_cellmass_map_file)

    # what is the organism key?
    organism_id_key = String(organism_id)

    # check - do we have this key?
    if (haskey(cellmass_file_dictionary, organism_id_key) == true)

        # what is the file name?
        file_name = cellmass_file_dictionary[organism_id_key]["cellmass_file_name"]

        # build path -
        path_to_cellmass_file = "$(path_to_package)/cobra/config/data/cellmass/$(file_name)"

        # return -
        return path_to_cellmass_file
    end

    return "$(path_to_package)/cobra/config/data/cellmass/Default-Biomass.csv"
end

function constrain_measured_metabolites(data_dictionary::Dict{String,Any}, path_to_measurements_file::String)

    # from the data dictionary, get the list of metabolite symbols -
    list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

    # load the measurements file -
    measurements_dictionary = JSON.parsefile(path_to_measurements_file)

    # get the stoichiometric_matrix -
    stoichiometric_matrix = copy(data_dictionary["stoichiometric_matrix"])

    # get the metabolite measurement information -
    absolute_metabolite_measurement_dictionary = measurements_dictionary["absolute_metabolite_measurements"]
    for (metabolite_key, local_measurement_dict) in absolute_metabolite_measurement_dictionary

        # ok, so we have a metabolite_key - find the index of this metabolite in the metabolite list -
        idx_metabolite_match = findall(list_of_metabolite_symbols .== metabolite_key)

        # if we have this metabolite, update the STM -
        if (isempty(idx_metabolite_match) == false)

            # get experimental values -
            mean_ratio_value = parse(Float64,local_measurement_dict["mean_value"])
            lower_bound = parse(Float64,local_measurement_dict["lower_bound"])
            upper_bound = parse(Float64,local_measurement_dict["upper_bound"])

            # compute a value -
            r_value = rand()
            value = (1 - r_value)*lower_bound + r_value*(upper_bound)

            # update the stm -> the growth rate is *always at the end, last reaction*
            metabolite_index = (getindex(idx_metabolite_match))[1]
            old_value = stoichiometric_matrix[metabolite_index,end]
            stoichiometric_matrix[metabolite_index,end] = (old_value - value)
        end
    end

    data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
    return data_dictionary
end

function generate_test_data_dictionary(organism_id::Symbol)

    # Hardcode my path information -
    #path_to_cobra_mat_file = "$(path_to_package)/cobra/config/matlab_cobra_files/modelCore.mat"
    #model_name = "modelCore"

    # path_to_cobra_mat_file = "$(path_to_package)/cobra/config/matlab_cobra_files/CoreCancerModel_v1.mat"
    # model_name = "CoreCancerModel_v1"

    path_to_cobra_mat_file = "$(path_to_package)/cobra/config/matlab_cobra_files/Test_v1.mat"
    model_name = "Test_v1"

    # load the biophysical_constants dictionary -
    default_biophysical_dictionary = load_default_biophysical_dictionary(organism_id)

    # TODO: check if string is a legit path to the cobra file, and the model name is ok
    # load the *.mat code from the cobra code folder -
    file = matopen(path_to_cobra_mat_file)
    cobra_dictionary = read(file, model_name)
    close(file)

    # get some stuff that we may need later ...
    # get the species symbol list -
    list_of_metabolite_symbols = cobra_dictionary["mets"]

    # Setup: the stoichiometric matrix -
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # we need to add the biomass equation for cellmass -
    # which cellmass defn are we using
    path_to_cellmass_file = generate_path_to_cellmass_file(organism_id)
    biomass_eqn_data = readdlm(path_to_cellmass_file,',')
    biomass_eqn = zeros(number_of_species)

    # how many biomass components do we have?
    number_of_biomass_species = length(biomass_eqn_data[:,1])
    for biomass_component_index = 1:number_of_biomass_species

        # get the component symbol and stcoeff -
        component_symbol = biomass_eqn_data[biomass_component_index,1]
        component_stoichiometric_coeff = biomass_eqn_data[biomass_component_index,2]

        # ok, so now we need to find the index of this component -
        index_of_species = find_index_of_species(cobra_dictionary["mets"],component_symbol)

        # check - do we have a value for the species index?
        if (index_of_species == nothing)
            throw(ArgumentException("missing $(component_symbol) when constructing the biomass equation"))
        end

        # update the entry in the biomass_eqn -
        biomass_eqn[index_of_species] = component_stoichiometric_coeff
    end

    # add the new *column* to the stoichometric array -
    stoichiometric_matrix = [stoichiometric_matrix biomass_eqn]

    # what is the *final* size of the system?
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # objective coefficient array -
    objective_coefficient_array = zeros(number_of_reactions)
    objective_coefficient_array_palsson = cobra_dictionary["c"]
    for (index,objective_value) in enumerate(objective_coefficient_array_palsson)
        objective_coefficient_array[index] = objective_value
    end

    # Add growh to rev list (0 = not reversible, 1 = reversible)
    # reversible_reactions = cobra_dictionary["rev"]
    # reversible_reactions = [reversible_reactions ; 0.0]
    # (number_of_fluxes, nc) = size(reversible_reactions)
    # reversible_reaction_array = zeros(number_of_fluxes)
    # for index = 1:number_of_fluxes
    #     reversible_reaction_array[index] = reversible_reactions[index][1]
    # end

    # # calculate the reaction name -
    reaction_name_array_tmp = cobra_dictionary["rxns"]
    reaction_name_array = String[]
    for rxn_name in reaction_name_array_tmp
        push!(reaction_name_array,rxn_name)
    end

    # add growth to reaction name list -
    push!(reaction_name_array,"growth")

    # get list of gene symbols -
    list_of_gene_symbols = cobra_dictionary["genes"]

    # default flux bounds array -
    default_flux_bounds_array = zeros(number_of_reactions,2)
    lb = cobra_dictionary["lb"] # lower bound -
    ub = cobra_dictionary["ub"] # upper bound -
    for reaction_index = 1:number_of_reactions - 1
        default_flux_bounds_array[reaction_index,1] = lb[reaction_index]
        default_flux_bounds_array[reaction_index,2] = ub[reaction_index]
    end

    # add default growth rate constraint?
    default_flux_bounds_array[end,1] = 0.0
    default_flux_bounds_array[end,2] = 1.0

    # build reversible reaction flag array -
    reversible_reaction_flag_array = ones(number_of_reactions)
    for reaction_index = 1:number_of_reactions

        if (default_flux_bounds_array[reaction_index,1] >= 0.0)
            reversible_reaction_flag_array[reaction_index] = 0.0
        end
    end
    exchange_flux_index_array = generate_exchange_flux_index_array(reaction_name_array,"EX_")

    # all EX_ are reversible -
    for exchange_index in exchange_flux_index_array
        reversible_reaction_flag_array[exchange_index] = 1.0
    end

    # calculate the default vamx (based upon biophysical_constants/literature) -
    (default_vmax, default_enzyme_concentration) = calculate_default_vmax(default_biophysical_dictionary)

    # load the ec_number file -
    path_to_ec_file::String = "$(path_to_package)/cobra/config/data/ec_numbers.dat"
    ec_number_dictionary =  load_ec_mapping_file(path_to_ec_file)

    # load the gene order array -
    path_to_gene_file = "$(path_to_package)/cobra/config/data/gene_list.dat"
    gene_order_array = load_gene_order_file(path_to_gene_file)

    # calculate the vamx array -
    path_to_brenda_data = "$(path_to_package)/cobra/config/data/ECN-CoreCancerModel-Palsson-SciReports-2017.json"
    total_vmax_array = generate_vmax_array(path_to_brenda_data, default_biophysical_dictionary, gene_order_array, ec_number_dictionary)

    # correct the total_vmax_array for this model, using the rules -
    local_data_dictionary = Dict{String,Any}()
    local_data_dictionary["default_vmax_value"] = default_vmax
    #model_vmax_array = calculate_rules_vector(local_data_dictionary, total_vmax_array)

    # update the default bounds array w/our "default" biophysical_constants -
    # flux_bounds_array = update_default_flux_bounds_array(default_flux_bounds_array, model_vmax_array, reversible_reaction_flag_array)

    # What sense do we do? (by default we min)
    is_minimum_flag = true

    # construct the reaction string list -
    list_of_chemical_reaction_strings = reconstruct_reaction_string_list(cobra_dictionary)

    # setup default additional constaints array -
    number_of_additional_constraints = 1
    species_bounds_array = zeros((number_of_species + number_of_additional_constraints),2)
    species_bounds_array[end,1] = 0.0
    species_bounds_array[end,2] = (1.0/0.0040);

    # =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
    data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
    data_dictionary["flux_bounds_array"] = default_flux_bounds_array
    data_dictionary["species_bounds_array"] = species_bounds_array
    data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
    data_dictionary["list_of_gene_symbols"] = list_of_gene_symbols
    data_dictionary["list_of_reaction_name_strings"] = reaction_name_array
    data_dictionary["list_of_chemical_reaction_strings"] = list_of_chemical_reaction_strings
	data_dictionary["is_minimum_flag"] = is_minimum_flag
    data_dictionary["reversible_reaction_flag_array"] = reversible_reaction_flag_array
    #data_dictionary["model_vmax_array"] = model_vmax_array

    # stuff for rules -
    data_dictionary["default_vmax"] = default_vmax

    # in case we need something later -
    data_dictionary["cobra_dictionary"] = cobra_dictionary

    # <a> -
    data_dictionary["average_crowding_constraint_value"] = 0.0040;  # hardcode the PNAS value -
    # ========================================================================================= #
    return data_dictionary
end

"""
TODO: Fill me in with some stuff ...
"""
function generate_default_data_dictionary(organism_id::Symbol)

    # Hardcode my path information -
    path_to_cobra_mat_file = "$(path_to_package)/cobra/config/matlab_cobra_files/CoreCancerModel_v1.mat"
    model_name = "CoreCancerModel_v1"

    # load the biophysical_constants dictionary -
    default_biophysical_dictionary = load_default_biophysical_dictionary(organism_id)

    # TODO: check if string is a legit path to the cobra file, and the model name is ok
    # load the *.mat code from the cobra code folder -
    file = matopen(path_to_cobra_mat_file)
    cobra_dictionary = read(file, model_name)
    close(file)

    # get some stuff that we may need later ...
    # get the species symbol list -
    list_of_metabolite_symbols = cobra_dictionary["mets"]

    # Setup: the stoichiometric matrix -
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # we need to add the biomass equation for cellmass -
    # which cellmass defn are we using
    path_to_cellmass_file = generate_path_to_cellmass_file(organism_id)
    biomass_eqn_data = readdlm(path_to_cellmass_file,',')
    biomass_eqn = zeros(number_of_species)

    # how many biomass components do we have?
    number_of_biomass_species = length(biomass_eqn_data[:,1])
    for biomass_component_index = 1:number_of_biomass_species

        # get the component symbol and stcoeff -
        component_symbol = biomass_eqn_data[biomass_component_index,1]
        component_stoichiometric_coeff = biomass_eqn_data[biomass_component_index,2]

        # ok, so now we need to find the index of this component -
        index_of_species = find_index_of_species(cobra_dictionary["mets"],component_symbol)

        # check - do we have a value for the species index?
        if (index_of_species == nothing)
            throw(ArgumentException("missing $(component_symbol) when constructing the biomass equation"))
        end

        # update the entry in the biomass_eqn -
        biomass_eqn[index_of_species] = component_stoichiometric_coeff
    end

    # add the new *column* to the stoichometric array -
    stoichiometric_matrix = [stoichiometric_matrix biomass_eqn]

    # what is the *final* size of the system?
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)

    # objective coefficient array -
    objective_coefficient_array = zeros(number_of_reactions)
    objective_coefficient_array_palsson = cobra_dictionary["c"]
    for (index,objective_value) in enumerate(objective_coefficient_array_palsson)
        objective_coefficient_array[index] = objective_value
    end

    # Add growh to rev list (0 = not reversible, 1 = reversible)
    # reversible_reactions = cobra_dictionary["rev"]
    # reversible_reactions = [reversible_reactions ; 0.0]
    # (number_of_fluxes, nc) = size(reversible_reactions)
    # reversible_reaction_array = zeros(number_of_fluxes)
    # for index = 1:number_of_fluxes
    #     reversible_reaction_array[index] = reversible_reactions[index][1]
    # end

    # # calculate the reaction name -
    reaction_name_array_tmp = cobra_dictionary["rxns"]
    reaction_name_array = String[]
    for rxn_name in reaction_name_array_tmp
        push!(reaction_name_array,rxn_name)
    end

    # add growth to reaction name list -
    push!(reaction_name_array,"growth")

    # get list of gene symbols -
    list_of_gene_symbols = cobra_dictionary["genes"]

    # default flux bounds array -
    default_flux_bounds_array = zeros(number_of_reactions,2)
    lb = cobra_dictionary["lb"] # lower bound -
    ub = cobra_dictionary["ub"] # upper bound -
    for reaction_index = 1:number_of_reactions - 1
        default_flux_bounds_array[reaction_index,1] = lb[reaction_index]
        default_flux_bounds_array[reaction_index,2] = ub[reaction_index]
    end

    # add default growth rate constraint?
    default_flux_bounds_array[end,1] = 0.0
    default_flux_bounds_array[end,2] = 1.0

    # build reversible reaction flag array -
    reversible_reaction_flag_array = ones(number_of_reactions)
    for reaction_index = 1:number_of_reactions

        if (default_flux_bounds_array[reaction_index,1] >= 0.0)
            reversible_reaction_flag_array[reaction_index] = 0.0
        end
    end
    exchange_flux_index_array = generate_exchange_flux_index_array(reaction_name_array,"EX_")

    # all EX_ are reversible -
    for exchange_index in exchange_flux_index_array
        reversible_reaction_flag_array[exchange_index] = 1.0
    end

    # calculate the default vamx (based upon biophysical_constants/literature) -
    (default_vmax, default_enzyme_concentration) = calculate_default_vmax(default_biophysical_dictionary)

    # load the ec_number file -
    path_to_ec_file::String = "$(path_to_package)/cobra/config/data/ec_numbers.dat"
    ec_number_dictionary =  load_ec_mapping_file(path_to_ec_file)

    # load the gene order array -
    path_to_gene_file = "$(path_to_package)/cobra/config/data/gene_list.dat"
    gene_order_array = load_gene_order_file(path_to_gene_file)

    # calculate the vamx array -
    path_to_brenda_data = "$(path_to_package)/cobra/config/data/ECN-CoreCancerModel-Palsson-SciReports-2017.json"
    total_vmax_array = generate_vmax_array(path_to_brenda_data, default_biophysical_dictionary, gene_order_array, ec_number_dictionary)

    # correct the total_vmax_array for this model, using the rules -
    local_data_dictionary = Dict{String,Any}()
    local_data_dictionary["default_vmax_value"] = default_vmax
    #model_vmax_array = calculate_rules_vector(local_data_dictionary, total_vmax_array)

    # update the default bounds array w/our "default" biophysical_constants -
    #flux_bounds_array = update_default_flux_bounds_array(default_flux_bounds_array, model_vmax_array, reversible_reaction_flag_array)
    (number_of_bounds, number_of_cols) = size(default_flux_bounds_array)
    flux_bounds_array = zeros(number_of_bounds,2)
    for bound_index = 1:number_of_bounds
        
        lower_bound = default_flux_bounds_array[bound_index,1]
        upper_bound = default_flux_bounds_array[bound_index,2]

        if (lower_bound!=0.0)
            flux_bounds_array[bound_index,1] = sign(lower_bound)*default_vmax
        end

        if (upper_bound!=0.0)
            flux_bounds_array[bound_index,2] = sign(upper_bound)*default_vmax
        end
    end


    # What sense do we do? (by default we min)
    is_minimum_flag = true

    # construct the reaction string list -
    list_of_chemical_reaction_strings = reconstruct_reaction_string_list(cobra_dictionary)

    # setup default additional constaints array -
    number_of_additional_constraints = 1
    species_bounds_array = zeros((number_of_species + number_of_additional_constraints),2)
    species_bounds_array[end,1] = 0.0
    species_bounds_array[end,2] = (1.0/0.0040);

    # =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
    data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
    data_dictionary["flux_bounds_array"] = flux_bounds_array
    data_dictionary["species_bounds_array"] = species_bounds_array
    data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
    data_dictionary["list_of_gene_symbols"] = list_of_gene_symbols
    data_dictionary["list_of_reaction_name_strings"] = reaction_name_array
    data_dictionary["list_of_chemical_reaction_strings"] = list_of_chemical_reaction_strings
	data_dictionary["is_minimum_flag"] = is_minimum_flag
    data_dictionary["reversible_reaction_flag_array"] = reversible_reaction_flag_array
    #data_dictionary["model_vmax_array"] = model_vmax_array

    # stuff for rules -
    data_dictionary["default_vmax"] = default_vmax

    # in case we need something later -
    data_dictionary["cobra_dictionary"] = cobra_dictionary
    # ========================================================================================= #
    return data_dictionary
end

# - PRIVATE HELPER METHODS -------------------------------------------------------------------- #
function load_gene_order_file(path_to_gene_file)

    # gene array -
    gene_symbol_array = String[]

    # open the gene file -
    # TODO: We need to check to see if this is a legit path -
    open(path_to_gene_file) do f

        # load the file into the buffer -
        buffer = read(f, String)

        # split along new line -
        list_of_records = split(buffer,'\n')

        # how many records do we have?
        number_of_records = length(list_of_records)
        for record_index = 1:number_of_records

            # convert the record to a string -
            record = string(list_of_records[record_index])

            # is the record empty?
            if (isempty(record) == false)

                # split into fields -
                field_array = split(record,',')

                # grab the ec number -
                order_index = field_array[1]
                gene_symbol = field_array[2]

                # add a record -
                push!(gene_symbol_array, gene_symbol)
            end
        end
    end

    # return -
    return gene_symbol_array
end

function load_ec_mapping_file(path_to_ec_file::String)

    # initialize the ec number array -
    ec_number_dictionary = Dict{String,String}()

    # open the ec file -
    open(path_to_ec_file) do f

        # load the file into the buffer -
        buffer = read(f, String)

        # split along new line -
        list_of_records = split(buffer,'\n')

        # how many records do we have?
        number_of_records = length(list_of_records)
        for record_index = 1:number_of_records

            # convert the record to a string -
            record = string(list_of_records[record_index])

            # is the record empty?
            if (isempty(record) == false)

                # split into fields -
                field_array = split(record,'\t')

                # grab the ec number -
                bg_number = field_array[1]
                ec_number = field_array[2]

                # add a record -
                ec_number_dictionary[bg_number] = ec_number

            end
        end
    end

    # return -
    return ec_number_dictionary
end

function generate_vmax_array(path_to_brenda_data::String, biophysical_dictionary::Dict{String,Any}, gene_symbol_array::Array{String,1}, ec_number_dictionary::Dict{String,String})

    # TODO: is the path_to_brenda_data legit?
    kcat_measurements_dictionary = JSON.parsefile(path_to_brenda_data)

    # initialize -
    vmax_array = Float64[]

    # calculate a default value -
    (default_vmax, default_enzyme_concentration) = calculate_default_vmax(biophysical_dictionary)

    # for each gene symbol, lookup an ec number, then get a range of measured kcats -
    for gene_symbol in gene_symbol_array

        # do we have this gene key?
        if (haskey(ec_number_dictionary, gene_symbol) == true)

            # what is the ec number for this key?
            ec_number_key = ec_number_dictionary[gene_symbol]

            # ok, so we have an ec number. Do we have a kcat measurement for this ec number?
            if (haskey(kcat_measurements_dictionary, ec_number_key) == true)


                # ok, we have an ec_number w/data - get the data
                # this will come back as an array of dictionaries
                local_vmax_array = Float64[]
                array_data_dictionaries = kcat_measurements_dictionary[ec_number_key]
                for local_data_dictionary::Dict{String,Any} in array_data_dictionaries

                    # get the turnoverNumber -
                    turnoverNumber = parse(Float64,local_data_dictionary["turnoverNumber"])*(3600) # hr^-1

                    # there is a "strange" data issue w/the Brenda data. They use -999 to mean no value
                    if (turnoverNumber>0)
                        # cache -
                        push!(local_vmax_array, turnoverNumber*default_enzyme_concentration)
                    end
                end

                if isempty(local_vmax_array) == true
                    push!(local_vmax_array, default_vmax)
                end

                # what is the min/max -
                min_value = minimum(local_vmax_array)
                max_value = maximum(local_vmax_array)

                # calculate a value (uniform random between the measurements)
                sampled_vmax_value = min_value + (max_value - min_value)*rand()

                # yes - we have an experimental value -
                push!(vmax_array, sampled_vmax_value)
            else

                # nope - bummer. Use the default value -
                push!(vmax_array, default_vmax)
            end
        else

            # nope - bummer. Use the default value -
            push!(vmax_array, default_vmax)
        end
    end

    # return -
    return vmax_array
end

function calculate_default_vmax(biophysical_dictionary::Dict{String,Any})

    # ok, so we need to get some stuff from the dictionary -
    # TODO: Check these keys are contained in the dictionary
    default_turnover_number = parse(Float64,biophysical_dictionary["biophysical_constants"]["default_turnover_number"]["value"])              # convert to h^-1
    default_enzyme_concentration = parse(Float64,biophysical_dictionary["biophysical_constants"]["default_enzyme_concentation"]["value"])     # mumol/gDW

    # calculate the default VMax -
    default_vmax = (default_turnover_number)*(default_enzyme_concentration)*(3600)

    return (default_vmax, default_enzyme_concentration)
end

function update_default_flux_bounds_array(flux_bounds_array::Array{Float64,2}, vmax_array::Array{Float64,1}, reversible_flag_array::Array{Float64,1})

    # how many bounds do we have?
    (number_of_bounds, number_of_columns) = size(flux_bounds_array)

    # initialize and update the array -
    updated_flux_bounds_array = zeros(number_of_bounds,number_of_columns)
    for flux_index = 1:number_of_bounds

        # get vmax value -
        vmax_value = vmax_array[flux_index]

        # setup new bounds -
        updated_flux_bounds_array[flux_index,1] = 0.0
        updated_flux_bounds_array[flux_index,2] = vmax_value
        if (reversible_flag_array[flux_index] == 1.0)
            updated_flux_bounds_array[flux_index,1] = -1*vmax_value
        end
    end

    # return -
    return updated_flux_bounds_array
end

function load_default_biophysical_dictionary(organism_id::Symbol)

    # load the organism_id mapping file -
    top_level_path = pwd()
    path_to_mapping_file = "$(path_to_package)/flux/config/Mapping.json"

    # TODO: Check do we have a mapping file?
    organism_id_map = JSON.parsefile(path_to_mapping_file)

    # TODO: Throw an error if we are missing this organism
    local_path_to_biophysical_dictionary = organism_id_map[String(organism_id)]["path_to_default_constants"]

    # TODO: check - does this path point to a file?
    global_path_to_biophysical_dictionary = "$(path_to_package)"*"$(local_path_to_biophysical_dictionary)"
    default_biophysical_dictionary = JSON.parsefile(global_path_to_biophysical_dictionary)

    # return -
    return default_biophysical_dictionary
end

function reconstruct_reaction_string_list(cobra_dictionary, path_to_reaction_file)

    # initialize -
    reaction_string_buffer = String[]

    # load the file -
    # file = matopen(path_to_model_file)
    # cobra_dictionary = read(file,data_name)
    # close(file)

    # get the stoichiometric array -
    stoichiometric_matrix = Matrix(cobra_dictionary["S"])
    list_of_reaction_tags = cobra_dictionary["rxns"]
    list_of_species = cobra_dictionary["mets"]
    list_of_reversible_reactions = cobra_dictionary["rev"]

    # what is the size?
    (number_of_species,number_of_reactions) = size(stoichiometric_matrix)
    for reaction_index = 1:number_of_reactions

        # initialize empty buffer -
        buffer = ""

        # get the reaction tag -
        reaction_tag_string = list_of_reaction_tags[reaction_index]

        # add the tag to the buffer -
        buffer *= "$(reaction_index),$(reaction_tag_string),"

        # find the reactants -
        idx_reactants = findall(stoichiometric_matrix[:,reaction_index].<0.0)
        if (isempty(idx_reactants) == true)
            buffer *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_reactants)
            counter = 1
            for index in idx_reactants

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != -1.0)
                    # add data to the buffer -
                    buffer *= "$(abs(stcoeff))*$(metabolite_string)"
                else
                    # add data to the buffer -
                    buffer *= "$(metabolite_string)"
                end



                # do we have more?
                if (counter < number_of_species)
                    buffer *= "+"
                else
                    buffer *= ","
                end

                counter = counter + 1
            end
        end

        # find the products -
        idx_products = findall(stoichiometric_matrix[:,reaction_index].>0.0)
        if (isempty(idx_products) == true)
            buffer *= "[],"
        else

            # how many species do we have?
            number_of_species = length(idx_products)
            counter = 1
            for index in idx_products

                # get the metabolite -
                metabolite_string = list_of_species[index]
                stcoeff = stoichiometric_matrix[index,reaction_index]

                if (stcoeff != 1.0)
                    # add data to the buffer -
                    buffer *= "$(stcoeff)*$(metabolite_string)"
                else
                    # add data to the buffer -
                    buffer *= "$(metabolite_string)"
                end

                # do we have more?
                if (counter < number_of_species)
                    buffer *= "+"
                else
                    buffer *= ","
                end

                counter = counter + 1
            end
        end

        # is this reaction reversible?
        rev_value = list_of_reversible_reactions[reaction_index]
        if (rev_value == 1.0)
            buffer *= "-inf,inf"
        else
            buffer *= "0,inf"
        end

        # add buffer to list of strings -
        push!(reaction_string_buffer,buffer)
    end

    # write to disk -
    open(path_to_reaction_file,"w") do f

        for statement in reaction_string_buffer
            write(f,"$(statement)\n")
        end
    end
end
