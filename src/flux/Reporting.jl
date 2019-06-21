function generate_csv_flux_report(path_to_report_file::String,flux_distribution::Array{Float64,2},data_dictionary::Dict{String,Any}; number_of_cases::Int = 2)

    # TODO: check - does path_to_report_file point to a legit directory?

    # initialize -
    report_buffer = Array{String,1}()

    # how many fluxes do we have?
    (number_of_fluxes,number_of_cols) = size(flux_distribution)

    # growth string -
    growth_record = "$(number_of_fluxes+1),growth,0,inf"

    # get the chemical reaction strings -
    list_of_chemical_reaction_strings = data_dictionary["list_of_chemical_reaction_strings"]
    list_of_chemical_reaction_strings = [list_of_chemical_reaction_strings ; growth_record]

    # build -
    for flux_index = 1:number_of_fluxes

        # init blank line -
        line = ""
        line *= list_of_chemical_reaction_strings[flux_index]

        for col_index = 1:number_of_cols

            line *= ","
            line *= "$(flux_distribution[flux_index, col_index])"
        end

        # push -
        +(report_buffer, line)
    end

    # write -
    write_file_to_path(path_to_report_file, report_buffer)
end
