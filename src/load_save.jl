using JLD

function save_output(file_path::AbstractString,output::DIC_Output)
    jldopen(file_path, "w") do file
        addrequire(file, DynamicPolynomials)
        write(file, "output", output)
    end
end