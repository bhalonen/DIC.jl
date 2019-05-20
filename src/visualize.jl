using FixedPointNumbers

make_heat_map(image::Matrix{T}, polynomial::Polynomial, time_frame::Real, result::DIC_Output) where T<:Gray = make_heat_map(size(image), polynomial, time_frame, result)

function make_heat_map(size_image::Tuple, polynomial::Polynomial, time_frame::Real, result::DIC_Output) where T<:Gray
    [in_roi(idxY, idxX, time_frame, result) ? polynomial(idxX,idxY,time_frame) : NaN
        for idxY in 1:size_image[1], idxX in 1:size_image[2]
    ]
end
