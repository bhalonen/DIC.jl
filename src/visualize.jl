using FixedPointNumbers

function make_heat_map(image::Matrix{T}, polynomial::Polynomial, time_frame::Real, result::DIC_Output) where T<:Gray
    [in_roi(idxY, idxX, time_frame, result) ? polynomial(idxX,idxY,time_frame) : NaN
        for idxY in 1:size(image)[1], idxX in 1:size(image)[2]
    ]
end
