using Colors
using FixedPointNumbers
function get_color(color_range, maximum::Real, minimum::Real, value::Real)
    maximum == minimum && return color_range[1]
    percentile = (value-minimum)/(maximum-minimum)
    return color_range[clamp(Int64(floor(percentile * (length(color_range)-1) )),1,length(color_range))]
end
function get_extrema(polynomial::Polynomial,roi::Rect_ROI,time_frame::Real)
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y
    values = vec([polynomial(idxX,idxY,time_frame) for idxX in x_range, idxY in y_range ])
    extrema(values)
end
function get_red_green_color_range()
    range(colorant"green",colorant"red")
end 
function make_heat_map(image::Matrix{T}, polynomial::Polynomial,roi::Rect_ROI, time_frame::Real) where T<:Gray
    color_range = get_red_green_color_range()
    color_image = map(px->RGB{N0f8}(px),image)
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y

    min_val,max_val = get_extrema(polynomial,roi,time_frame)
    for idxX in x_range, idxY in y_range
        color_image[idxX,idxY]=weighted_color_mean(.5,color_image[idxX,idxY],get_color(color_range,min_val,max_val, polynomial(idxX,idxY,time_frame)) )
    end
    return color_image
end 
function make_heat_map(image::Matrix{T}, polynomial::Polynomial, time_frame::Real, result::DIC_Output{Lagrangian}) where T<:Gray
    roi = result.roi
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y

    heat_map = [idxX in x_range && idxY in y_range ?  polynomial(idxX,idxY,time_frame) : NaN
        for idxX in 1:size(image)[1], idxY in 1:size(image)[2]
    ]
    return heat_map
end

function make_heat_map(image::Matrix{T}, polynomial::Polynomial, time_frame::Real, result::DIC_Output{Eulerian}) where T<:Gray

    heat_map = [
        in_roi(idxX,idxY,time_frame, result) ?  polynomial(idxX,idxY,time_frame) : NaN
        for idxX in 1:size(image)[1], idxY in 1:size(image)[2]
    ]
    return heat_map
end