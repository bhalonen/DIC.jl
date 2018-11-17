using Colors
using FixedPointNumbers
function get_color(color_range, maximum::Real, minimum::Real, value::Real)
    percentile = (value-minimum)/(maximum-minimum)
    return color_range[Int64(floor(percentile * (length(color_range)-1) )) + 1]
end
function get_extrema(polynomial::Polynomial,roi::Rect_ROI)
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y
    values = Vector{Float64}()
    for idxX in x_range, idxY in y_range
        push!(values,polynomial(idxX,idxY))
    end
    extrema(values)
end
function make_heat_map(image::Matrix{T}, polynomial::Polynomial,roi::Rect_ROI) where T<:Gray
    green = colorant"green"
    red = colorant"red"
    color_range = range(green,red)
    color_image = map(px->RGB{N0f8}(px),image)
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y

    min_val,max_val = get_extrema(polynomial,roi)
    for idxX in x_range, idxY in y_range
        color_image[idxX,idxY]=weighted_color_mean(.5,color_image[idxX,idxY],get_color(color_range,min_val,max_val, polynomial(idxX,idxY)) )
    end
    return color_image
end 