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
    color_range = get_red_green_color_range()
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y

    min_val,max_val = get_extrema(polynomial,roi,time_frame)
    color_image = [idxX in x_range && idxY in y_range ? weighted_color_mean(.5,RGB{N0f8}(image[idxX,idxY]),get_color(color_range,min_val,max_val, polynomial(idxX,idxY,time_frame))) :
                        RGB{N0f8}(image[idxX,idxY])
        for idxX in 1:size(image)[1], idxY in 1:size(image)[2]
    ]
    return color_image
end
# will be replaced by reverse transform.
function build_transformed_roi(irregular_grid, size_image)
    roi_transformed = zeros(Bool,size_image)
    foreach(irregular_grid) do image_coordinate
        roi_transformed[Int(round(image_coordinate[1])),Int(round(image_coordinate[2]))] = true
    end
    erode(dilate(roi_transformed))
end 

function make_heat_map(image::Matrix{T}, polynomial::Polynomial, time_frame::Real, result::DIC_Output{Eulerian}) where T<:Gray
    roi = result.roi
    color_range = get_red_green_color_range()
    x_range = roi.corner1.x:roi.corner2.x
    y_range = roi.corner1.y:roi.corner2.y

    min_val,max_val = get_extrema(polynomial,roi,time_frame)
    irregular_grid = vec([position_bounded_to_image(position(idxX, idxY, time_frame, result)..., size(image)) for idxX in x_range, idxY in y_range ])
    values = vec([polynomial(idxX,idxY,time_frame) for idxX in x_range, idxY in y_range ])
    transformed_roi = build_transformed_roi(irregular_grid, size(image))
    transformed_interp = Spline2D(map(val->val.x,irregular_grid), map(val->val.y,irregular_grid), values; kx=3, ky=3, s=length(values))
    color_image = [
        transformed_roi[idxX,idxY] ?  weighted_color_mean(.5,RGB{N0f8}(image[idxX,idxY]),get_color(color_range,min_val,max_val, evaluate(transformed_interp,idxX,idxY))) :
                                     RGB{N0f8}(image[idxX,idxY])
        for idxX in 1:size(image)[1], idxY in 1:size(image)[2]
    ]
    return color_image
end