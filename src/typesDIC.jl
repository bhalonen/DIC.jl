struct Point
    y::Real
    x::Real
end 
abstract type ROI end
struct Rect_ROI<:ROI
    corner1::Point
    corner2::Point
end 
struct Image_ROI<:ROI
    map::Matrix{Bool}
end 
struct DIC_Setting 
    sample_count::Int
end 
struct DIC_Run_Parameters
    scale_factor::Real
    dic_setting::DIC_Setting
    u_model::MonomialVector
    v_model::MonomialVector
end 
struct DIC_Input{T}
    images::Vector{Matrix}
    time_table::Vector{T}
    roi::ROI
    dic_run_params::DIC_Run_Parameters
end 
@enum Perspective begin
    Eulerian
    Lagrangian
end 
struct DIC_Output{T}
    u_transform::Polynomial
    v_transform::Polynomial
    roi::ROI
end
function displacement( y::Real, x::Real, t::Real, result::DIC_Output)
    (Δx = result.u_transform(x, y ,t), Δy = result.v_transform(x, y ,t)) 
end 

function position( y::Real, x::Real, t::Real, result::DIC_Output)
    transform = displacement(y, x, t, result)
    return Point( y + transform.Δy, x + transform.Δx)
end 
function in_roi(y::Real, x::Real, t::Real, result::DIC_Output{Lagrangian})
    in_roi(y, x, result)
end
function in_roi( y::Real, x::Real, result::DIC_Output{Lagrangian})
    in_roi(result.roi, Point(y,x))
end
function in_roi(y::Real, x::Real, t::Real, result::DIC_Output{Eulerian})
    in_roi(result.roi, position(y, x, t, result))
end

function position_bounded_to_image( y_position::Real, x_position::Real, size_image::Tuple)
    x_point = max(1,min(x_position, size_image[2]))
    y_point = max(1,min(y_position, size_image[1]))
    return (x = x_point, y = y_point)
end 

