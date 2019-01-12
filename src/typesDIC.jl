struct Point
    x::Int32
    y::Int32
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
function displacement(x::Real, y::Real, t::Real, result::DIC_Output)
    (Δx = result.u_transform(x, y ,t), Δy = result.v_transform(x, y ,t)) 
end 
function position(x::Real, y::Real, t::Real, result::DIC_Output)
    transfrom = displacement(x,y,t,result)
    return (x = x + transfrom.Δx, y = y + transfrom.Δy)
end 
function position_bounded_to_image(x_position::Real, y_position::Real, size_image::Tuple)
    x_point = max(1,min(x_position, size_image[1]))
    y_point = max(1,min(y_position, size_image[2]))
    return (x = x_point, y = y_point)
end 

