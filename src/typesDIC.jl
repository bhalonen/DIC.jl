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
abstract type Interpolation_Type end
struct Cubic_Keys<: Interpolation_Type end
struct Linear <: Interpolation_Type end 

abstract type Sub_Region end
struct Circle <: Sub_Region end 
struct Square <: Sub_Region 
end
struct DIC_Setting 
    sample_count::Int
end 
struct DIC_Run_Parameters
    scale_factor::Real
    sub_region::Sub_Region
    radius::Real
    dic_setting::DIC_Setting
    u_model::MonomialVector
    v_model::MonomialVector
end 
struct DIC_Input
    images::Vector{Matrix}
    roi::ROI
    dic_run_params::DIC_Run_Parameters
end 


struct Displacement
    u::Real
    v::Real
end 
abstract type Perspective end 
struct Eulerian <: Perspective end
struct Lapalacian <: Perspective end 
struct DIC_Output
    u::Polynomial 
    v::Polynomial
    perspective::Perspective
end

