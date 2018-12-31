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
struct DIC_Output
    u_frames::Polynomial
    v_frames::Polynomial
    perspective::Perspective
end

