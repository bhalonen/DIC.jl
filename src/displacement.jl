function __init__()
    global sobol_seq = SobolSeq(2)
end
struct Sample 
    point::CartesianIndex
    frame::Integer
    local_area_coords::Vector{CartesianIndex}
    local_area::Vector{<:Real}
end 
function get_local_area_samples(sample::CartesianIndex{2})
    vec([CartesianIndex(sample[1]+idxY, sample[2]+idxX) for idxY in -3:3, idxX in -3:3])
end
function get_roi_point(roi::Rect_ROI)
    sample = next!(sobol_seq)
    x = roi.corner1.x + Int(floor(sample[1]*(roi.corner2.x-roi.corner1.x)))
    y = roi.corner1.y + Int(floor(sample[2]*(roi.corner2.y-roi.corner1.y)))
    Point(y,x)
end
function get_sample_point(roi::Rect_ROI, original_image::Matrix{<:Real}, frame::Integer)
    point = get_roi_point(roi)
    center_sample = CartesianIndex(point.y,point.x)
    grid_samples = get_local_area_samples(center_sample)
    original_values = map(sample->original_image[sample],grid_samples)
    return Sample(center_sample, frame, grid_samples, original_values)
end
function DIC_analysis(dic_input::DIC_Input)
    reference_image = dic_input.images[1]
    deformed_images = dic_input.images[2:end]
    initial_guess_u = (rand(length(dic_input.dic_run_params.u_model)).-.5)/2000
    initial_guess_v = (rand(length(dic_input.dic_run_params.v_model)).-.5)/2000
    original_image = map(reference_image) do px
        Float32(px.val)
    end 
    roi_samples=map(1:dic_input.dic_run_params.dic_setting.sample_count) do idx
        get_sample_point(dic_input.roi,original_image,Int32(idx%(length(deformed_images))+1))
    end
    
    deformed_images_plain = map(deformed_images) do image
            size(reference_image)==size(image) || error("all input images must be the same size $(size(reference_image)) versus $(size(image))")
            map(image) do px
                Float32(px.val)
            end 
        end 
    deformed_images_itps = map(image->interpolate(image, BSpline(Cubic(Line(OnGrid())))),deformed_images_plain)
    cost_function(polynomial_coeff) = cost_function_builder(deformed_images_itps,
                                    roi_samples,
                                    dic_input.dic_run_params,
                                    polynomial_coeff,
                                    dic_input.time_table)
    @time result = optimize(cost_function, cat(dims=1,initial_guess_u,initial_guess_v), NelderMead())
    u_params=Optim.minimizer(result)[1:length(initial_guess_u)]
    v_params=Optim.minimizer(result)[length(initial_guess_u)+1:end]
    DIC_Output{Lagrangian}(Polynomial(u_params,dic_input.dic_run_params.u_model),Polynomial(v_params,dic_input.dic_run_params.v_model),dic_input.roi)
end
function cost_function_builder( deformed_image_itps::Vector{<:AbstractInterpolation},
                                roi_samples::Vector{<:Sample},
                                dic_run_params::DIC_Run_Parameters, 
                                polynomial_coeff::Vector{<:Real},
                                time_table::Vector{<:Real})
    u = Polynomial(polynomial_coeff[1:length(dic_run_params.u_model)],dic_run_params.u_model)
    v = Polynomial(polynomial_coeff[length(dic_run_params.u_model)+1:end],dic_run_params.v_model)
    sum(roi_samples) do sample 
        get_local_distance(deformed_image_itps[sample.frame], sample, u, v, time_table[sample.frame+1])
    end
end 
function get_transformed_point(u::Polynomial,v::Polynomial,test_point::CartesianIndex{2},size_image::Tuple, time_value::Real)
    Δx = u(test_point[2], test_point[1], time_value)
    Δy = v(test_point[2], test_point[1], time_value)
    return position_bounded_to_image( Δy + test_point[1], Δx + test_point[2], size_image)
end 

function get_local_distance(image_itp::AbstractInterpolation, sample::Sample, u::Polynomial, v::Polynomial, time_value::Real) 
    transformed_values = map(sample.local_area_coords) do local_sample_point
        point = get_transformed_point(u, v, local_sample_point, size(image_itp), time_value)
        point.x == size(image_itp)[2] || point.y == size(image_itp)[1]|| point.x ==0 || point.y ==0 && return 99999999
        image_itp(point.y, point.x)
    end
    norm(transformed_values - sample.local_area)
end
