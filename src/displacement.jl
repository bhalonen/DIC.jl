
function __init__()
    global sobol_seq = SobolSeq(2)
end
function get_sample_point(roi::Rect_ROI)
    sample = next!(sobol_seq)
    x = roi.corner1.x + Int(floor(sample[1]*(roi.corner2.x-roi.corner1.x)))
    y = roi.corner1.y + Int(floor(sample[2]*(roi.corner2.y-roi.corner1.y)))
    return CartesianIndex(x,y)
end
function DIC_analysis(dic_input::DIC_Input)
    reference_image = dic_input.images[1]
    deformed_images = dic_input.images[2:end]
    initial_guess_u=zeros(length(dic_input.dic_run_params.u_model))
    initial_guess_v=zeros(length(dic_input.dic_run_params.v_model))
    size(reference_image)==size(deformed_images[1]) || error("all input images must be the same size $(size(reference_image)) versus $(size(deformed_image))")
    roi_samples=map(1:dic_input.dic_run_params.dic_setting.sample_count) do idx
        (Int32(idx%(length(deformed_images))+1),get_sample_point(dic_input.roi))
    end
    original_values = map(sample -> reference_image[sample[2][1],sample[2][2]].val,roi_samples)
    deformed_images_plain = map(deformed_images) do image
            map(image) do px
                Float32(px.val)
            end 
        end 
    deformed_images_itps = map(image->interpolate(image, BSpline(Cubic(Line(OnGrid())))),deformed_images_plain)
    cost_function(polynomial_coeff) = cost_function_builder(deformed_images_itps,
                                    roi_samples,
                                    dic_input.dic_run_params,
                                    polynomial_coeff,
                                    dic_input.time_table,
                                    original_values)
    @time result = optimize(cost_function, cat(dims=1,initial_guess_u,initial_guess_v), NelderMead())
    u_params=Optim.minimizer(result)[1:length(initial_guess_u)]
    v_params=Optim.minimizer(result)[length(initial_guess_u)+1:end]
    DIC_Output(Polynomial(u_params,dic_input.dic_run_params.u_model),Polynomial(v_params,dic_input.dic_run_params.v_model),Lagrangian)
end

function cost_function_builder( deformed_image_itps::Vector{<:AbstractInterpolation},
                                roi_samples::Vector{Tuple{Int32,CartesianIndex{2}}},
                                dic_run_params::DIC_Run_Parameters, 
                                polynomial_coeff::Vector{<:Real},
                                time_table::Vector{<:Real}, 
                                original_values::Vector{Normed{UInt8,8}})
    u = Polynomial(polynomial_coeff[1:length(dic_run_params.u_model)],dic_run_params.u_model)
    v = Polynomial(polynomial_coeff[length(dic_run_params.u_model)+1:end],dic_run_params.v_model)
    transformed_values = pmap(roi_samples) do sample
        get_transformed_value(deformed_image_itps[sample[1]],sample[2], u, v, time_table[sample[1]+1])
    end
    1-cor(original_values,transformed_values)
end 
function get_transformed_point(u::Polynomial,v::Polynomial,test_point::CartesianIndex{2},size_image::Tuple, time_value::Real)
    x_translation = u(test_point[1], test_point[2], time_value)
    y_translation = v(test_point[1], test_point[2], time_value)
    x_point = max(1,min(x_translation + test_point[1], size_image[1]))
    y_point = max(1,min(y_translation + test_point[2], size_image[2]))
    return (x = x_point, y = y_point)
end 
function get_transformed_value(image_itp::AbstractInterpolation, sample::CartesianIndex{2}, u::Polynomial, v::Polynomial,time_value::Real) 
    point = get_transformed_point(u,v,sample, size(image_itp),time_value)
    image_itp(point.x,point.y)
end
