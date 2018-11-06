using Optim
using Sobol
using ImageFeatures, Images
using DynamicPolynomials
using FixedPointNumbers
using ColorTypes
using LinearAlgebra: norm
using Distributed
using ProgressMeter: @showprogress
function __init__()
    global sobol_seq = SobolSeq(2)

end
function DIC_analysis(dic_input::DIC_Input)
    length(dic_input.images)>=2 || error("must have more than 2 images")
    initial_guess_u=zeros(length(dic_input.dic_run_params.u_model))
    initial_guess_v=zeros(length(dic_input.dic_run_params.v_model))
    u_frames=Vector{Polynomial}()
    v_frames=Vector{Polynomial}()
    roi_samples=map(idx->get_sample_point(dic_input.roi),1:dic_input.dic_run_params.dic_setting.sample_count)
    @showprogress map(dic_input.images[2:end]) do image
        @time u,v=FuncDIC(dic_input.images[1],image,roi_samples,dic_input.dic_run_params,initial_guess_u,initial_guess_v)
        initial_guess_u = u.a 
        initial_guess_v = v.a 
        push!(u_frames,u)
        push!(v_frames,v)
        nothing
    end
    DIC_Output(u_frames,v_frames,Lagrangian)
end
function FuncDIC(reference_image::Matrix{T},deformed_image::Matrix{T},roi_samples::Vector{CartesianIndex{2}},dic_run_params::DIC_Run_Parameters,initial_guess_u,initial_guess_v) where T<:Gray
    size(reference_image)==size(deformed_image) || error("all input images must be the same size $(size(reference_image)) versus $(size(deformed_image))")
    dic_run_params.radius<5 && error("radius must be 5 or greater")
    cost_function(polynomial_coeff) = cost_function_builder(reference_image,deformed_image,roi_samples,dic_run_params,polynomial_coeff)
    result = optimize(cost_function, cat(dims=1,initial_guess_u,initial_guess_v), NelderMead())
    u_params=Optim.minimizer(result)[1:length(initial_guess_u)]
    v_params=Optim.minimizer(result)[length(initial_guess_u)+1:end]
    Polynomial(u_params,dic_run_params.u_model),Polynomial(v_params,dic_run_params.v_model)
end
function get_sample_point(roi::Rect_ROI)
    sample = next!(sobol_seq)
    x = roi.corner1.x + Int(floor(sample[1]*(roi.corner2.x-roi.corner1.x)))
    y = roi.corner1.y + Int(floor(sample[2]*(roi.corner2.y-roi.corner1.y)))
    return CartesianIndex(x,y)
end
function get_transformed_point(u::Polynomial,v::Polynomial,test_point::CartesianIndex,size_image::Tuple, radius::Int)
    x_translation = u(test_point[1], test_point[2])
    y_translation = v(test_point[1], test_point[2])
    x_point = max(radius+1, min(size_image[1]-radius-1, test_point[1] + Int(floor(x_translation))))
    y_point = max(radius+1, min(size_image[2]-radius-1, test_point[2] + Int(floor(y_translation))))
    return CartesianIndex(x_point, y_point)
end 
function cost_function_builder(reference_image::Matrix{T},deformed_image::Matrix{T},roi_samples::Vector{CartesianIndex{2}},dic_run_params::DIC_Run_Parameters, polynomial_coeff::Vector) where T<:Gray
    u = Polynomial(polynomial_coeff[1:length(dic_run_params.u_model)],dic_run_params.u_model)
    v = Polynomial(polynomial_coeff[length(dic_run_params.u_model)+1:end],dic_run_params.v_model)
    # parallelize here
    @distributed (+) for idx =1:length(roi_samples)
        test_point = roi_samples[idx]
        compare_point = get_transformed_point(u,v,test_point, size(reference_image), dic_run_params.radius)
        orb_similarity(reference_image,deformed_image,test_point,compare_point,dic_run_params.radius*2)
    end
end 
function get_sub_image(image::Matrix{T},center_point::CartesianIndex, size::Int) where T<:Gray
    image[center_point[1]-size÷2:center_point[1]+size÷2,center_point[2]-size÷2:center_point[2]+size÷2]
end 
function orb_similarity(reference_image::Matrix{T},deformed_image::Matrix{T},test_point::CartesianIndex,compare_point::CartesianIndex,size_box::Int) where T<:Gray
    sub_img_orig = get_sub_image(reference_image,test_point,size_box)
    sub_img_compare = get_sub_image(deformed_image,compare_point,size_box)
    edges_orig,count_orig = imhist(lbp(sub_img_orig, lbp_rotation_invariant),40)
    edges_compare,count_compare = imhist(lbp(sub_img_compare, lbp_rotation_invariant),edges_orig)
    norm(count_orig-count_compare)
end
