using Optim
using Sobol
using ImageFeatures, Images
using DynamicPolynomials
using FixedPointNumbers
using ColorTypes
function __init__()
    global sobol_seq = SobolSeq(2)
end
function DIC_analysis(dic_input::DIC_Input)
    length(dic_input.images)>=2 || error("must have more than 2 images")
    initial_guess_u=zeros(length(dic_input.dic_run_params.u_model))
    initial_guess_v=zeros(length(dic_input.dic_run_params.v_model))
    u_frames=Vector{Polynomial}()
    v_frames=Vector{Polynomial}()
    foreach(dic_input.images[2:end]) do image
        u,v=FuncDIC(dic_input.images[1],image,dic_input.roi,dic_input.dic_run_params,initial_guess_u,initial_guess_v)
        initial_guess_u = u.a 
        initial_guess_v = v.a 
        push!(u_frames,u)
        push!(v_frames,v)
    end
    DIC_Output(u_frames,v_frames,Lapalacian())
end
function FuncDIC(reference_image::Matrix{T},deformed_image::Matrix{T},roi::Rect_ROI,dic_run_params::DIC_Run_Parameters,initial_guess_u,initial_guess_v) where T<:Gray
    size(reference_image)==size(deformed_image) || error("all input images must be the same size $(size(reference_image)) versus $(size(deformed_image))")
    dic_run_params.radius<5 && error("radius must be 5 or greater")
    cost_function(polynomial_coeff) = cost_function_builder(reference_image,deformed_image,roi,dic_run_params,polynomial_coeff)
    result = optimize(cost_function, cat(dims=1,initial_guess_u,initial_guess_v), BFGS())
    u_params=Optim.minimizer(result)[1:length(initial_guess_u)]
    v_params=Optim.minimizer(result)[length(initial_guess_u)+1:end]
    Polynomial(u_params,dic_run_params.u_model),Polynomial(v_params,dic_run_params.v_model)
end
function get_sample_point(roi::Rect_ROI)
    sample = next!(sobol_seq)
    x = roi.corner1.x + Int(floor(sample[1]*(roi.corner2.x-roi.corner1.x)))
    y = roi.corner1.y + Int(floor(sample[2]*(roi.corner2.y-roi.corner1.y)))
    return Point(x,y)
end
function get_transformed_point(u::Polynomial,v::Polynomial,test_point::Point,size_image::Tuple, radius::Int)
    x_translation = u(test_point.x, test_point.y)
    y_translation = v(test_point.x, test_point.y)
    x_point = max(radius+1, min(size_image[1]-radius-1, test_point.x + Int(floor(x_translation))))
    y_point = max(radius+1, min(size_image[2]-radius-1, test_point.y + Int(floor(y_translation))))
    return Point(x_point, y_point)
end 
function cost_function_builder(reference_image::Matrix{T},deformed_image::Matrix{T},roi::Rect_ROI,dic_run_params::DIC_Run_Parameters, polynomial_coeff::Vector) where T<:Gray
    u = Polynomial(polynomial_coeff[1:length(dic_run_params.u_model)],dic_run_params.u_model)
    v = Polynomial(polynomial_coeff[length(dic_run_params.u_model)+1:end],dic_run_params.v_model)
    # parallelize here
    sum(1:dic_run_params.dic_setting.sample_count) do idx
      test_point = get_sample_point(roi)
      compare_point = get_transformed_point(u,v,test_point, size(reference_image), dic_run_params.radius)
        orb_similarity(reference_image,deformed_image,test_point,compare_point,dic_run_params.radius*2)
    end
end 
function get_sub_image(image::Matrix{T},center_point::Point, size::Int) where T<:Gray
    image[center_point.x-size÷2:center_point.x+size÷2,center_point.y-size÷2:center_point.y+size÷2]
end 
function orb_similarity(reference_image::Matrix{T},deformed_image::Matrix{T},test_point::Point,compare_point::Point,size::Int) where T<:Gray
    sub_img_orig=get_sub_image(reference_image,test_point,size)
    sub_img_compare = get_sub_image(deformed_image,compare_point,size)
    brief_params = BRIEF(size = 8, window = 3, seed = 123)
    keypoints_1 = Keypoints(fastcorners(sub_img_orig, 1, 0.35))
    keypoints_2 = Keypoints(fastcorners(sub_img_orig, 1, 0.35))
    desc_1, ret_keypoints_2 = create_descriptor(sub_img_orig, keypoints_1, brief_params)
    desc_2, ret_keypoints_2 = create_descriptor(sub_img_compare, keypoints_1, brief_params)
    hamming_distance(desc_1[1], desc_2[1])
end
