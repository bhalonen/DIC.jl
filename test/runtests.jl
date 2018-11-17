using Distributed
#addprocs(6)
@everywhere using DIC
using Images
using FileIO
using ImageView
using Colors
using DynamicPolynomials

test_dir=string(@__DIR__)

images=map(readdir(joinpath(test_dir,"test_images"))) do image_file_name
    map(pixel->Gray(pixel),load(joinpath(test_dir,"test_images",image_file_name)))
end
images=images[1:3]
@polyvar x y 
u_model=MonomialVector([x*x,y*y,x*y,x,y,1])
v_model=MonomialVector([x*x,y*y,x*y,x,y,1])
roi=Rect_ROI(Point(150,50),Point(200,100))
dic_run_params=DIC_Run_Parameters(1,Square(),5 ,DIC_Setting(30), u_model, v_model)

@show result = DIC_analysis(DIC_Input(images, roi, dic_run_params))
#uncomment to see plots
#imshow(make_heat_map(images[2],result.u_frames[1], roi))

# mapped_roi = Matrix{Bool}(undef,size(images)...)

# for idxX in 1:size(mapped_roi)[1], idxY in 1:size(mapped_roi)[2]
#     mapped_roi[idxX,idxY] = roi_contains(roi,Point(idxX,idxY))
# end

