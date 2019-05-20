using Distributed
using DIC
using Images
using FileIO
using ImageView
using Colors
using DynamicPolynomials
using Plots 
plotlyjs()



test_dir=string(@__DIR__)

images=map(readdir(joinpath(test_dir,"test_images"))) do image_file_name
    map(pixel->Gray(pixel),load(joinpath(test_dir,"test_images",image_file_name)))
end
# images=images[1:3]
@polyvar x y t
u_model=monomials((x, y, t), [0, 1,2, 3,4], mono -> (degree(mono,x)<4 && degree(mono,y)<4 && 0<degree(mono,t)<3))#MonomialVector([x*x*t,x*x*t*t,y*y*t,y*y*t*t,x*y*t,x*y*t*t,x*t,y*t,t,1])
v_model=monomials((x, y, t), [0, 1,2, 3,4], mono -> (degree(mono,x)<4 && degree(mono,y)<4 && 0<degree(mono,t)<3))#MonomialVector([x*x*t,x*x*t*t,y*y*t,y*y*t*t,x*y*t,x*y*t*t,x*t,y*t,t,1])
#monomials((x, y, t), [0, 1,2, 3,4])
roi=Rect_ROI(Point(150,50),Point(200,100))
dic_run_params=DIC_Run_Parameters(1,DIC_Setting(800), u_model, v_model)
time_table = collect(0:(length(images)-1))
@show result = DIC_analysis(DIC_Input{Float32}(images, time_table, roi,  dic_run_params))
#uncomment to see plots
eulerian_result = find_eulerian(result, time_table)
plot(images[7],aspect_ratio=1, axis = false,size=(800,800))
heat_map_exx = make_heat_map(images[7],eulerian_result.u_transform, time_table[7], eulerian_result)
heatmap!(heat_map_exx, alpha=.5, clims=extrema(filter(!isnan,heat_map_exx)), colorbar=true,colorbar_title="u")
plot(images[11],aspect_ratio=1, axis = false,size=(800,800))
heat_map_exx = make_heat_map(images[11],result.v_transform, time_table[11], result)
heatmap!(heat_map_exx, alpha=.5, clims=extrema(filter(!isnan,heat_map_exx)), colorbar=true,colorbar_title="v")
# plot(images[1],aspect_ratio=1, axis = false,size=(800,800))
# heat_map_exy = exy_heatmap(size(images[1]),  result, time_table[1], x, y)
# heatmap!(heat_map_exy, alpha=.5, clims=extrema(filter(!isnan,heat_map_exy)), colorbar=true,colorbar_title="v")

# imshow(make_heat_map(images[7],result.v_transform, time_table[7], result))
# imshow(make_heat_map(images[7],eulerian_result.v_transform, time_table[7], eulerian_result))

# mapped_roi = Matrix{Bool}(undef,size(images)...)

# for idxX in 1:size(mapped_roi)[1], idxY in 1:size(mapped_roi)[2]
#     mapped_roi[idxX,idxY] = roi_contains(roi,Point(idxX,idxY))
# end

plot(images[11],aspect_ratio=1, axis = false,size=(800,800))
heat_map_exy = exy_heatmap(size(images[11]),  eulerian_result,time_table[11], x, y)
heatmap!(heat_map_exy, alpha=.5, clims=extrema(filter(!isnan,heat_map_exy)), colorbar=true,colorbar_title="v")
