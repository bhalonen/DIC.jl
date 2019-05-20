
function exx_heatmap(size_image::Tuple, result::DIC_Output, time_frame,x, y)
    @show e_xx = differentiate(result.u_transform,x) + (differentiate(result.u_transform,x)^2 + differentiate(result.v_transform,x)^2)/2
    make_heat_map(size_image,e_xx, time_frame, result)
end
function eyy_heatmap(size_image::Tuple, result::DIC_Output, time_frame,x, y)
    e_yy = differentiate(result.v_transform,y) + (differentiate(result.u_transform,y)^2 + differentiate(result.v_transform,y)^2)/2
    make_heat_map(size_image, e_yy, time_frame, result)
end
function exy_heatmap(size_image::Tuple, result::DIC_Output, time_frame, x, y)
    e_xy =  (differentiate(result.u_transform,y) +
        differentiate(result.v_transform,x) +
        differentiate(result.u_transform,y)*differentiate(result.u_transform,x) + 
        differentiate(result.v_transform,y)*differentiate(result.v_transform,x))/2
    make_heat_map(size_image,differentiate(result.v_transform,y), time_frame, result)
end