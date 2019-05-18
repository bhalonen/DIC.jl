import Base: <=
<=(point1::Point,point2::Point) = point1.x<=point2.x && point1.y<=point2.y

function in_roi(roi::Rect_ROI,point::Point)
    roi.corner1 <= point && point <= roi.corner2
end 
function downscale_image(image::Matrix,scale_factor::Real)
    size_x,size_y=size(image)/scale_factor
    step_x=size(image)[1]/size_x
    step_y=size(image)[2]/size_y
    image[1:step_x:size(image)[1], 1:step_y:size(image)[2]]

end