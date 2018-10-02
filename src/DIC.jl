module DIC
using DynamicPolynomials

include("typesDIC.jl")
include("utils.jl")
include("displacement.jl")


export roi_contains, DIC_Types
export Point, 
        ROI, 
        Rect_ROI, 
        Interpolation_Type, 
        Cubic_Keys, 
        Linear, 
        Sub_Region, 
        Circle, 
        Square,
        DIC_Setting, 
        No_Update, 
        Keep_Most_Points, 
        Remove_Bad_Points,
        DIC_Input,
        DIC_Run_Parameters,
        Displacement,
        Perspective,
        Eulerian,
        DIC_Output,
        DIC_analysis
end # module
