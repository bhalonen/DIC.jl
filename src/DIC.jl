module DIC
using DynamicPolynomials
using Optim
using Sobol
using Images
using DynamicPolynomials
using FixedPointNumbers
using ColorTypes
using LinearAlgebra: norm
using ProgressMeter: @showprogress
using Interpolations

include("typesDIC.jl")
include("utils.jl")
include("displacement.jl")
include("eulerian.jl")
include("visualize.jl")
include("strains.jl")

export roi_contains, DIC_Types
export Point, 
        ROI, 
        Rect_ROI, 
        DIC_Setting, 
        DIC_Input,
        DIC_Run_Parameters,
        Displacement,
        Perspective,
        Eulerian,
        DIC_Output,
        DIC_analysis,
        make_heat_map,
        find_eulerian,
        exx_heatmap,
        exy_heatmap
end # module
