module DIC
using DynamicPolynomials
using Optim
using Sobol
using ImageFeatures, Images
using DynamicPolynomials
using FixedPointNumbers
using ColorTypes
using LinearAlgebra: norm
using Distributed
using ProgressMeter: @showprogress
using Statistics: cor
using Interpolations

include("typesDIC.jl")
include("utils.jl")
include("displacement.jl")
include("visualize.jl")


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
        make_heat_map
end # module
