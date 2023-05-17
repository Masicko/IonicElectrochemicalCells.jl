abstract type AbstractCell end

mutable struct CommonCell <: AbstractCell
    parameters#::ModelParameters.Model
    system#::VoronoiFVM.AbstractSystem
    U#::VoronoiFVM.AbstractUnknown
    Uold#::VoronoiFVM.AbstractUnknown
end

function biastest!(cell :: AbstractCell)
    for bias ∈ collect(0.0:-0.1:-1.0)
        stationary_update!(cell , Dict(:bias => bias))
    end
end

