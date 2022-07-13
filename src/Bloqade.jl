module Bloqade

using Reexport

@reexport using Yao
@reexport using BloqadeMIS
@reexport using BloqadeODE
@reexport using BloqadeExpr
@reexport using BloqadeKrylov
@reexport using YaoSubspaceArrayReg
@reexport using BloqadeWaveforms

# partially reexport
@reexport using Measurements: ±, Measurement
@reexport using BloqadeLattices:
    BloqadeLattices,
    AbstractLattice,
    GeneralLattice,
    HoneycombLattice,
    SquareLattice,
    TriangularLattice,
    ChainLattice,
    LiebLattice,
    KagomeLattice,
    GeneralLattice,
    RectangularLattice,
    AtomList,
    # interfaces
    generate_sites,
    offset_axes,
    rescale_axes,
    MaskedGrid,
    make_grid,
    random_dropout,
    clip_axes,
    lattice_sites,
    lattice_vectors,
    make_kdtree,
    grouped_nearest,
    collect_atoms

export rydberg_density, rydberg_corr, bitstring_hist, bitstring_hist!

using PythonCall
const plt = PythonCall.pynew()
const braket = PythonCall.pynew()

# println("wtf_outside")
# using Pkg
# Pkg.add("CondaPkg")
# using CondaPkg
# CondaPkg.add("amazon-braket-sdk")

# [conda]
# packages = ["python>=3.7.2", "amazon-braket-sdk"]
# channels = ["conda-forge"]


function __init__()
    # copied from PyPlotCall.jl
    println("wtf")
    PythonCall.pycopy!(plt, pyimport("matplotlib.pyplot"))
    PythonCall.pycopy!(braket, pyimport("amazon-braket-sdk"))
    return
end

include("plots/plots.jl")
include("observables.jl")
include("precompile.jl")

end
