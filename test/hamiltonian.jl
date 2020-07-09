using Test
using RydbergEmulator
using LightGraphs: SimpleGraph, add_edge!
using SparseArrays
using OrderedCollections
using CUDA
using LuxurySparse
using Yao
using Yao.ConstGate: P0, P1

function rydinteract(C, atoms)
    n = length(atoms)
    terms = []
    for i in 1:n, j in 1:n
        if i != j
            push!(terms, C/(2 * RydbergEmulator.distance(atoms[i], atoms[j])^6) * kron(n, i=>P1, j=>P1))
        end
    end
    return sum(terms)
end

@testset "simple graph hamiltonian subspace" begin
    subspace = Subspace(test_graph)
    @test collect(keys(subspace)) == sort(test_subspace_v)
    @test collect(values(subspace)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    # test hamiltonian creation
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]

    target = create_test_hamiltonian(Δ, Ω, ϕ)
    h = XTerm(Ω, ϕ) + ZTerm(Δ)
    @test SparseMatrixCSC(h, subspace) ≈ target

    H = SparseMatrixCSC(h, subspace)
    @test H ≈ update_term!(copy(H), h, subspace)

    @testset "cuda" begin
        if CUDA.functional()
            using CUDA
            using CUDA.CUSPARSE
            dΩ = cu(Ω)
            dϕ = cu(ϕ)
            dΔ = cu(Δ)
            h = XTerm(dΩ, dϕ) + ZTerm(dΔ)
            dH = CuSparseMatrixCSR(H)
            ds = cu(subspace)
            update_term!(dH, h, ds)
            @test isapprox(SparseMatrixCSC(dH), H; rtol=1e-7)
        end
    end
end

@testset "X term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.X) for k in 1:5]))
    h = XTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
end

@testset "Z term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.Z) for k in 1:5]))
    h = ZTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
end

@testset "rydberg interact term" begin
    atoms = RydbergEmulator.square_lattice(4, 0.8)
    H = mat(rydinteract(2.0, atoms))
    h = RydInteract(2.0, atoms)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
end

@testset "composite term" begin
    atoms = RydbergEmulator.square_lattice(5, 0.8)
    h = XTerm(5, 2.0) + RydInteract(2.0, atoms) + ZTerm(5, 1.0)
    H = rydinteract(2.0, atoms) +
        sum([1.0 * kron(5, k=>Yao.Z) for k in 1:5]) +
        sum([2.0 * kron(5, k=>Yao.X) for k in 1:5])
    H = mat(H)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
end

@testset "XTerm subspace" begin
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]
    h = XTerm(Ω, ϕ)

    subspace = Subspace(test_graph)
    H = SparseMatrixCSC(h, subspace)
    @test update_term!(copy(H), h, subspace) ≈ H
end
