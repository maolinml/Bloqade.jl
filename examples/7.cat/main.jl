
# import pacakges 

using Bloqade
using PythonCall
using StatsBase
using Distributed
using BitBasis

matplotlib = pyimport("matplotlib")
plt = pyimport("matplotlib.pyplot")


ncores = Sys.CPU_THREADS
if nprocs()==1
    addprocs(ncores; exeflags="--project")
else
    println("You already have quite some processes")
end


# Experiment 1: Preparing the GHZ state
# Here we focus on an atomic chain of even sites, and generate the state of following (not sure why it is called GHZ..)

#  $
#  |\text{GHZ}⟩ = 1/√2(|0101...⟩ + |1010...⟩)
#  $

# (There seems no easy to generate such state if the number of sites is odd)

Δmin = -6 * 2*pi
Δmax = 10 * 2*pi
Ωmax = 5 * 2*pi
Ωmin = 0 

t1 = 0.5 # μs
t2 = 4.0 # μs
t3 = 0.5 # μs
tq = 0.25 # μs
dt = 1e-3 # μs

a = 5.5 # 7.0
N = 10 # 17

# Local detuning
δ1 = -4.5 * 2*pi # For the edges
δ2 = -1.5 * 2*pi # For the third site from the edges

δq = -3.8 * 2*pi


atoms = generate_sites(ChainLattice(), N, scale=a)



function get_global_pulse()
    
    
    Δ = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Δmin, Δmin, Δmax, Δmax])
    Ω = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Ωmin, Ωmax, Ωmax, Ωmin])
    
    return Δ, Ω
end

function get_local_pulse_1()
    
    Δ = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Δmin, Δmin, Δmax, Δmax])
    Ω = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Ωmin, Ωmax, Ωmax, Ωmin])
    

    Δlocal = map(1:length(atoms)) do idx
                 if idx == 1 || idx == N
                     constant(duration=t1+t2+t3, value=δ1)
                 elseif idx == 3 || idx == N-2
                     constant(duration=t1+t2+t3, value=δ2)
                 else
                     constant(duration=t1+t2+t3, value=0)
                 end
             end        
    
    Δ = Δ .+ Δlocal
    
    return Δ, Ω
end



function getbarplot(Δ, Ω)
    h = rydberg_h(atoms; Δ=Δ, Ω=Ω)
    reg = zero_state(length(atoms))

    # We can then simulate the time evolution of the quantum state using an ODE solver:
    duration = Ω.duration
    prob = SchrodingerProblem(reg, duration, h);
    integrator = init(prob, Vern8());
    # Then, we measure the real-time expectation value of the Rydberg density and entanglement entropy: 

    # entropy = Float64[]
    densities = []
    for _ in TimeChoiceIterator(integrator, 0.0:1e-3:duration)
        push!(densities, rydberg_density(reg))
        # rho = density_matrix(reg, (1,2,3,4,5)) # calculate the reduced density matrix
        # push!(entropy, von_neumann_entropy(rho)) # compute entropy from the reduced density matrix
    end

    bitstring_hist(reg; nlargest = 20)

end


Δ, Ω = get_global_pulse()

fig, ax = plt.subplots(1, 1, figsize = (10,4))
Bloqade.plot!(ax, Ω)
Bloqade.plot!(ax, Δ)
ax.grid()

ax.legend(["Ω", "Δ"])

fig


Δ, Ω = get_local_pulse_1()

fig, ax = plt.subplots(1, 1, figsize = (10,4))
for Δp in Δ
    Bloqade.plot!(ax, Δp)
end

ax.grid()

ax.legend(range(1,length(Δ)))

fig


## Experiment 2: Control and measure the phase of the GHZ state

# We shall measure the phase in the GHZ state

#  $
#  |\text{GHZ}⟩ = 1/√2(|0101...⟩ + e^{i\phi}|1010...⟩)
#  $


# function get_local_pulse_2(tq, δq)
    
#     Δglobal = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Δmin, Δmin, Δmax, Δmax])    
#     Δ = []
#     tq2 = pi/(sqrt(2)*Ωmax)
#     for i in range(1, N)
#         if i==1
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ1)
#             Δtemp = append(Δtemp, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==N
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ1)
#             Δtemp = append(Δtemp, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==3
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ2)
#             Δtemp = append(Δtemp, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==N-2
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ2)
#             Δtemp = append(Δtemp, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                        
#         elseif i % 2 == 1
#             Δtemp = append(Δglobal, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                        
#         elseif i % 2 == 0            
#             Δtemp = append(Δglobal, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                                    
#         end
#         push!(Δ, Δtemp)
#     end
    
    
#     Ω = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Ωmin, Ωmax, Ωmax, Ωmin])    
    
#     Ω = append(append(Ω, constant(duration=tq, value=0)), constant(duration=tq2, value=Ωmax))
    
#     return Δ, Ω
# end


function get_local_pulse_2(tq, δq)
    tq2 = pi/(sqrt(2)*Ωmax)    
    Δ1 = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Δmin, Δmin, Δmax, Δmax])    
    Δ1 = append(Δ1, constant(duration=tq+tq2, value=0))
    
    Δ2 = []
    Δ3 = []
    Δ4 = []
    Δ5 = []    
    
    for i in range(1, N)
        if i==1 || i==N
            push!(Δ2, append(constant(duration=t1+t2+t3, value=δ2), constant(duration=tq+tq2, value=0)))
        else
            push!(Δ2, constant(duration=t1+t2+t3+tq+tq2, value=0))
        end
        
        if i==3 || i==N-2
            push!(Δ3, append(constant(duration=t1+t2+t3, value=δ2), constant(duration=tq+tq2, value=0)))
        else
            push!(Δ3, constant(duration=t1+t2+t3+tq+tq2, value=0))
        end
        
        if i % 2 == 1
            push!(Δ4, append(
                    constant(duration=t1+t2+t3, value=0), 
                    constant(duration=tq, value=δq),
                    constant(duration=tq2, value=0)                    
                    )
            )
        else
            push!(Δ5, append(
                    constant(duration=t1+t2+t3, value=0), 
                    constant(duration=tq, value=-δq),
                    constant(duration=tq2, value=0)                    
                    )
            )
        end
        
                
                
        
#         if i==1
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ1)
#             Δtemp = append(Δtemp, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==N
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ1)
#             Δtemp = append(Δtemp, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==3
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ2)
#             Δtemp = append(Δtemp, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))            
#         elseif i==N-2
#             Δtemp = Δglobal .+ constant(duration=Δglobal.duration, value=δ2)
#             Δtemp = append(Δtemp, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                        
#         elseif i % 2 == 1
#             Δtemp = append(Δglobal, constant(duration=tq, value=δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                        
#         elseif i % 2 == 0            
#             Δtemp = append(Δglobal, constant(duration=tq, value=-δq))
#             Δtemp = append(Δtemp, constant(duration=tq2, value=0))                                    
#         end
#         push!(Δ, Δtemp)
#     end
    
    
    Ω = piecewise_linear(clocks=[0.0, t1, t1+t2, t1+t2+t3], values=[Ωmin, Ωmax, Ωmax, Ωmin])    
    
    Ω = append(append(Ω, constant(duration=tq, value=0)), constant(duration=tq2, value=Ωmax))
    
    return Δ1, Δ2, Δ3, Δ4, Δ5, Ω
end

Δ, Ω = get_local_pulse_2(tq, δq)

fig, ax = plt.subplots(2,1, figsize = (10,4))

Bloqade.plot!(ax[0], Ω)

for Δp in Δ
# for Δp in [Δ[10]]
    Bloqade.plot!(ax[1], Δp)
end
ax[0].grid()
ax[1].grid()

ax[1].legend(range(1,length(Δ)))

fig


function get_parity(confs, freqs)
    """Calculate the parity for a given set of configurations and the corresponding frequencies,
    at a given time
    """
    parity = 0
    confs = collect(confs)
    for i in range(1, length(confs))
        numof1 = length(baddrs(confs[i]))
        parity += freqs[i] * (-1)^numof1
    end
    return parity
end



function get_counts(Δ, Ω, shots)
    h = rydberg_h(atoms; Δ=Δ, Ω=Ω)
    N = length(atoms)
    reg = zero_state(N)
    duration = Ω.duration

    prob = SchrodingerProblem(reg, duration, h);
    integrator = init(prob, Vern8());
    # Then, we measure the real-time expectation value of the Rydberg density and entanglement entropy: 

    distributions = []

    for _ in TimeChoiceIterator(integrator, 0.0:dt:duration)
        push!(distributions, measure(reg; nshots=shots))
    end
    startind = Int((t1+t2+t3)/dt)
    distributions = distributions[startind:end]
    
    conflist = []
    freqlist = []
    paritylist = []
    
    for dist in distributions # loop through the dist at each time 
        temp = countmap(dist) # The dict with different configs and frequencies
        temp = sort(temp, byvalue=true, rev=true)
        confs = keys(temp)
        freqs = values(temp)
                
        freqs = [freq/shots for freq in freqs]
        push!(conflist, confs)
        push!(freqlist, freqs)         
        
        parity = get_parity(confs, freqs)
        
        print(freqs, confs, parity, '\n')   

        push!(paritylist, parity)


    end    
    
    return conflist, freqlist, paritylist

end


shots = 200
conflist, freqlist, paritylist = get_counts(Δ, Ω, shots) ;