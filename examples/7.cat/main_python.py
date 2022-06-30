import numpy as np
pi = np.pi
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import multiprocessing  
import tqdm
import time

import qutip
from pulser import Pulse, Sequence, Register
from pulser.simulation import Simulation
from pulser.waveforms import RampWaveform, ConstantWaveform, CompositeWaveform
from pulser.devices import Chadoq2
from pulser.devices import MockDevice
backend = MockDevice # We use MockDevice to avoid some constraints in Chadoq2
sampling_rate = 1
shots = 5000
nsteps = 3000


Δmin = -6 * 2*pi
Δmax = 10 * 2*pi
Ωmax = 5 * 2*pi
Ωmin = 0 

t1 = 0.5 * 1e3
t2 = 4.0 * 1e3 # 4.0
t3 = 0.5 * 1e3

a = 5.5 # 7.0
N = 10 # 17

# Local detuning
δ1 = -4.5 * 2*pi # For the edges
δ2 = -1.5 * 2*pi # For the third site from the edges



def get_global_pulse(N):
    reg = Register.rectangle(1, N, spacing=a)
    seq = Sequence(reg, backend)
    seq.declare_channel('ising', 'rydberg_global')    

    seq1 = Pulse.ConstantDetuning(RampWaveform(t1, Ωmin, Ωmax), Δmin, 0.0) 
    seq2 = Pulse.ConstantAmplitude(Ωmax, RampWaveform(t2, Δmin, Δmax), 0.0)   
    seq3 = Pulse.ConstantDetuning(RampWaveform(t3, Ωmax, Ωmin), Δmax, 0.0) 
    seq.add(seq1, 'ising')   
    seq.add(seq2, 'ising')
    seq.add(seq3, 'ising')
    
    return seq

def get_counts(seq, sampling_rate, shots, nsteps):
    """
    Run the pulse sequence and collect the results
    """
    simul = Simulation(seq, sampling_rate=sampling_rate)
    results = simul.run(nsteps=nsteps)
    count = results.sample_final_state(N_samples=shots)

#     # If the backend is a simulator, then we could do
#     results = Experiment(seq).run(shots=shots)
#     count = results.count()
    
    counts = {k:v for k,v in count.items()}
#     print(sum(counts.values()))
    freq = list(counts.values()) / sum(counts.values())
    bits = list(counts.keys())
    
    conf = [[int(i) for i in item] for item in bits]
    ind = np.argsort(freq)[::-1]
    freq = freq[ind] # The most frequent state in the front
    conf = [conf[i] for i in ind]
    bits = [bits[i] for i in ind]
    return conf, freq

def plot_hist(conf, freq):

    plt.figure(figsize=(15,4))
    conf2 = [ "".join([str(i) for i in item]) for item in conf  ]
    plt.bar(conf2, freq)
    plt.xlabel("Bit string", fontsize=22)
    plt.ylabel("Probability", fontsize=22)    
    plt.xticks(rotation=45)
    plt.show() 

seq_global = get_global_pulse(N)
seq_global.draw()


conf, freq = get_counts(seq_global, sampling_rate, shots, nsteps)
plot_hist(conf, freq)


def get_pulse_GHZ(δ1, δ2, N):
    seq = get_global_pulse(N)
    
    seq.declare_channel('local1', 'rydberg_local', initial_target = [0, N-1] )    
    seq.declare_channel('local2', 'rydberg_local', initial_target = [2, N-3] )    
         
    seq.add(Pulse.ConstantPulse(t1+t2+t3, 0, δ1, 0.0), 'local1', protocol = 'no-delay')
    seq.add(Pulse.ConstantPulse(t1+t2+t3, 0, δ2, 0.0), 'local2', protocol = 'no-delay')    
    return seq
seq_GHZ = get_pulse_GHZ(δ1, δ2, N)
seq_GHZ.draw()



conf, freq = get_counts(seq_GHZ, sampling_rate, shots, nsteps)
plot_hist(conf, freq)


δp = -3.8 * 2*pi
shots = 200
def get_pulse_for_phase(t, δ1, δ2, N):
    seq = get_pulse_GHZ(δ1, δ2, N) # Get the pulse for the state preparation. 
  
    # Apply the stagger field
    t = int(t)    
    seq.declare_channel('local3', 'rydberg_local', initial_target = range(0,N,2) ) 
    seq.declare_channel('local4', 'rydberg_local', initial_target = range(1,N,2) )        
    seq.add(Pulse.ConstantPulse(t, 0, δp, 0.0), 'local3')
    seq.add(Pulse.ConstantPulse(t, 0, -δp, 0.0), 'local4')
    
    # Apply Ux
    tq2 = int(pi/(np.sqrt(2)*Ωmax) * 1e3)
    seq5 = Pulse.ConstantPulse(tq2, Ωmax, 0, 0.0)
    seq.add(seq5, 'ising')     
    
    return seq

def job(para):
    t, N = para[0], para[1]
    seq = get_pulse_for_phase(t, δ1, δ2, N)
    conf, freq = get_counts(seq, sampling_rate, shots, nsteps)    
    return conf, freq

def get_parity(conf, freq):
    res = 0
    for i in range(len(conf)):
        print(sum(conf[i]))        
        res += freq[i] * (-1)**(sum(conf[i]))
        
    return res

    
tq = 0.25 * 1e3
seq = get_pulse_for_phase(tq, δ1, δ2, N)
seq.draw()


Nrange = [4]
trange = np.linspace(0, tq, 101)[1:]
paras = [(t, N) for N in Nrange for t in trange]

confs, freqs = [], []

for para in paras:
    conf, freq = job(para)
    confs.append(conf)
    freqs.append(freq)


phases = [get_parity(conf, freq) for conf, freq in zip(confs, freqs)]


plt.plot(trange, phases)