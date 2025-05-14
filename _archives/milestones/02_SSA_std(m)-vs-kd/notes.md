- we’ve tried with kd=0.01 first but then tried kf=0.001 since it seemed the constructive reactions overpowered the diffusion
- then we’re doing a log range for kd as to test for very small values near zero
- I’m increasing tau_max from 100 --> 1000 since we’re not hitting the steady state (also trying to increase memory since this last batch failed)
- then we realize there’s a bug in JLD2, fixed by upgrading the module
- t=500 works, but t=1000 still fails (JLD2 file is really big...)
- then I’m trying 2way diffusion, but I think that in Ensemble.jl the find_inflow_nodes scrambles the lattice—trying the 1way diffusion for now



H_lattice-2way-diffusion

- I think I’ve fixed the thing with inflow_nodes
- trying again 2way diffusion



I_lattice-2way-diffusion-wider

- the STD of (H) was too large, trying again with a new set of parameters: changing the forward rate to 5e-5 and the number of reactors to 25, and the diffusion 1e-6...1e2
- result: parameter sweep too large, exceeded 4h computation limit



J_lattice-2way-diffusion-wider-T=100

- reducing T to try again (I)
- ==(still running...)==



K_lattice-2way-diffusion-wider-1e1-T=100

- reducing diffusion parameter range because it’s still too long to compute
- result: means & std zero everywhere, is this because forward rate = 5e-5? that’s the value cole had used in his PPT figure but maybe there’s something else going on



L_lattice-2way-diffusion-increased-kd

- well here I’m increasing again the forward rate to what I had set initially, still keeping a wider range for the diffusion rate coefficient
- results: interesting fig



note: for the 51 pegasi b proposal we’re using the figs from `L_lattice-2way-diffusion-increased-kd` and `03_AI-vs-C` (data from 02L)
