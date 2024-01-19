'''
Example of using the nuSQuIDSDecoh class to calculate neutrino transition probabilties resulting from neutrino decoherence

A simple example showing the popagation of a 1 GeV neutrino over distance, with decoherence effects damping the oscillations over distance

Tom Stuttard 
'''

import matplotlib
# matplotlib.use('AGG') # Might need this depending on your installation (some matplotlib installation hang without this)
import matplotlib.pyplot as plt

import numpy as np
import nuSQuIDS as nsq

# Get nuSQuiDS units
units = nsq.Const()

# Create nuSQuIDSDecoh instance (single energy mode)
nusquids = nsq.nuSQUIDSDecoh(
    3, # Num neutrino states
    nsq.NeutrinoType.neutrino, # nu  (as opposed to nubar)
)

# Solver settings
error = 1.e-8
nusquids.Set_rel_error(error)
nusquids.Set_abs_error(error)
nusquids.Set_ProgressBar(False)

# Set standard osc params
nusquids.Set_MixingAngle( 0, 1, np.deg2rad(33.) )
nusquids.Set_MixingAngle( 0, 2, np.deg2rad(8.5) )
nusquids.Set_MixingAngle( 1, 2, np.deg2rad(49.) )
nusquids.Set_SquareMassDifference( 1, 7.5e-5*units.eV*units.eV )
nusquids.Set_SquareMassDifference( 2, 2.5e-3*units.eV*units.eV )
nusquids.Set_CPPhase( 0, 2, np.deg2rad(0.) )

# Solving in vacuum
nusquids.Set_Body(nsq.Vacuum())

# Define neutrino energy and path
E_GeV = 1.
L_km = 8000.
print(f"\nE = {E_GeV} GeV, L = {L_km} km")

# Pass E to nuSQuIDS
nusquids.Set_E(E_GeV*units.GeV)

# Define decoherence model
n = 0 # Energy-indepdent damping coefficients
E0_eV = 1e9 # 1 GeV (a common choice, but is arbitrary)
gamma0_eV = 1e-13 # Strength of the damping
D_eV = np.diag( [0.] + (8*[gamma0_eV] )) # "State selection" model from arxiv:2007.00068
print(f"\nGamma0 = {gamma0_eV} eV, E0 = {E0_eV} eV")
print("\nD matrix [eV]:")
print(D_eV)


# Create figure
fig, ax = plt.subplots(figsize=(6,4))

# Want to plot both with and without decoherence
for include_decoherence in [ False, True ] :

    # If including decoherence, pass the model parameters to nuSQuIDS
    if include_decoherence :
        nusquids.Set_DecoherenceGammaEnergyDependence(n)
        nusquids.Set_DecoherenceGammaEnergyScale(E0_eV*units.eV)
        nusquids.Set_DecoherenceGammaMatrix(D_eV*units.eV)

    # Get transition probability, vs distance
    L_values_km = np.linspace(0., L_km, num=500)
    Pmm = [] # Placeholder to store numu survival prob vs distance
    for i, L in enumerate(L_values_km) :
        track = nsq.Vacuum.Track(L*units.km)
        nusquids.Set_Track(track)
        nusquids.Set_initial_state( np.array([0.,1.,0.]), nsq.Basis.flavor ) # Start with a pure numu flux
        nusquids.EvolveState()
        Pmm.append( nusquids.EvalFlavor(1) ) # Get final numu flux
    Pmm = np.asarray(Pmm)

    # Plot
    color = "orange" if include_decoherence else "grey"
    label = "+ Decoherence" if include_decoherence else "Standard oscillations"
    ax.plot(L_values_km, Pmm, lw=3, linestyle="-", color=color, label=label)

# Format
ax.grid(True)
ax.legend()
ax.set_xlim(L_values_km[0], L_values_km[-1])
ax.set_ylim(-0.02, 1.02)
ax.set_xlabel(r"$L$ [km]")
ax.set_ylabel(r"$P(\nu_\mu \rightarrow \nu_\mu)$")
fig.tight_layout()

# Save fig
fig_path = "decoherence.png"
fig.savefig(fig_path)
print(f"\nSaved : {fig_path}")





# nusquids.Set_ProgressBar(False)


# #nusquids.Set_CPPhase( 0, 2, np.pi)


# #print "nuSQuIDS PMNS:"
# #nusquids.PrintTransformationMatrix()

# # Function for calculating nuSQuIDS osc probs 
# def get_nusq_osc_probs() :
#     prob = np.full( (len(L_km),len(flavs)), np.NaN )
#     for i_L,L in enumerate(L_km) :
# #        track = nsq.Vacuum.Track(L*units.km) if V is None else nsq.ConstantDensity.Track(L*units.km)
#         track = nsq.Vacuum.Track(L*units.km) if V is None else nsq.VariableDensity.Track(L*units.km)
#         nusquids.Set_Track(track)
#         nusquids.Set_initial_state( np.array([0.,1.,0.]), nsq.Basis.flavor ) # Start with a pure numu flux
#         nusquids.EvolveState()
#             prob[i_L,i_ff] = nusquids.EvalFlavor( flav, E_GeV*units.GeV, 0 ) if use_energy_nodes else nusquids.EvalFlavor( flav )
#     return prob

# if include_std_osc :
#     nusquids.Set_DecoherenceGammaMatrix( 0., 0., 0. )
#     std_osc_probs_nsq = get_nusq_osc_probs()

# if include_decoherence : 
#     nusquids.Set_DecoherenceGammaMatrix( gamma_eV[0], gamma_eV[1], gamma_eV[2] )
#     decoh_osc_probs_nsq = get_nusq_osc_probs()
