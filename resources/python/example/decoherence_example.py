'''
Simple example of using the neutrino decoherence implementation 

Tom Stuttard (thomas.stuttard@nbi.ku.dk)
'''

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import nuSQUIDSpy as nsq
from nuSQUIDSDecohPy import nuSQUIDSDecoh, nuSQUIDSDecohAtm

# Grab units
units = nsq.Const()

# Define flavours
flavs = [0, 1, 2]
flavs_tex = [r"\nu_{e}", r"\nu_{\mu}", r"\nu_{\tau}"]

# Create decoherence nusquids class instance
nusquids = nuSQUIDSDecoh(len(flavs), nsq.NeutrinoType.neutrino)

# Choose precision
error = 1.e-4
nusquids.Set_rel_error(error)
nusquids.Set_abs_error(error)

# Vacuum oscillations
nusquids.Set_Body(nsq.Vacuum())

# DeepCore-like values
E_GeV = 25.
nusquids.Set_E(E_GeV*units.GeV)
L_values_km = np.linspace(0., 12.7e3, num=100)

# Choose strength of decoherence
# 0 means standard oscillations
gamma_eV_values = [ 0., 1e-14 ]

# Create a figure for plotting
fig, ax = plt.subplots(3)

# Loop over cases
for gamma_eV, color, label in zip(gamma_eV_values, ["blue","red"], ["Standard", "Decoherence"]) :

    # Array to fill
    transition_prob = [ [] for flav in flavs ]

    # This is the matrix in equation 11 in arXiv:2007.00068
    nusquids.Set_DecoherenceGammaMatrixDiagonal(np.array([ 0., gamma_eV, gamma_eV, gamma_eV, gamma_eV, gamma_eV, gamma_eV, gamma_eV, gamma_eV ]))
    print("D matrix:")
    print(nusquids.Get_DecoherenceGammaMatrix())

    # Loop over ezch baseline to test
    for L_km in L_values_km :

        print(gamma_eV, L_km)

        # Set teh (vacuum) track
        nusquids.Set_Track( nsq.Vacuum.Track(L_km*units.km) )

        # Initially pure numu "beam"
        nusquids.Set_initial_state(np.array([0.,1.,0.]), nsq.Basis.flavor)

        # Evolve
        nusquids.EvolveState()

        # Get transition prob
        for i, flav in enumerate(flavs) :
            transition_prob[i].append( nusquids.EvalFlavor(flav) )

    # numpy-ify
    transition_prob = np.array(transition_prob)

    # Plot flavour transitions
    for i, flav in enumerate(flavs) :
        ax[i].plot(L_values_km, transition_prob[i], color=color, label=label)

# Format plot
for i, tex in enumerate(flavs_tex) :
    ax[i].set_ylabel( r"$P( \nu_{\mu} \rightarrow %s)$"%tex )
    ax[i].set_ylim(0., 1.)
    ax[i].set_xlim(L_values_km[0], L_values_km[-1])
    ax[i].grid(True)
    ax[i].legend()
ax[-1].set_ylabel( r"$E_{\nu}$ [GeV]" )

# Save fig
file_path = "decoherence_example.png"
fig.savefig(file_path)
print("Done! See %s" % file_path)

