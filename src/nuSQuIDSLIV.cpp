#include <math.h>
#include <nuSQuIDS/nuSQuIDSLIV.h>

namespace nusquids{

void nuSQUIDSLIV::AddToPreDerive(double x) {

    for(unsigned int ei = 0; ei < ne; ei++){

        // asumming same mass hamiltonian for neutrinos/antineutrinos
        squids::SU_vector h0 = H0( E_range[ei], 0 );

        LIVP_Eindep_evol[ei] = LIVP_Eindep.Evolve( h0, (x-Get_t_initial()) );
        LIVP_Edep_evol[ei] = LIVP_Edep.Evolve( h0, (x-Get_t_initial()) );

    }

}


squids::SU_vector nuSQUIDSLIV::HI(unsigned int ie, unsigned int irho) const {

    // Get the standard time-dependent Hamilonian
    squids::SU_vector potential = nuSQUIDS::HI(ie, irho);

    // Add SME effective Hamiltonian term
    double sign = 1;
    // if (NT==antineutrino){                                      // maybe add (irho == 1 and NT==both)   //TODO Implement this?
    //     // antineutrino matter potential flips sign
    //     sign*=(-1);
    potential += sign * ( E_range[ie] * LIVP_Edep_evol[ie] + LIVP_Eindep_evol[ie]);

    return potential;
}


void nuSQUIDSLIV::init() {

    // Allocate some matrices
    LIVP_Eindep_evol.resize(ne);
    LIVP_Edep_evol.resize(ne);

    for(unsigned int ei = 0; ei < ne; ei++){
        LIVP_Eindep_evol[ei] = squids::SU_vector(nsun);
        LIVP_Edep_evol[ei] = squids::SU_vector(nsun);
    }

}


void nuSQUIDSLIV::Set_LIVCoefficient(const marray<double,3>& a_mat,const marray<double,3>& c_mat,  double ra_rad, double dec_rad){

    // 
    // Assign LIV parameters:
    //

    // dmat is a 3D array of shape (3,3,3) containing a 3x3 matrix for each direction (x,y,z)
    // loop over directions and assign matrices to gsl matrices

    gsl_matrix_complex* a_eV_x = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* a_eV_y = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* a_eV_z = gsl_matrix_complex_calloc(3, 3);

    gsl_matrix_complex* c_tx = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* c_ty = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* c_tz = gsl_matrix_complex_calloc(3, 3);


    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            gsl_matrix_complex_set(a_eV_x, i, j, gsl_complex_rect(a_mat[0][i][j], 0));
            gsl_matrix_complex_set(a_eV_y, i, j, gsl_complex_rect(a_mat[1][i][j], 0));
            gsl_matrix_complex_set(a_eV_z, i, j, gsl_complex_rect(a_mat[2][i][j], 0));
            gsl_matrix_complex_set(c_tx, i, j, gsl_complex_rect(c_mat[0][i][j], 0));
            gsl_matrix_complex_set(c_ty, i, j, gsl_complex_rect(c_mat[1][i][j], 0));
            gsl_matrix_complex_set(c_tz, i, j, gsl_complex_rect(c_mat[2][i][j], 0));
        }
    }


    //
    // Define Coordinate system
    //

    // celestial colatitude and longitude
    double theta = M_PI/2 + dec_rad;
    double phi = ra_rad;

    // r vector
    double NX = sin(theta) * cos(phi);
    double NY = sin(theta) * sin(phi);
    double NZ = - cos(theta);



    //
    // Calculate LIV effective Hamiltonian with GSL matrix operations
    //

    // Amplitude equations: (cos(omega_sid L) amplitudes)
    // Ac0 = -NX * ax - NY * ay; 
    // Ac1 = 2 * NX * cxt + 2 * NY * cyt;
    // Const = NZ * az;

    // Declare Ac0, Ac1, and Const matrices
    gsl_matrix_complex* Ac0 = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* Ac1 = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* Const = gsl_matrix_complex_calloc(3, 3);

    // Calculate Ac0, Ac1, and Const matrices
    gsl_matrix_complex_scale(a_eV_x, gsl_complex_rect(-NX, 0));
    gsl_matrix_complex_scale(a_eV_y, gsl_complex_rect(-NY, 0));
    gsl_matrix_complex_memcpy(Ac0, a_eV_x);
    gsl_matrix_complex_add(Ac0, a_eV_y);

    gsl_matrix_complex_scale(c_tx, gsl_complex_rect(2 * NX, 0));
    gsl_matrix_complex_scale(c_ty, gsl_complex_rect(2 * NY, 0));
    gsl_matrix_complex_memcpy(Ac1, c_tx);
    gsl_matrix_complex_add(Ac1, c_ty);

    gsl_matrix_complex_memcpy(Const, a_eV_z);
    gsl_matrix_complex_scale(Const, gsl_complex_rect(NZ, 0));

    // Declare LIVP_Edep and LIVP_Eindep matrices
    gsl_matrix_complex* LIVP_Edep_GSL = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* LIVP_Eindep_GSL = gsl_matrix_complex_calloc(3, 3);

    // Calculate LIVP_Edep matrix
    gsl_matrix_complex_memcpy(LIVP_Edep_GSL, Ac1);

    // Calculate LIVP_Eindep matrix
    gsl_matrix_complex_memcpy(LIVP_Eindep_GSL, Ac0);
    gsl_matrix_complex_add(LIVP_Eindep_GSL, Const);

    // Heff =  Ac0 + Const + E * Ac1 = LIVP_Eindep + E * LIVP_Edep (see HI function above)
    LIVP_Edep = squids::SU_vector(LIVP_Edep_GSL);           // E dependent part of LIVP 
    LIVP_Eindep = squids::SU_vector(LIVP_Eindep_GSL);   // E independent part of LIVP


    //
    // Done
    //

    // free allocated matrix
    gsl_matrix_complex_free(a_eV_x);
    gsl_matrix_complex_free(a_eV_y);
    gsl_matrix_complex_free(a_eV_z);
    gsl_matrix_complex_free(c_tx);
    gsl_matrix_complex_free(c_ty);
    gsl_matrix_complex_free(c_tz);
    gsl_matrix_complex_free(Ac0);
    gsl_matrix_complex_free(Ac1);
    gsl_matrix_complex_free(Const);
    gsl_matrix_complex_free(LIVP_Edep_GSL);
    gsl_matrix_complex_free(LIVP_Eindep_GSL);

}

} // close namespace
