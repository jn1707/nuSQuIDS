#include <math.h>
#include <nuSQuIDS/nuSQuIDSLIV.h>

namespace nusquids {

void nuSQUIDSLIV::AddToPreDerive(double x) {
    for (unsigned int ei = 0; ei < ne; ei++) {
        squids::SU_vector h0 = H0(E_range[ei], 0);

        CPT_odd_Eindep_evol[ei] = CPT_odd_Eindep.Evolve(h0, (x - Get_t_initial()));
        CPT_odd_Edep_evol[ei] = CPT_odd_Edep.Evolve(h0, (x - Get_t_initial()));
        CPT_even_Eindep_evol[ei] = CPT_even_Eindep.Evolve(h0, (x - Get_t_initial()));
        CPT_even_Edep_evol[ei] = CPT_even_Edep.Evolve(h0, (x - Get_t_initial()));
    }
}

squids::SU_vector nuSQUIDSLIV::HI(unsigned int ie, unsigned int irho) const {
    squids::SU_vector potential = nuSQUIDS::HI(ie, irho);

    double sign = 1;
    if ((irho == 1 and NT == both) or NT == antineutrino) {
        sign *= -1;
    }

    // Add SME Hamiltonian term
    potential += sign * ( E_range[ie] * CPT_odd_Edep_evol[ie] + CPT_odd_Eindep_evol[ie] ); // CPT odd
    potential += (E_range[ie] * CPT_even_Edep_evol[ie] + CPT_even_Eindep_evol[ie]); // CPT even

    return potential;
}

void nuSQUIDSLIV::init() {
    CPT_odd_Eindep_evol.resize(ne);
    CPT_odd_Edep_evol.resize(ne);
    CPT_even_Eindep_evol.resize(ne);
    CPT_even_Edep_evol.resize(ne);

    for (unsigned int ei = 0; ei < ne; ei++) {
        CPT_odd_Eindep_evol[ei] = squids::SU_vector(nsun);
        CPT_odd_Edep_evol[ei] = squids::SU_vector(nsun);
        CPT_even_Eindep_evol[ei] = squids::SU_vector(nsun);
        CPT_even_Edep_evol[ei] = squids::SU_vector(nsun);
    }
}

void nuSQUIDSLIV::Set_LIVCoefficient(const marray<double, 3>& a_mat, const marray<double, 4>& c_mat, double ra_rad, double dec_rad) {

    // Init "a" matrix (N,N) per element
    gsl_matrix_complex* a_eV_t0 = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* a_eV_x = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* a_eV_y = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* a_eV_z = gsl_matrix_complex_calloc(nsun, nsun);

    // Init "a" matrix (N,N) per element
    gsl_matrix_complex* c_tt = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_tx = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_ty = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_tz = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_xx = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_xy = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_xz = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_yy = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_yz = gsl_matrix_complex_calloc(nsun, nsun);
    gsl_matrix_complex* c_zz = gsl_matrix_complex_calloc(nsun, nsun);

    // Set matrix elements based on the inputs
    for (size_t i = 0; i < nsun; i++) {
        for (size_t j = 0; j < nsun; j++) {

            gsl_matrix_complex_set(a_eV_t0, i, j, gsl_complex_rect(a_mat[0][i][j], 0)); // Note that setting imaginary part to zero here (not currently supported)
            gsl_matrix_complex_set(a_eV_x, i, j, gsl_complex_rect(a_mat[1][i][j], 0));
            gsl_matrix_complex_set(a_eV_y, i, j, gsl_complex_rect(a_mat[2][i][j], 0));
            gsl_matrix_complex_set(a_eV_z, i, j, gsl_complex_rect(a_mat[3][i][j], 0));

            gsl_matrix_complex_set(c_tt, i, j, gsl_complex_rect(c_mat[0][0][i][j], 0));
            gsl_matrix_complex_set(c_tx, i, j, gsl_complex_rect(c_mat[0][1][i][j], 0));
            gsl_matrix_complex_set(c_ty, i, j, gsl_complex_rect(c_mat[0][2][i][j], 0));
            gsl_matrix_complex_set(c_tz, i, j, gsl_complex_rect(c_mat[0][3][i][j], 0));
            gsl_matrix_complex_set(c_xx, i, j, gsl_complex_rect(c_mat[1][1][i][j], 0));
            gsl_matrix_complex_set(c_xy, i, j, gsl_complex_rect(c_mat[1][2][i][j], 0));
            gsl_matrix_complex_set(c_xz, i, j, gsl_complex_rect(c_mat[1][3][i][j], 0));
            gsl_matrix_complex_set(c_yy, i, j, gsl_complex_rect(c_mat[2][2][i][j], 0));
            gsl_matrix_complex_set(c_yz, i, j, gsl_complex_rect(c_mat[2][3][i][j], 0));
            gsl_matrix_complex_set(c_zz, i, j, gsl_complex_rect(c_mat[3][3][i][j], 0));

        }
    }
    // Print the matrices
    std::cout << "a_eV_t0: " << a_eV_t0 << std::endl;
    std::cout << "a_eV_x: " << a_eV_x << std::endl;
    std::cout << "a_eV_y: " << a_eV_y << std::endl;
    std::cout << "c_ty: " << c_ty << std::endl;

    // Convert celestial coordinates to spherical coordinates
    double theta = M_PI / 2 - dec_rad; // Colatitude
    double phi = ra_rad; // Longitude

    // Calculate the unit propagation vectors in spherical coordinates
    double p_x = -sin(theta) * cos(phi);
    double p_y = -sin(theta) * sin(phi);
    double p_z = -cos(theta);

    // Compute operators...
    gsl_matrix_complex* CPT_odd_Eindep_GSL = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* CPT_odd_Edep_GSL = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* CPT_even_Eindep_GSL = gsl_matrix_complex_calloc(3, 3);
    gsl_matrix_complex* CPT_even_Edep_GSL = gsl_matrix_complex_calloc(3, 3);

    // CPT-odd energy-independent part
    gsl_matrix_complex* a_term = gsl_matrix_complex_calloc(3, 3);

    gsl_matrix_complex_memcpy(a_term, a_eV_t0);
    gsl_matrix_complex_add(CPT_odd_Eindep_GSL, a_term);

    gsl_matrix_complex_memcpy(a_term, a_eV_x);
    gsl_matrix_complex_scale(a_term, gsl_complex_rect(-p_x, 0));
    gsl_matrix_complex_add(CPT_odd_Eindep_GSL, a_term);

    gsl_matrix_complex_memcpy(a_term, a_eV_y);
    gsl_matrix_complex_scale(a_term, gsl_complex_rect(-p_y, 0));
    gsl_matrix_complex_add(CPT_odd_Eindep_GSL, a_term);

    gsl_matrix_complex_memcpy(a_term, a_eV_z);
    gsl_matrix_complex_scale(a_term, gsl_complex_rect(-p_z, 0));
    gsl_matrix_complex_add(CPT_odd_Eindep_GSL, a_term);

    // CPT-even energy-dependent part
    gsl_matrix_complex* cterm = gsl_matrix_complex_calloc(3, 3);

    gsl_matrix_complex_memcpy(cterm, c_tt);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(0.5 * (3 - p_z * p_z), 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_zz);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(0.5 * (-1 + 3 * p_z * p_z), 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_tx);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(2 * p_x, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_ty);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(2 * p_y, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);
    
    gsl_matrix_complex_memcpy(cterm, c_tz);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-2 * p_z, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_xx);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-p_x * p_x, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_xy);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-2 * p_x * p_y, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_xz);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-2 * p_x * p_z, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_yy);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-p_y * p_y, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    gsl_matrix_complex_memcpy(cterm, c_yz);
    gsl_matrix_complex_scale(cterm, gsl_complex_rect(-2 * p_y * p_z, 0));
    gsl_matrix_complex_add(CPT_even_Edep_GSL, cterm);

    // Print the matrices
    std::cout << "CPT_odd_Eindep_GSL: " << CPT_odd_Eindep_GSL << std::endl;
    std::cout << "CPT_even_Edep_GSL: " << CPT_even_Edep_GSL << std::endl;
    CPT_even_Edep = squids::SU_vector(CPT_even_Edep_GSL);
    CPT_even_Eindep = squids::SU_vector(CPT_even_Eindep_GSL);
    CPT_odd_Eindep = squids::SU_vector(CPT_odd_Eindep_GSL);
    CPT_odd_Edep = squids::SU_vector(CPT_odd_Edep_GSL);

    gsl_matrix_complex_free(a_eV_t0);
    gsl_matrix_complex_free(a_eV_x);
    gsl_matrix_complex_free(a_eV_y);
    gsl_matrix_complex_free(a_eV_z);
    gsl_matrix_complex_free(c_tt);
    gsl_matrix_complex_free(c_tx);
    gsl_matrix_complex_free(c_ty);
    gsl_matrix_complex_free(c_tz);
    gsl_matrix_complex_free(c_xx);
    gsl_matrix_complex_free(c_xy);
    gsl_matrix_complex_free(c_xz);
    gsl_matrix_complex_free(c_yy);
    gsl_matrix_complex_free(c_yz);
    gsl_matrix_complex_free(c_zz);
    gsl_matrix_complex_free(a_term);
    gsl_matrix_complex_free(cterm);
    gsl_matrix_complex_free(CPT_odd_Eindep_GSL);
    gsl_matrix_complex_free(CPT_odd_Edep_GSL);
    gsl_matrix_complex_free(CPT_even_Eindep_GSL);
    gsl_matrix_complex_free(CPT_even_Edep_GSL);
}

}  // close namespace