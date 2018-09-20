/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/
#include <math.h>
#include <nuSQuIDS/nuSQuIDSDecoh.h>

//TODO Put this somewhere, proper const correctness, etc
#if 1
#include <iostream>
#include <iomanip>
void print_gsl_matrix(gsl_matrix_complex* matrix) {
  for (size_t i = 0; i < matrix->size1; i++) {
    for (size_t j = 0; j < matrix->size2; j++) {
      //std::cout << std::fixed << std::setprecision(2) << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
      std::cout << std::scientific << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
    }
    std::cout << std::endl;
  }
}
#endif

namespace nusquids{


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // For now want to calculate in the mass basis as is simpler
  // WARNING: There is no splined interpolation between energy nodes (maybe also coszen?) in the mass basis
  //Set_Basis(mass); //TODO Seems to be a problem here, results don't match

  // Currently 3 flavors (could extend this if needed, but not supported yet)
  if ( ! ( (numneu == 2) ||  (numneu == 3) ) )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Only currently supporting 2 or 3 neutrino flavors");

  // Enable decoherence term in the nuSQuIDS numerical solver
  Set_DecoherenceTerms(true);

  // Initialise the decoherence matrix //TODO Is there a simpler way to do this? If stick with SU_vector can use the SetComponentsToZero function
  marray<double,2> tmp_dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      tmp_dmat[i][j] = 0.;
    }
  }
  Set_DecoherenceGammaMatrix(tmp_dmat);
  
  // Init interaction basis gamma matrices
#if 0
  decoherence_gamma_matrix_evol.resize(ne);
  for(int ei = 0; ei < ne; ei++){
    decoherence_gamma_matrix_evol[ei] = squids::SU_vector(nsun);
  }
#endif

}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(const marray<double,2>& dmat) {//, Basis basis = flavor) {  //TODO specify basis

  // Check dimensions (should be NxN, where N is number neutrino flavors)
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
    if( dmat.extent(dim) != numneu)
      throw std::runtime_error("nuSQUIDSDecoh::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(numneu) + ")" );
  }

  //TODO Do I need to use this GSL stuff here? I took this from the NSI example, ask Carlos about this...
  // Defining a complex matrix M to contain our decoherence parameterisation // Does it need to be complex?
  gsl_matrix_complex * M = gsl_matrix_complex_calloc(numneu,numneu);

  // Loop over matrix elements (2D)
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {

      // Write this element to the matrix
      gsl_complex c {{ dmat[i][j] , 0.0 }}; //Only using real part right now
      gsl_matrix_complex_set(M,i,j,c); // TODO there is a gsl_complex_conjugate method, could be useful...

    }
  }

  decoherence_gamma_gsl_matrix.reset(M);

  // Update the SU(N) vector from the matrix   
  decoherence_gamma_matrix = squids::SU_vector(M);

  // TODO rotate to basis
  //decoherence_gamma_matrix.RotateToB1(params);
  
  // Done with the matrix now
  //gsl_matrix_complex_free(M); //TODO

}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(double Gamma) {

  if ( numneu != 2 )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Function only valid for 2 neutrino flavors");

  marray<double,2> dmat{2,2};
  dmat[0][0] = 0.;
  dmat[0][1] = Gamma;
  dmat[1][0] = Gamma;
  dmat[1][1] = 0.;
  Set_DecoherenceGammaMatrix(dmat);
}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(double Gamma21,double Gamma31,double Gamma32) {//, Basis basis = flavor) {  //TODO specify basis

  if ( numneu != 3 )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Function only valid for 3 neutrino flavors");

  //TODO Is there a more efficient way to do this (is part of the PISA minimizer so is called multiple times)?
  marray<double,2> dmat{3,3};
  dmat[0][0] = 0.;
  dmat[0][1] = Gamma21;
  dmat[0][2] = Gamma31;
  dmat[1][0] = Gamma21;
  dmat[1][1] = 0.;
  dmat[1][2] = Gamma32;
  dmat[2][0] = Gamma31;
  dmat[2][1] = Gamma32;
  dmat[2][2] = 0.;

  Set_DecoherenceGammaMatrix(dmat);

}


// Get the current value of the decoherence matrix
// This is not used as part of the main solver, just for user scripts
marray<double,2> nuSQUIDSDecoh::Get_DecoherenceGammaMatrix() const {
  //auto decoherence_gamma_gsl_matrix = decoherence_gamma_matrix.GetGSLMatrix(); //TODO use this instead, including imaginary
  marray<double,2> dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) { //TODO Do this in a more efficient way
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      dmat[i][j] = decoherence_gamma_matrix[k];
    }
  }
  return dmat;
}



void nuSQUIDSDecoh::Set_EnergyDependence(double n)  {
  n_energy = n; 
}

double nuSQUIDSDecoh::Get_EnergyDependence() const {
  return n_energy; 
}


squids::SU_vector nuSQUIDSDecoh::D_Rho(unsigned int ei,unsigned int index_rho, double t) const {

  // This function returns the D[rho] dissipation/decoheence term
  // It gets used like: drho/dt = -i[H,rho] - D[rho]
  // D[rho] is a NxN matrix (N is num neutrino flavors) where each element is: D[rho]_ij = rho_i,j * Gamma_i,j

  //
  // Prepare for calculation
  //

  // Not currently supporting calculating in the flavor basis TODO Move this check somewhere that is not called constantly...
  if(basis == flavor)
    throw std::runtime_error("nuSQUIDSDecoh::Error::Flavor basis is not supported for the calculation");

  if(basis == mass)
    throw std::runtime_error("nuSQUIDSDecoh::Error::Mass basis is not supported for the calculation");

  // Define shifts between mass and interaction basis
  double int_to_mass_basis_time_shift = Get_t_initial() - t; // Backwards in time from t to t0
  double mass_to_int_basis_time_shift = t - Get_t_initial(); // Forwards in time from t0 to t

  // Get the energy dependence
  double n_energy = Get_EnergyDependence();

  // Get energy conversion
  double energy_conversion = pow(10,-9);


  //
  // Determine D[rho] from the gamma matrix and rho
  //

#define SOLVE_IN_MASS_BASIS 0
#define SOLVE_IN_INT_BASIS 0

  // Get the components of rho
  auto rho_int_basis = estate[ei].rho[index_rho];
#if SOLVE_IN_MASS_BASIS
  // Convert rho from int to mass basis
  squids::SU_vector rho_mass_basis = estate[ei].rho[index_rho].Evolve(H0_array[ei],int_to_mass_basis_time_shift);
  auto rho = rho_mass_basis.GetComponents();
#else
  auto rho = rho_int_basis.GetComponents();
#endif

  // Get the components of gamma
#if SOLVE_IN_INT_BASIS
  // Convert gamma from mass to int basis
  squids::SU_vector gamma_int_basis = decoherence_gamma_matrix.Evolve(H0_array[ei],mass_to_int_basis_time_shift);
  auto gamma = gamma_int_basis.GetComponents();
#else
  auto gamma = decoherence_gamma_matrix.GetComponents();
#endif

  // Element-wise SU vector multiplication to get D[rho] //TODO replace with Chris' new function
  std::vector<double> D_rho_vect(nsun*nsun); //Cannot use a member variable as the buffer since the function is constant
  for( unsigned int i = 0 ; i < (nsun*nsun) ; ++i ) {
      D_rho_vect[i] = rho[i] * gamma[i] * pow( E_range[ei] * energy_conversion, n_energy);
  }

  // Convert to SU vector
#if SOLVE_IN_MASS_BASIS
  // Convert back to interaction basis
  squids::SU_vector D_rho_mass_basis(D_rho_vect);
  squids::SU_vector D_rho_val = D_rho_mass_basis.Evolve(H0_array[ei],mass_to_int_basis_time_shift);
#else
  squids::SU_vector D_rho_val(D_rho_vect);
#endif

  return D_rho_val;

}


#if 0
void nuSQUIDSDecoh::AddToPreDerive(double x){
  for(int ei = 0; ei < ne; ei++){
    // Get the gamma matrix in the interaction basis
    decoherence_gamma_matrix_evol[ei] = decoherence_gamma_matrix.Evolve(H0_array[ei],(x-Get_t_initial())); // mass basis -> interaction basis
  }
}
#endif



void nuSQUIDSDecoh::PrintTransformationMatrix() const {
  auto transformation_matrix = GetTransformationMatrix().get();
  for(unsigned int i = 0 ; i < nsun ; ++i) {
    for(unsigned int j = 0 ; j < nsun ; ++j) {
      std::cout << GSL_REAL(gsl_matrix_complex_get(transformation_matrix,i,j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(transformation_matrix,i,j)) << "i   "; 
    }
    std::cout << std::endl;
  }
}


void nuSQUIDSDecoh::PrintState() const {
  auto rho_gsl = estate[0].rho[0].GetGSLMatrix();
//  print_gsl_matrix(rho_gsl.get());
}


} // close namespace
