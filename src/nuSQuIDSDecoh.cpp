/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/
#include <math.h>
#include <nuSQuIDS/nuSQuIDSDecoh.h>

//TODO Put this somewhere generic, with proper const correctness, etc
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

namespace nusquids{


void nuSQUIDSDecoh::EnableDecoherence(bool enable) {
  Set_DecoherenceTerms(enable);
}


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // Performing calculation in the interaction basi (which is the default)
  // This is the only basis for which splined interpolation between nodes is implemented
  // Could optionally try solving in e.g. the mass basis (Set_Basis(mass)), which would be simpler for debugging, but have discovered an issue with mass/falvor bassi calculation that needs resolving first

  // Currently 3 flavors (could extend this if needed, but not supported yet)
  if ( ! ( (numneu == 2) || (numneu == 3) ) )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Only currently supporting 2 or 3 neutrino flavors");

  // Enable decoherence term in the SQuIDS numerical solver
  EnableDecoherence(true);

  // Initialise the decoherence matrix (default is no decoherence)
  decoherence_gamma_matrix.SetAllComponents(0.);

  // Choose some defaults
  Set_DecoherenceBasis(mass);
  Set_DecoherenceGammaEnergyDependence(0);
  Set_DecoherenceGammaEnergyScale(1.*units.GeV);
  
}


void nuSQUIDSDecoh::Set_DecoherenceBasis(Basis basis) {
  decoherence_basis = basis;
}


Basis nuSQUIDSDecoh::Get_DecoherenceBasis() const {
  return decoherence_basis;
}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(const marray<double,2>& dmat) {

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
      gsl_complex c{{ dmat[i][j] , 0.0 }}; //Only using real part right now
      gsl_matrix_complex_set(M,i,j,c); // TODO there is a gsl_complex_conjugate method, could be useful...

    }
  }

  // Store the temporary matrix in the member variable
  // This takes over the memory for M, so no need to free(M)
  decoherence_gamma_gsl_matrix.reset(M);

  // Update the SU(N) vector from the matrix   
  decoherence_gamma_matrix = squids::SU_vector(M);

  // Handle mass vs flavor basis decoherence
  if( decoherence_basis == flavor ) {
      decoherence_gamma_matrix.RotateToB0(params);
  }

  
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

  //TODO Is there a more efficient way to do this? This is called at each step in the fit so is called many multiple times)?
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


// Get the current value of the decoherence gamma matrix
// Return the value as the marray type, which can be converted to a numpy array in pyhton
// This is not used as part of the main solver, just for user scripts
//TODO Currently this doesn't support imaginary components, but we don't support writing imaginary components right now either
marray<double,2> nuSQUIDSDecoh::Get_DecoherenceGammaMatrix() const {
  auto decoherence_gamma_gsl_matrix = decoherence_gamma_matrix.GetGSLMatrix();
  marray<double,2> dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      dmat[i][j] = GSL_REAL(gsl_matrix_complex_get(decoherence_gamma_gsl_matrix.get(),i,j));
    }
  }
  return dmat;
}


void nuSQUIDSDecoh::Set_DecoherenceGammaEnergyDependence(double n)  {
  gamma_energy_dep_index = n; 
}

double nuSQUIDSDecoh::Get_DecoherenceGammaEnergyDependence() const {
  return gamma_energy_dep_index; 
}

void nuSQUIDSDecoh::Set_DecoherenceGammaEnergyScale(double energy)  {
  gamma_energy_scale = energy; 
}

double nuSQUIDSDecoh::Get_DecoherenceGammaEnergyScale() const {
  return gamma_energy_scale; 
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

  // At time of writing there are known issues with the mass basis calculation
  if(basis == mass)
    throw std::runtime_error("nuSQUIDSDecoh::Error::Mass basis is not supported for the calculation");

  // Define shifts between mass and interaction basis
  //double int_to_mass_basis_time_shift = Get_t_initial() - t; // Backwards in time from t to t0
  //double mass_to_int_basis_time_shift = t - Get_t_initial(); // Forwards in time from t0 to t

  //
  // Determine D[rho] from the gamma matrix and rho
  //


  // Get the energy dependence of the Gamma damping parameters for this energy
  // Note that we define the Gamma value relative to some user-defined scale
  double energy_dependence = pow( (E_range[ei] / gamma_energy_scale) , gamma_energy_dep_index);

  // Perform element-wise SU vector multiplication to get D[rho] operator (as a NxN GSL matrix)
  // Include the energy dependence
  squids::SU_vector D_rho_val( squids::ElementwiseProduct( estate[ei].rho[index_rho], decoherence_gamma_matrix * energy_dependence ) );

  return D_rho_val;

}


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
 print_gsl_matrix(rho_gsl.get());
}


} // close namespace
