#include <math.h>
#include <nuSQuIDS/nuSQuIDSBSMMatter.h>


namespace nusquids{


// Common initialisation tasks
void nuSQUIDSBSMMatter::init() {

  // Currently 3 flavors (could extend this if needed, but not supported yet)
  if ( numneu != 3 )
    std::cout << "WARNING : Code has only been tested for 3 neutrino flavors" << std::endl; 

  // Initialise NSI matrices
  vector_nsi_matrix_sun = squids::SU_vector(nsun);
  vector_nsi_matrix_sun.SetAllComponents(0.);

  scalar_nsi_matrix_sun = squids::SU_vector(nsun);
  scalar_nsi_matrix_sun.SetAllComponents(0.);

  // Init evolution vectors
  // Not current required as disabling interaction picture evolution
  // vector_nsi_matrix_sun_evol.resize(ne);
  // scalar_nsi_matrix_sun_evol.resize(ne);
  // for(int ei = 0; ei < ne; ei++){
  //   vector_nsi_matrix_sun_evol[ei] = squids::SU_vector(nsun);
  //   scalar_nsi_matrix_sun_evol[ei] = squids::SU_vector(nsun);
  // }

  // Pre-calculate factors
  V_prefactor = params.sqrt2 * params.GF * params.Na * pow(params.cm,-3);

}


void nuSQUIDSBSMMatter::Set_VectorNSIMatrix(const marray<double,2>& dmat) {

  // Check dimensions
  // Should be N x N
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
      if( dmat.extent(dim) != numneu)
        throw std::runtime_error("nuSQUIDSBSMMatter::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(numneu) + ")" );
  }

  // Create a GSL complex matrix
  gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
      for( unsigned int j = 0 ; j < numneu ; ++j ) {

          // Create the complex number to write
          gsl_complex c{{ dmat[i][j] , 0.0 }}; //Only using real part right now //TODO

          // Write to the matrix element
          gsl_matrix_complex_set(M, i, j, c);

      }
  }

  // Convert to SU(N) vector represenation
  // vector_nsi_matrix_sun.reset(gsl_matrix_complex);
  vector_nsi_matrix_sun = squids::SU_vector(M);

  // Rotate to mass basis
  vector_nsi_matrix_sun.RotateToB1(params);

  // Free memory from the temporary GSL matrix
  gsl_matrix_complex_free(M);

}



void nuSQUIDSBSMMatter::Set_ScalarNSIMatrix(const marray<double,2>& dmat) {

  // Check dimensions
  // Should be N x N
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
      if( dmat.extent(dim) != numneu)
        throw std::runtime_error("nuSQUIDSBSMMatter::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(numneu) + ")" );
  }

  // Create a GSL complex matrix
  gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
      for( unsigned int j = 0 ; j < numneu ; ++j ) {

          // Create the complex number to write
          gsl_complex c{{ dmat[i][j] , 0.0 }}; //Only using real part right now //TODO

          // Write to the matrix element
          gsl_matrix_complex_set(M, i, j, c);

      }
  }

  // Convert to SU(N) vector represenation
  scalar_nsi_matrix_sun = squids::SU_vector(M);

  // Rotate to mass basis
  scalar_nsi_matrix_sun.RotateToB1(params);

  // Free memory from the temporary GSL matrix
  gsl_matrix_complex_free(M);

}


marray<double,2> nuSQUIDSBSMMatter::Get_VectorNSIMatrix() const {

    auto gsl_matrix = vector_nsi_matrix_sun.GetGSLMatrix();

    marray<double,2> dmat{numneu, numneu};
    for( unsigned int i = 0 ; i < numneu ; ++i ) {
        for( unsigned int j = 0 ; j < numneu ; ++j ) {
          dmat[i][j] = GSL_REAL(gsl_matrix_complex_get(gsl_matrix.get(),i,j)); //TODO imaginary components
        }
    }

    return dmat;
}


marray<double,2> nuSQUIDSBSMMatter::Get_ScalarNSIMatrix() const {

    auto gsl_matrix = scalar_nsi_matrix_sun.GetGSLMatrix();

    marray<double,2> dmat{numneu, numneu};
    for( unsigned int i = 0 ; i < numneu ; ++i ) {
        for( unsigned int j = 0 ; j < numneu ; ++j ) {
          dmat[i][j] = GSL_REAL(gsl_matrix_complex_get(gsl_matrix.get(),i,j)); //TODO imaginary components
        }
    }

    return dmat;
}


void nuSQUIDSBSMMatter::iniH0(){

  // Disable H0_array, as I am handling H0 in HI
  // if(ienergy){ //TODO private -> protected in nuSQUIDS
  for(unsigned int ei = 0; ei < ne; ei++){
    // H0_array[ei] = H0(E_range[ei],0);
    auto null = squids::SU_vector(nsun);
    null.SetAllComponents(0.);
    H0_array[ei] = null;
  }
  // }

}


squids::SU_vector nuSQUIDSBSMMatter::H0(double Enu, unsigned int irho) const{

  // Disable H0, as I am handling H0 in HI

  // Define conventional mass matrix
  auto H0 = squids::SU_vector(nsun);
  H0.SetAllComponents(0.);
  return H0;

}


void nuSQUIDSBSMMatter::AddToPreDerive(double x){

  // If using interaction picture, need to evolve matrices
  // Howeber, since I am manually including H0 in HI I can forget this.
  // if (basis == interaction) { 
  //   for(int ei = 0; ei < ne; ei++){
  //     // asumming same hamiltonian for neutrinos/antineutrinos
  //     auto dt = x - Get_t_initial();
  //     vector_nsi_matrix_sun_evol[ei] = vector_nsi_matrix_sun.Evolve(H0_array[ei], dt);
  //     scalar_nsi_matrix_sun_evol[ei] = scalar_nsi_matrix_sun.Evolve(H0_array[ei], dt);
  //   }
  // }
}


squids::SU_vector nuSQUIDSBSMMatter::HI(unsigned int ei, unsigned int index_rho) const{

  // Note: Have included HI in here as well as HI as a hacky was of solving in the mass basis, not interaction basis
  // This is solving directly in the mass basis is broken in nuSQUIDSAtm (only Atm, looks to be in the dedicated 
  // EvalFlavor method, which always Evolves, unlike non-atmo case).
  // The reason I need to sole in mass basis is that H0 is no longer time indepdent in the case of scalar NSI (TODO depdends on how M and dM are combined, maybe can revert)

  // Considering coupling to quarks
  double num_quarks_per_electron = 3.;


  //
  // Mass matrix
  //

  // Define conventional mass matrix
  auto M_mass = squids::SU_vector(nsun);
  M_mass.SetAllComponents(0.);
  for(unsigned int i = 1; i < nsun; i++){
      M_mass += (b0_proj[i])*params.GetEnergyDifference(i);
  }

  // Construct scalar NSI delta M matrix in flavor basis
  double M_prefactor = params.GetEnergyDifference(2) * current_density * current_ye * num_quarks_per_electron;
  // squids::SU_vector M_bsm( M_prefactor * scalar_nsi_matrix_sun_evol[ei] );
  squids::SU_vector M_bsm( M_prefactor * scalar_nsi_matrix_sun );

  // Get overall mass matrix
  auto M_tot = squids::SU_vector( ( M_mass + M_bsm ) * ( 0.5 / E_range[ei] ) ); //TODO combine matrices proeprly



  //
  // Matter potential
  //

  // Get std matter potential
  // Note that H0_array is taken into accout here too
  auto V_std = nuSQUIDS::HI(ei,index_rho); //TODO anything else in there?

  // Construct vector NSI potential in flavor basis
  double CC = V_prefactor * current_density * current_ye * num_quarks_per_electron;
  // squids::SU_vector V_bsm( CC * vector_nsi_matrix_sun_evol[ei] );
  squids::SU_vector V_bsm( CC * vector_nsi_matrix_sun );

  // Negative for antineutrinos
  if ((index_rho == 1 and NT==both) or NT==antineutrino){
      V_bsm = (-1.0) * std::move(V_bsm);
  }

  // Get overall potential
  auto V_tot = squids::SU_vector(V_std + V_bsm); //TODO


  //
  // Hamiltonian
  //

  // Combine to get full Hamiltonian
  auto H = M_tot + V_tot;

  return H;

}


} // close namespace
