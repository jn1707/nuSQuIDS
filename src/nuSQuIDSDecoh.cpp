/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/

#include <nuSQuIDS/nuSQuIDSDecoh.h>


//TODO Put this somewhere, proper const correctness, etc
void print_gsl_matrix(gsl_matrix_complex* matrix) {
  for (size_t i = 0; i < matrix->size1; i++) {
    for (size_t j = 0; j < matrix->size2; j++) {
      std::cout << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
    }
    std::cout << std::endl;
  }
}


namespace nusquids{


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // TODO REMOVE???
  // For now want to calculate in the mass basis as is simpler
  // WARNING: There is no splined interpolation between energy nodes (maybe also coszen?) in the mass basis
  Set_Basis(mass);

  // Currently 3 flavors (could extend this if needed, but not supported yet)
  if ( numneu != 3 )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Only currently supporting 3 neutrino flavors");

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
      std::cout << "Setting [" << i << "," << j << "] = " << GSL_REAL(c) << std::endl;
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

#if 0
  // Update decoherence matrix // TODO Is ther a more efficient way to set the elements of a SU_vector?
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      decoherence_gamma_matrix[k] = dmat[i][j];
    }
  }
#endif

  //if(print) std::cout << "+++ decoherence_gamma_matrix = " << std::endl; //TODO REMOVE


}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(double Gamma21,double Gamma31,double Gamma32) {//, Basis basis = flavor) {  //TODO specify basis
  //TODO Is there a more efficient way to do this?
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


squids::SU_vector nuSQUIDSDecoh::D_Rho(unsigned int ei,unsigned int index_rho, double t) const {

  // This function returns the D[rho] dissipation/decoheence term
  // It gets used like: drho/dt = -i[H,rho] - D[rho]
  // D[rho] is a NxN matrix (N is num neutrino flavors) where each element is: D[rho]_ij = rho_i,j * Gamma_i,j

  // Not currently supporting calculating in the flavor basis TODO Move this check somewhere that is not called constantly...
  if(basis == flavor)
    throw std::runtime_error("nuSQUIDSDecoh::Error::Flavor basis is not supported fo the calculation");


  // Define shifts between mass and interaction basis
  double int_to_mass_basis_time_shift = Get_t_initial() - t; // Backwards in time from t to t0
  double mass_to_int_basis_time_shift = t - Get_t_initial(); // Forwards in time from t0 to t
//  double int_to_mass_basis_time_shift = t - Get_t_initial();
//  double mass_to_int_basis_time_shift = Get_t_initial() - t;

//TODO need to decide how to do this: This version rotates rho to the mass basis for the calculation, the D[rho] back to the interaction basis
#if 0
  std::vector<double> D_rho_buffer_mass_basis(nsun*nsun); //TODO Member

  // Get the Gamma matrix elements (mass basis)
  auto gamma_mass_basis = decoherence_gamma_matrix.GetComponents(); //TODO Is this a waste of a memory allocation?

  // Interaction basis case
  // TODO Make smarter toggles with less code reproduction here...
  if(basis == interaction) {

    // Get rho and rotate from interaction to mass basis
    //squids::SU_vector rho_evol = estate[ei].rho[index_rho].Evolve(H0_array[ei],int_to_mass_basis_time_shift);
    //auto rho_mass_basis = rho_evol.GetComponents();
    auto rho_mass_basis = estate[ei].rho[index_rho];
/*
    if(false) {
      std::cout << std::endl << "t = " << t << std::endl;
      std::cout << "rho:" << std::endl;
      print_gsl_matrix(rho_evol.GetGSLMatrix().get());
    }
*/
    // Multiply the Gamma and rho matrices element-wise to get D[rho]
    //TODO Does this SU vector element-wise multiplication do the same thing multiplying the elements in the standard 3x3 GSL matrix?
    for( unsigned int i = 0 ; i < (nsun*nsun) ; ++i ) {
        D_rho_buffer_mass_basis[i] = rho_mass_basis[i] * gamma_mass_basis[i];
    }
    squids::SU_vector D_rho_mass_basis(D_rho_buffer_mass_basis);
    return D_rho_mass_basis;

    // Rotate D[rho] to from the mass basis to the interaction basis //TODO Is this definitely mathemetically valid?
    squids::SU_vector D_rho_int_basis = D_rho_mass_basis.Evolve(H0_array[ei],mass_to_int_basis_time_shift);

    return D_rho_int_basis;

  }
  else if(basis == mass) {

    std::cout << "Mass basis" << std::endl;

    auto rho_mass_basis = estate[ei].rho[index_rho].GetComponents();

    for( unsigned int i = 0 ; i < (nsun*nsun) ; ++i ) {
        D_rho_buffer_mass_basis[i] = rho_mass_basis[i] * gamma_mass_basis[i];
    }
    squids::SU_vector D_rho_mass_basis(D_rho_buffer_mass_basis);
    return D_rho_mass_basis;

  }
#endif



//TODO I don't think this should work as I'm not handling the bases at all here, but somehow it seems to exactly match my solver. What is going on!?!
//TODO Perhaps my gamma matrix is automatically converted to the correct basis when I create the SU_vector? 
#if 1
  std::vector<double> D_rho_buffer_tmp(nsun*nsun); //TODO Member

  // Get the Gamma matrix elements (mass basis)
  auto gamma = decoherence_gamma_matrix.GetComponents(); //TODO Is this a waste of a memory allocation?

  auto rho = estate[ei].rho[index_rho].GetComponents();

  for( unsigned int i = 0 ; i < (nsun*nsun) ; ++i ) {
      D_rho_buffer_tmp[i] = rho[i] * gamma[i];
  }
  squids::SU_vector D_rho_val(D_rho_buffer_tmp);
  return D_rho_val;

#endif




//TODO need to decide how to do this: This version rotates gamma to the interaction basis for the calculation
#if 0
  std::vector<double> D_rho_buffer_int_basis(nsun*nsun); //TODO Member

  // Get rho in the interactio  basis (it already is)
  auto rho_int_basis = estate[ei].rho[index_rho];

  // Rotate Gamma matrix to interaction basis
  //TODO Do this in `PreDerive` instead for efficiency
  squids::SU_vector gamma_evol = decoherence_gamma_matrix.Evolve(H0_array[ei],mass_to_int_basis_time_shift);
  auto gamma_int_basis = gamma_evol.GetComponents();

  // Multiply the Gamma and rho matrices element-wise to get D[rho]
  //TODO Does this SU vector element-wise multiplication do the same thing multiplying the elements in the standard 3x3 GSL matrix?
  for( unsigned int i = 0 ; i < (nsun*nsun) ; ++i ) {
      D_rho_buffer_int_basis[i] = rho_int_basis[i] * gamma_int_basis[i];
  }
  squids::SU_vector D_rho_int_basis(D_rho_buffer_int_basis);
#endif



// This is how to do the rho * Gamma elementwise calculation using GSL amtrix instead of SU_vector //TODO Remove this once sure don't need it
#if 0
  //TODO state or estate?
  // Convert rho to mass basis (from interaction basis) //TODO PreDerive, e.g. only once
  //squids::SU_vector state_reversal = estate[ei].rho[index_rho].Evolve(H0_array[ei],first_shift);
  squids::SU_vector& state_reversal = estate[ei].rho[index_rho]; //TODO REMOVE

  auto rho_matrix_mass_basis = state_reversal.GetGSLMatrix();
  if(print) std::cout << "rho :" << std::endl;
  if(print) print_gsl_matrix(rho_matrix_mass_basis.get());

  auto D_rho_matrix_mass_basis = state_reversal.GetGSLMatrix();

  //TODO Use a GSL buffer D_Rho
  for( unsigned int i = 0 ; i < numneu ; ++i ) { //TODO Do this in a more efficient way
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      gsl_complex D_rho_ij{{ GSL_REAL(gsl_matrix_complex_get(decoherence_gamma_gsl_matrix.get(),i,j)) * GSL_REAL(gsl_matrix_complex_get(rho_matrix_mass_basis.get(),i,j)), 0. }};
      gsl_matrix_complex_set(D_rho_matrix_mass_basis.get(),i,j,D_rho_ij);
    }
  }
  if(print) std::cout << "D[rho] :" << std::endl;
  if(print) print_gsl_matrix(D_rho_matrix_mass_basis.get());

  squids::SU_vector D_rho(D_rho_matrix_mass_basis.get());
#endif


  //return D_rho_int_basis;
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
  print_gsl_matrix(rho_gsl.get());
}




} // close namespace


