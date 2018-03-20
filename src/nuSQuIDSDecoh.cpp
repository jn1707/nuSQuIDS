/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/

#include <nuSQuIDS/nuSQuIDSDecoh.h>

namespace nusquids{


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // Currently 3 flavors (could extend this if neeed, but not supported yet)
  if ( numneu != 3 )
    throw std::runtime_error("nuSQUIDSDecoh::Error::Only currently supporting 3 neutrino flavors");

  // Enable decoherence term in the nuSQuIDS numerical solver
  Set_DecoherenceTerms(true);

  // Initialise the decoherence matrix //TODO Is thre a simpelr way to do this?
  marray<double,2> tmp_dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      tmp_dmat[i][j] = 0.;
    }
  }
  Set_DecoherenceGammaMatrix(tmp_dmat);

}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(const marray<double,2>& dmat) {//, Basis basis = flavor) {  //TODO specify basis

  // Check dimensions (should be NxN, where N is number neutrino flavors)
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
    if( dmat.extent(dim) != numneu)
      throw std::runtime_error("nuSQUIDSDecoh::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(numneu) + ")" );
  }

  //TODO Do I need to use this GSL stuff here? I took this from the NSI example, ask Carlos about this...
#if 0
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

  // Update the SU(N) vector from the matrix   
  decoherence_gamma_matrix = squids::SU_vector(M);

  for(int i = 0 ; i< 9 ; i++ ) {
    std::cout << "Found in SU(3) [" << i << "] = " << decoherence_gamma_matrix[i] << std::endl;
  }

  // TODO SU_Vector seems to be missing data here????

  // Done with the matrix now
  gsl_matrix_complex_free(M);
#endif

#if 1
  // Update decoherence matrix // TODO Is ther a more efficient way to set the elements of a SU_vector?
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      decoherence_gamma_matrix[k] = dmat[i][j];
    }
  }
#endif

  // TODO rotate to basis
  //decoherence_gamma_matrix.RotateToB1(params);

}


void nuSQUIDSDecoh::Set_DecoherenceGammaMatrix(double Gamma21,double Gamma31,double Gamma32) {//, Basis basis = flavor) {  //TODO specify basis
  // TODO Instead call Set_DecoherenceGammaMatrix(array) so that anythign clever that function does with memory is used here
  decoherence_gamma_matrix[0] = 0.; //Is there a more efficient way to do this?
  decoherence_gamma_matrix[1] = Gamma21;
  decoherence_gamma_matrix[2] = Gamma31;
  decoherence_gamma_matrix[3] = Gamma21;
  decoherence_gamma_matrix[4] = 0.;
  decoherence_gamma_matrix[5] = Gamma32;
  decoherence_gamma_matrix[6] = Gamma31;
  decoherence_gamma_matrix[7] = Gamma32;
  decoherence_gamma_matrix[8] = 0.;
}


// Get the current value of the decoherence matrix
marray<double,2> nuSQUIDSDecoh::Get_DecoherenceGammaMatrix() const {
  marray<double,2> dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) { //TODO Do this in a more efficient way
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      dmat[i][j] = decoherence_gamma_matrix[k];
    }
  }
  return dmat;
}


squids::SU_vector nuSQUIDSDecoh::DecohGamma(unsigned int ei,unsigned int index_rho, double t) const {
  return decoherence_gamma_matrix;
}


} // close namespace


