/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/

//TODO Create a .cpp file

#include <nuSQuIDS/nuSQuIDSDecoh.h>

namespace nusquids{


// Common initialisation tasks
void nuSQUIDSDecoh::init() {

  // Enable deocherence term in the nuSQuIDS numerical solver
  Set_DecoherenceTerms(true);

  // By default use an empty deocherene matrix
  decoherence_matrix.SetAllComponents(0.); 

}


// Set the decoherence matrix elements (set each individually)
void nuSQUIDSDecoh::Set_DecoherenceMatrix(const marray<double,2>& dmat) {//, Basis basis = flavor) {  //TODO specify basis?

  // Check dimensions (should be NxN, where N is number neutrino flavors)
  for( unsigned int dim = 0 ; dim < 2 ; ++dim ) {
    if( dmat.extent(dim) != numneu)
      throw std::runtime_error("nuSQUIDSDecoh::Error:Input decoherence matrix has wrong dimensions in dimension " + std::to_string(dim) + " (found " + std::to_string(dmat.extent(dim))+ ", should be " + std::to_string(numneu) + ")" );
  }

  // Update decoherence matrix /TDO Is ther a more efficient way to set the elements of a SU_vector?
  for( unsigned int i = 0 ; i < numneu ; ++i ) {
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      decoherence_matrix[k] = dmat[i][j];
    }
  }

}


// Set the decoherence matrix elements using the Gamma21,31,32 model
// See https://arxiv.org/pdf/1708.05495.pdf equation 12
void nuSQUIDSDecoh::Set_DecoherenceMatrix(double Gamma21,double Gamma31,double Gamma32) {//, Basis basis = flavor) {  //TODO specify basis?
  decoherence_matrix[0] = 0.; //Is there a more efficient way to do this?
  decoherence_matrix[1] = Gamma21;
  decoherence_matrix[2] = Gamma31;
  decoherence_matrix[3] = Gamma21;
  decoherence_matrix[4] = 0.;
  decoherence_matrix[5] = Gamma32;
  decoherence_matrix[6] = Gamma31;
  decoherence_matrix[7] = Gamma32;
  decoherence_matrix[8] = 0.;
}


// Get the current value of the decoherence matrix
marray<double,2> nuSQUIDSDecoh::Get_DecoherenceMatrix() const {
  marray<double,2> dmat{numneu,numneu};
  for( unsigned int i = 0 ; i < numneu ; ++i ) { //TODO Do this in a more efficient way
    for( unsigned int j = 0 ; j < numneu ; ++j ) {
      unsigned int k = (i*numneu) + j;
      dmat[i][j] = decoherence_matrix[k];
    }
  }
  return dmat;
}




// Function to return the Gamma decoherence matrix (basically D[rho] without the rho part)
// TODO Can I access the rho part here?
squids::SU_vector nuSQUIDSDecoh::DecohGamma(unsigned int ei,unsigned int index_rho, double t) const {
  //TODO Should this be rotated to the flavor basis?????
  //squids::SU_vector decoherenceMatrix{ std::vector<double>{ decoherenceParam,0.,0., 0.,decoherenceParam,0., 0.,0.,decoherenceParam } };
  //squids::SU_vector decoherenceMatrix{ std::vector<double>{ 0.,decoherenceParam,decoherenceParam, decoherenceParam,0.,decoherenceParam, decoherenceParam,decoherenceParam,0. } };
  return decoherence_matrix;
}


} // close namespace


