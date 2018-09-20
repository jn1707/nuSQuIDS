/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/

//TODO Create a .cpp file


#ifndef NUSQUIDSDECOH_H
#define NUSQUIDSDECOH_H

#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids{


/*
  nuSQuIDS extedned to include decohernece
*/

class nuSQUIDSDecoh: public nuSQUIDS {


  public:

    //Default void constructor (required by boost python).
    nuSQUIDSDecoh(){}

    // Constructor : Single energy mode
    nuSQUIDSDecoh(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu,NT)
      , decoherence_gamma_matrix(numneu)
      , n_energy(0.)
    { init(); }

    // Constructor : Multi-energy mode
    nuSQUIDSDecoh(marray<double,1> E_vector, unsigned int numneu, NeutrinoType NT = both,
       bool iinteraction = false, std::shared_ptr<NeutrinoCrossSections> ncs = nullptr)
      : nuSQUIDS(E_vector,numneu,NT,iinteraction,ncs)
      , decoherence_gamma_matrix(numneu)
      , n_energy(0.)
    { init(); }

    // TODO Copy and move constructors...
  
    // Set the Gamma decoherence matrix elements individually
    void Set_DecoherenceGammaMatrix(const marray<double,2>& dmat);

    // Set the Gamma decoherence matrix elements using the Gamma21,31,32 model
    // See https://arxiv.org/pdf/1708.05495.pdf equation 12
    void Set_DecoherenceGammaMatrix(double Gamma); // 2 flavor
    void Set_DecoherenceGammaMatrix(double Gamma21,double Gamma31,double Gamma32); // 3 flavor

    // Get the current value of the decoherence matrix
    marray<double,2> Get_DecoherenceGammaMatrix() const;

    // Set the energy dependence using the form:  Gamma = Gamma_0 * (E/E_0) ^ n_energy
    // See https://arxiv.org/abs/1803.04438.pdf equation 2.13
    void Set_EnergyDependence(double n_energy);

    // Get the current value of the energy depedence 
    double Get_EnergyDependence() const;

    //TODO REMOVE, these are for debugging
    void PrintTransformationMatrix() const;
    void PrintState() const;


  protected:

    // Common initialisation tasks
    void init();

  private:

    // The matrix containing the Gamma coefficients in D[rho]
    // This is the general case of nsun*nsun independent elements,
    // which also supports the Gamma21,Gamma31,Gamma32 model when the correct setter is used
    squids::SU_vector decoherence_gamma_matrix; //TODO Is this the correct type to store this as?
    std::unique_ptr<gsl_matrix_complex> decoherence_gamma_gsl_matrix;

    // Function to return D[rho] (result of decohernce operator on rho)
    squids::SU_vector D_Rho(unsigned int ei,unsigned int index_rho, double t) const override;

    // Buffer for holding energy power (n_energy)
    double n_energy;

};


/*
  nuSQUIDSAtm extended to include decoherence
*/

class nuSQUIDSAtmDecoh : public nuSQUIDSAtm<nuSQUIDSDecoh> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSDecoh>::nuSQUIDSAtm;

    // TODO Copy and move constructors...

    // Wrap debug setter
    void Set_Debug(bool debug) {
      for(nuSQUIDSDecoh& nsq : nusq_array) nsq.Set_Debug(debug);
    } 

    // Wrap decoherence Gamma matrix getters/setters
    void Set_DecoherenceGammaMatrix(const marray<double,2>& dmat) {
      for(nuSQUIDSDecoh& nsq : nusq_array) nsq.Set_DecoherenceGammaMatrix(dmat);
    }

    void Set_DecoherenceGammaMatrix(double Gamma21,double Gamma31,double Gamma32) {
      for(nuSQUIDSDecoh& nsq : nusq_array) nsq.Set_DecoherenceGammaMatrix(Gamma21,Gamma31,Gamma32);
    }

    marray<double,2> Get_DecoherenceGammaMatrix() const {
      return nusq_array[0].Get_DecoherenceGammaMatrix();
    }

    // Wrap decoherence energy dependence getters/setters
    void Set_EnergyDependence(double n_energy) {
      for(nuSQUIDSDecoh& nsq : nusq_array) nsq.Set_EnergyDependence(n_energy);
    }

    double Get_EnergyDependence() const {
      return nusq_array[0].Get_EnergyDependence();
    }


};



} // close namespace

#endif //NUSQUIDSDECOH_H

