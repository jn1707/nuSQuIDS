/*
Implementing BSM matter effects in SQuIDS/nuSQuIDS

Tom Stuttard (Niels Bohr Institute)
*/


#ifndef NUSQUIDSDECOH_H
#define NUSQUIDSDECOH_H

#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids{


/*
  nuSQuIDS extended to include NSI
*/

class nuSQUIDSBSMMatter: public nuSQUIDS {

  // template<typename,typename>
  // friend class nuSQUIDSAtmDecoh;

  public:

    //Default void constructor (required by boost python).
    nuSQUIDSBSMMatter(){} //TODO remove now?

    // Move constructor
    // nuSQUIDS(nuSQUIDS&&); //TODO

    // Constructor : Single energy mode
    nuSQUIDSBSMMatter(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu,NT)
      , units()
    { init(); }

    // Constructor : Multi-energy mode
    nuSQUIDSBSMMatter(marray<double,1> E_vector, unsigned int numneu, NeutrinoType NT = both,
       bool iinteraction = false, std::shared_ptr<NeutrinoCrossSections> ncs = nullptr)
      : nuSQUIDS(E_vector,numneu,NT,iinteraction,ncs)
      , units()
    { init(); }

    // Constructor: Load from serialized file (not yet implemented, but required to exist for pybindings)
    nuSQUIDSBSMMatter(std::string hdf5_filename, std::string grp = "/",
      std::shared_ptr<InteractionStructure> int_struct = nullptr)
    {
      throw std::runtime_error("nuSQUIDSBSMMatter::Error::HDF5 serialization not implemented");
    }
  
    // Enable/disable the decoherence D[rho] operator
    void EnableDecoherence(bool enable);

    // Set the scalar NSI matrix elements
    void Set_ScalarNSIMatrix(const marray<double,2>& dmat);
    void Set_VectorNSIMatrix(const marray<double,2>& dmat);

    // Get the current value of the decoherence Gamma matrix
    marray<double,2> Get_ScalarNSIMatrix() const;
    marray<double,2> Get_VectorNSIMatrix() const;

  public: //TODO protected, but fix py binding issue

    // Common initialisation tasks
    void init();

    void AddToPreDerive(double x);

    squids::SU_vector HI(unsigned int ei, unsigned int index_rho) const;

    void iniH0();
    squids::SU_vector H0(double Enu, unsigned int irho) const;

  private:

    double V_prefactor;

    squids::SU_vector vector_nsi_matrix_sun;
    // std::vector<squids::SU_vector> vector_nsi_matrix_sun_evol;

    squids::SU_vector scalar_nsi_matrix_sun;
    // std::vector<squids::SU_vector> scalar_nsi_matrix_sun_evol;

    // Units container
    squids::Const units;

};


/*
  nuSQUIDSAtm extended to include decoherence
*/

class nuSQUIDSBSMMatterAtm : public nuSQUIDSAtm<nuSQUIDSBSMMatter> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSBSMMatter>::nuSQUIDSAtm;

    // TODO Copy and move constructors...

    // Wrap all the getters/setters
    void Set_Debug(bool debug) {
      for(nuSQUIDSBSMMatter& nsq : this->GetnuSQuIDS()) nsq.Set_Debug(debug);
    } 

    void Set_VectorNSIMatrix(const marray<double,2>& dmat) {
      for(nuSQUIDSBSMMatter& nsq : this->GetnuSQuIDS()) nsq.Set_VectorNSIMatrix(dmat);
    }

    marray<double,2> Get_VectorNSIMatrix() {
      return this->GetnuSQuIDS(0).Get_VectorNSIMatrix();
    }

    void Set_ScalarNSIMatrix(const marray<double,2>& dmat) {
      for(nuSQUIDSBSMMatter& nsq : this->GetnuSQuIDS()) nsq.Set_ScalarNSIMatrix(dmat);
    }

    marray<double,2> Get_ScalarNSIMatrix() {
      return this->GetnuSQuIDS(0).Get_ScalarNSIMatrix();
    }

};



} // close namespace

#endif //NUSQUIDSDECOH_H

