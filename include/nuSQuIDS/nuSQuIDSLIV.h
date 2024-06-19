#ifndef NUSQUIDLIV_H
#define NUSQUIDLIV_H

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include <iomanip>
#include <math.h>
#include <cmath>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>


namespace nusquids {


//
// nuSQuIDS LIV model
//

class nuSQUIDSLIV: public nuSQUIDS {

  private:

    squids::SU_vector CPT_odd_Eindep_evol;
    std::vector<squids::SU_vector> CPT_odd_Eindep_evol;

    squids::SU_vector CPT_odd_Edep_evol;
    std::vector<squids::SU_vector> CPT_odd_Edep_evol;

    squids::SU_vector CPT_even_Eindep_evol;
    std::vector<squids::SU_vector> CPT_even_Eindep_evol;

    squids::SU_vector CPT_even_Edep_evol;
    std::vector<squids::SU_vector> CPT_even_Edep_evol;

    // Override nuSQuIDS functions to add SME terms
    void AddToPreDerive(double x);
    squids:: SU_vector HI(unsigned int ie,unsigned int irho) const;

    // Common init function, called by various constructors
    void init() ;


  public:


    //
    // Constructors
    //

    //Default void constructor (never used but required by boost python, specifically RegisterBasicNuSQuIDSPythonBindings, not sure why).
    nuSQUIDSLIV() {}

    /// Multiple energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(marray<double,1> E_vector,unsigned int numneu, NeutrinoType NT = both, bool iinteraction = false, std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
      : nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs)
      // : nuSQUIDS(E_vector, numneu, NT, false, nullptr)
    { init(); }

    /// Single energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu, NT)
    { init(); }

    // Load from serialized file (not yet implemented, but required to exist for pybindings)
    nuSQUIDSLIV(std::string hdf5_filename, std::string grp = "/", std::shared_ptr<InteractionStructure> int_struct = nullptr)
      // : nuSQUIDS(hdf5_filename, grp, int_struct)
    {
      throw std::runtime_error("nuSQUIDSLIV::Error::HDF5 serialization not implemented");
    }



    //
    // SME param getters/setters
    //

     void Set_LIVCoefficient(const marray<double,3>& a_mat, const marray<double,3>& c_mat, double ra_rad, double dec_rad);

};


/*
  nuSQUIDSAtm extended to include LIV
*/

class nuSQUIDSLIVAtm : public nuSQUIDSAtm<nuSQUIDSLIV> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSLIV>::nuSQUIDSAtm;

    // Wrap all the getters/setters
    void Set_LIVCoefficient(const marray<double,3>& a_mat, const marray<double,3>& c_mat, double ra_rad, double dec_rad){
      for(nuSQUIDSLIV& nsq : this->GetnuSQuIDS()) nsq.Set_LIVCoefficient(a_mat, c_mat, ra_rad, dec_rad);
    } 

};

} // close namespace

#endif //NUSQUIDLIV_H
