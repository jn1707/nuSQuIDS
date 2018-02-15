/*
Implementing a neutrino decoherence model in SQuIDS/nuSQuIDS
Tom Stuttard, Mikkel Jensen (Niels Bohr Institute)
*/


#ifndef NUSQUIDSDECOH_H
#define NUSQUIDSDECOH_H

#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

class nuSQUIDSDecoh: public nuSQUIDS {

  public:

    nuSQUIDSDecoh(NeutrinoType NT = neutrino)
      : nuSQUIDS(3,NT)
      , decoherenceParam(0.)
      {

      //init(); //TODO Think is is already called by the nuSQUIDS constructor so should be able to remove it

      //TODO Add interactions

      Set_DecoherenceTerms(true);

    }
  
    void Set_DecoherenceParam(double opt) { decoherenceParam = opt; }


  private:

    double decoherenceParam;

    // Function to return the Gamma decoherence matrix (basically D[rho] without the rho part)
    squids::SU_vector DecohGamma(unsigned int ei,unsigned int index_rho, double t) const override {
      //TODO Should this be rotated to the flavor basis?????
      //squids::SU_vector decoherenceMatrix{ std::vector<double>{ decoherenceParam,0.,0., 0.,decoherenceParam,0., 0.,0.,decoherenceParam } };
      squids::SU_vector decoherenceMatrix{ std::vector<double>{ 0.,decoherenceParam,decoherenceParam, decoherenceParam,0.,decoherenceParam, decoherenceParam,decoherenceParam,0. } };
      return decoherenceMatrix;
    }

};

#endif //NUSQUIDSDECOH_H

