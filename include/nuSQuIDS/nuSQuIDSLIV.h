#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>



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



namespace nusquids {

class nuSQUIDSLIV: public nuSQUIDS {

  private:

    // The parameters of the LIV model
    int LIV_n;
    double LIV_cft;

    squids::SU_vector LIVP;
    std::vector<squids::SU_vector> LIVP_evol;

    void AddToPreDerive(double x){
      //std::cout << "+++ AddToPreDerive : start" << std::endl;
      for(unsigned int ei = 0; ei < ne; ei++){
        // asumming same mass hamiltonian for neutrinos/antineutrinos
        //std::cout << "+++ AddToPreDerive : " << ei << " : Before h0" << std::endl;
        squids::SU_vector h0 = H0( E_range[ei], 0 );
        //std::cout << "+++ AddToPreDerive : " << ei << " : After h0, before evolve" << std::endl;
        LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );
        //std::cout << "+++ AddToPreDerive : " << ei << " : After evolve" << std::endl;
      }
      //std::cout << "+++ AddToPreDerive : end" << std::endl;
    }


    // squids::SU_vector H0(double E, unsigned int irho) const {
    //   squids::SU_vector h0 = nuSQUIDS::H0( E, irho );
    //   squids::SU_vector h0_liv = h0 + ( LIV_cft * pow(E, LIV_n) );
    //   return h0_liv;
    // }


    squids:: SU_vector HI(unsigned int ie,unsigned int irho) const {
      //std::cout << "+++ HI : start" << std::endl;
      squids::SU_vector potential = nuSQUIDS::HI(ie, irho);
      double sign = 1;
      // if ((irho == 1 and NT==both) or NT==antineutrino){
      //     // antineutrino matter potential flips sign
      //     sign*=(-1);
      // }
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      potential += sign*pow(E_range[ie], LIV_n) * LIVP_evol[ie]; // <- super important line here is where all the physics is set
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      //std::cout << "+++ HI : end" << std::endl;
      return potential;
    }

  void init() {
    /*
      Common init function, called by various constructors
    */

    //std::cout << "+++ nuSQUIDSLIV::init" << std::endl;

    // Init LIV params to null values
    LIV_n = 0;
    LIV_cft = 0;

    // Allocate some matrices
    LIVP_evol.resize(ne);
    for(unsigned int ei = 0; ei < ne; ei++){
      LIVP_evol[ei] = squids::SU_vector(nsun);
    }

  }


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

    void Set_LIVCoefficient(double cft){
      
      gsl_complex c{{ cft , 0.0 }}; //Only using real part right now


       // defining a complex matrix M which will contain our flavor
       // violating flavor structure.
       gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3); //TODO check num nu 
       // gsl_matrix_complex_set(M, 0, 0, c);
       // gsl_matrix_complex_set(M, 1, 1, c);
       // gsl_matrix_complex_set(M, 2, 2, c);

       // gsl_matrix_complex_set(M, 0, 0, c); //TODO This needs fixing - what is correct structure?
       // gsl_matrix_complex_set(M, 1, 1, c);
       gsl_matrix_complex_set(M, 2, 2, c);

       // gsl_matrix_complex_set(M, 1, 0, c);
       // gsl_matrix_complex_set(M, 0, 1, c);
       // gsl_matrix_complex_set(M, 2, 1, c);
       // gsl_matrix_complex_set(M, 1, 2, c);

       LIVP = squids::SU_vector(M);

       std::cout << "Matrix:" << std::endl;
       print_gsl_matrix(M);

       // rotate from flavor to mass basis
       // LIVP.RotateToB1(params);

       // free allocated matrix
       gsl_matrix_complex_free(M);


       // c_params = lv_params;
       // LIV_params_set = true;
    }

    // void Set_LIVCoefficient(double cft){ //TODO complex number?
    //   LIV_cft = cft;
    // }

    void Set_LIVEnergyDependence(int n){
      LIV_n = n;
    }

};

} // close nusquids namespace

#endif //nusquidslv_h

