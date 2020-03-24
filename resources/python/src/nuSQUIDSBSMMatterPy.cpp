#include "nuSQUIDSpy.h"
#include <nuSQuIDS/nuSQuIDSBSMMatter.h>

using namespace nusquids;

//TODO use use standard atm version?
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSBSMMatterAtm_Set_initial_state,Set_initial_state,1,2)
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSBSMMatterAtm_EvalFlavor_overload,EvalFlavor,3,5)


BOOST_PYTHON_MODULE(nuSQUIDSBSMMatterPy)
{

    //
    // nuSQUIDSBSMMatter
    //

    // Register all standard nuSQuIDS and nuSQuIDS atmospheric functions for the user class
    auto nusquids_decoh_register = RegisterBasicNuSQuIDSPythonBindings<nuSQUIDSBSMMatter>("nuSQUIDSBSMMatter");

    // Register additional functions or members of the user class
    auto nusquids_decoh_class_object = nusquids_decoh_register.GetClassObject();

    nusquids_decoh_class_object->def("Set_VectorNSIMatrix",(void(nuSQUIDSBSMMatter::*)(const marray<double,2>&))&nuSQUIDSBSMMatter::Set_VectorNSIMatrix);
    nusquids_decoh_class_object->def("Get_VectorNSIMatrix",&nuSQUIDSBSMMatter::Get_VectorNSIMatrix);

    nusquids_decoh_class_object->def("Set_ScalarNSIMatrix",(void(nuSQUIDSBSMMatter::*)(const marray<double,2>&))&nuSQUIDSBSMMatter::Set_ScalarNSIMatrix);
    nusquids_decoh_class_object->def("Get_ScalarNSIMatrix",&nuSQUIDSBSMMatter::Get_ScalarNSIMatrix);

    //
    // nuSQUIDSBSMMatterAtm
    //

    //TODO Ideally would use `RegisterBasicAtmNuSQuIDSPythonBindings` for this, but having trouble with it (either due to 
    // return type of the derived class, or non-bound vector<nuSQuIDSNSI>). Am pushed for time, so going quick and dirty
    // and directly implenting the bindings for my class here (very poor OAOO beahviour...)

    // auto nusquids_decoh_atm_register = RegisterBasicAtmNuSQuIDSPythonBindings<nuSQUIDSBSMMatterAtm>("nuSQUIDSBSMMatterAtm");
    // auto nusquids_decoh_atm_class_object = nusquids_decoh_atm_register.GetClassObject();

    class_<nuSQUIDSBSMMatterAtm, boost::noncopyable, std::shared_ptr<nuSQUIDSBSMMatterAtm> >("nuSQUIDSBSMMatterAtm", no_init)
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType>(args("CosZenith_vector","E_vector","numneu","NT")))
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction")))
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<NeutrinoCrossSections>>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction","ncs")))
        .def(init<std::string>(args("filename")))
        .def("EvolveState",&nuSQUIDSBSMMatterAtm::EvolveState)
        .def("Set_TauRegeneration",&nuSQUIDSBSMMatterAtm::Set_TauRegeneration)
        //.def("EvalFlavor",&nuSQUIDSBSMMatterAtm::EvalFlavor)
        // .def("EvalFlavor",(double(nuSQUIDSBSMMatterAtm::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSBSMMatterAtm::EvalFlavor, nuSQUIDSAtm_EvalFlavor_overload(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"), "Reads an HDF5 file and loads the contents into the current object."))
        .def("EvalFlavor",(double(nuSQUIDSBSMMatterAtm::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSBSMMatterAtm::EvalFlavor, nuSQUIDSBSMMatterAtm_EvalFlavor_overload<nuSQUIDSBSMMatterAtm>(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"), "nuSQuIDSAtm evaluate flux.."))
        .def("Set_EvalThreads",&nuSQUIDSBSMMatterAtm::Set_EvalThreads)
        .def("Get_EvalThreads",&nuSQUIDSBSMMatterAtm::Get_EvalThreads)
        .def("Set_EarthModel",&nuSQUIDSBSMMatterAtm::Set_EarthModel)
        .def("Set_Body",&nuSQUIDSBSMMatterAtm::Set_Body)
        // .def("WriteStateHDF5",&nuSQUIDSBSMMatterAtm::WriteStateHDF5) // Haven't implemented this
        // .def("ReadStateHDF5",&nuSQUIDSBSMMatterAtm::ReadStateHDF5) // Haven't implemented this
        .def("Set_MixingAngle",&nuSQUIDSBSMMatterAtm::Set_MixingAngle)
        .def("Set_CPPhase",&nuSQUIDSBSMMatterAtm::Set_CPPhase)
        .def("Set_SquareMassDifference",&nuSQUIDSBSMMatterAtm::Set_SquareMassDifference)
        .def("Set_ProgressBar",&nuSQUIDSBSMMatterAtm::Set_ProgressBar)
        .def("Set_MixingParametersToDefault",&nuSQUIDSBSMMatterAtm::Set_MixingParametersToDefault)
        // .def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP)
        .def("Set_rel_error",(void(nuSQUIDSBSMMatterAtm::*)(double))&nuSQUIDSBSMMatterAtm::Set_rel_error)
        .def("Set_rel_error",(void(nuSQUIDSBSMMatterAtm::*)(double, unsigned int))&nuSQUIDSBSMMatterAtm::Set_rel_error)
        .def("Set_abs_error",(void(nuSQUIDSBSMMatterAtm::*)(double))&nuSQUIDSBSMMatterAtm::Set_abs_error)
        .def("Set_abs_error",(void(nuSQUIDSBSMMatterAtm::*)(double, unsigned int))&nuSQUIDSBSMMatterAtm::Set_abs_error)
        .def("GetNumE",&nuSQUIDSBSMMatterAtm::GetNumE)
        .def("GetNumCos",&nuSQUIDSBSMMatterAtm::GetNumCos)
        .def("GetNumNeu",&nuSQUIDSBSMMatterAtm::GetNumNeu)
        .def("GetNumRho",&nuSQUIDSBSMMatterAtm::GetNumRho)
        //.def("EvalMass",(double(nuSQUIDS::*)(unsigned int,double,unsigned int) const)&nuSQUIDS::EvalMass)
        .def("GetnuSQuIDS",(nuSQUIDSBSMMatter&(nuSQUIDSBSMMatterAtm::*)(unsigned int))&nuSQUIDSBSMMatterAtm::GetnuSQuIDS,boost::python::return_internal_reference<>())
        // .def("Set_initial_state",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,3>&, Basis))&nuSQUIDSBSMMatterAtm::Set_initial_state,nuSQUIDSBSMMatterAtm_Set_initial_state())
        // .def("Set_initial_state",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,4>&, Basis))&nuSQUIDSBSMMatterAtm::Set_initial_state,nuSQUIDSBSMMatterAtm_Set_initial_state())
        .def("Set_initial_state",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,3>&, Basis))&nuSQUIDSBSMMatterAtm::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSBSMMatterAtm>())
        .def("Set_initial_state",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,4>&, Basis))&nuSQUIDSBSMMatterAtm::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSBSMMatterAtm>())
        .def("GetERange",&nuSQUIDSBSMMatterAtm::GetERange)
        .def("GetCosthRange",&nuSQUIDSBSMMatterAtm::GetCosthRange)
        .def("Set_VectorNSIMatrix",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,2>&))&nuSQUIDSBSMMatterAtm::Set_VectorNSIMatrix)
        .def("Get_VectorNSIMatrix",&nuSQUIDSBSMMatterAtm::Get_VectorNSIMatrix)
        .def("Set_ScalarNSIMatrix",(void(nuSQUIDSBSMMatterAtm::*)(const marray<double,2>&))&nuSQUIDSBSMMatterAtm::Set_ScalarNSIMatrix)
        .def("Get_ScalarNSIMatrix",&nuSQUIDSBSMMatterAtm::Get_ScalarNSIMatrix)
    ;

}