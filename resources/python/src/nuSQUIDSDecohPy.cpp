#include "nuSQUIDSpy.h"
#include <nuSQuIDS/nuSQuIDSDecoh.h>

using namespace nusquids;

//TODO use use standard atm version?
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSDecohAtm_Set_initial_state,Set_initial_state,1,2)
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSDecohAtm_EvalFlavor_overload,EvalFlavor,3,5)


BOOST_PYTHON_MODULE(nuSQUIDSDecohPy)
{

    //
    // nuSQUIDSDecoh
    //

    // Register all standard nuSQuIDS and nuSQuIDS atmospheric functions for the user class
    auto nusquids_decoh_register = RegisterBasicNuSQuIDSPythonBindings<nuSQUIDSDecoh>("nuSQUIDSDecoh");

    // Register additional functions or members of the user class
    auto nusquids_decoh_class_object = nusquids_decoh_register.GetClassObject();

    nusquids_decoh_class_object->def("Set_DecoherenceGammaMatrix",(void(nuSQUIDSDecoh::*)(const marray<double,2>&))&nuSQUIDSDecoh::Set_DecoherenceGammaMatrix);
    nusquids_decoh_class_object->def("Set_DecoherenceGammaMatrixDiagonal",(void(nuSQUIDSDecoh::*)(const marray<double,1>&))&nuSQUIDSDecoh::Set_DecoherenceGammaMatrixDiagonal);
    nusquids_decoh_class_object->def("Get_DecoherenceGammaMatrix",&nuSQUIDSDecoh::Get_DecoherenceGammaMatrix);

    nusquids_decoh_class_object->def("EnableDecoherence",&nuSQUIDSDecoh::EnableDecoherence);

    nusquids_decoh_class_object->def("Set_DecoherenceGammaEnergyDependence",&nuSQUIDSDecoh::Set_DecoherenceGammaEnergyDependence);
    nusquids_decoh_class_object->def("Get_DecoherenceGammaEnergyDependence",&nuSQUIDSDecoh::Get_DecoherenceGammaEnergyDependence);

    nusquids_decoh_class_object->def("Set_DecoherenceGammaEnergyScale",&nuSQUIDSDecoh::Set_DecoherenceGammaEnergyScale);
    nusquids_decoh_class_object->def("Get_DecoherenceGammaEnergyScale",&nuSQUIDSDecoh::Get_DecoherenceGammaEnergyScale);


    // nusquids_decoh_atm_class_object->def("Set_DecoherenceGammaEnergyScale",&nuSQUIDSDecohAtm::Get_DecoherenceGammaEnergyScale);


    //
    // nuSQUIDSDecohAtm
    //

    //TODO Ideally would use `RegisterBasicAtmNuSQuIDSPythonBindings` for this, but having trouble with it (either due to 
    // return type of the derived class, or non-bound vector<nuSQuIDSDecoh>). Am pushed for time, so going quick and dirty
    // and directly implenting the bindings for my class here (very poor OAOO beahviour...)

    // auto nusquids_decoh_atm_register = RegisterBasicAtmNuSQuIDSPythonBindings<nuSQUIDSDecohAtm>("nuSQUIDSDecohAtm");
    // auto nusquids_decoh_atm_class_object = nusquids_decoh_atm_register.GetClassObject();

    class_<nuSQUIDSDecohAtm, boost::noncopyable, std::shared_ptr<nuSQUIDSDecohAtm> >("nuSQUIDSDecohAtm", no_init)
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType>(args("CosZenith_vector","E_vector","numneu","NT")))
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction")))
        .def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<NeutrinoCrossSections>>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction","ncs")))
        .def(init<std::string>(args("filename")))
        .def("EvolveState",&nuSQUIDSDecohAtm::EvolveState)
        .def("Set_TauRegeneration",&nuSQUIDSDecohAtm::Set_TauRegeneration)
        //.def("EvalFlavor",&nuSQUIDSDecohAtm::EvalFlavor)
        // .def("EvalFlavor",(double(nuSQUIDSDecohAtm::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSDecohAtm::EvalFlavor, nuSQUIDSAtm_EvalFlavor_overload(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"), "Reads an HDF5 file and loads the contents into the current object."))
        .def("EvalFlavor",(double(nuSQUIDSDecohAtm::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSDecohAtm::EvalFlavor, nuSQUIDSDecohAtm_EvalFlavor_overload<nuSQUIDSDecohAtm>(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"), "nuSQuIDSAtm evaluate flux.."))
        .def("Set_EvalThreads",&nuSQUIDSDecohAtm::Set_EvalThreads)
        .def("Get_EvalThreads",&nuSQUIDSDecohAtm::Get_EvalThreads)
        .def("Set_EarthModel",&nuSQUIDSDecohAtm::Set_EarthModel)
        .def("Set_Body",&nuSQUIDSDecohAtm::Set_Body)
        // .def("WriteStateHDF5",&nuSQUIDSDecohAtm::WriteStateHDF5) // Haven't implemented this
        // .def("ReadStateHDF5",&nuSQUIDSDecohAtm::ReadStateHDF5) // Haven't implemented this
        .def("Set_MixingAngle",&nuSQUIDSDecohAtm::Set_MixingAngle)
        .def("Set_CPPhase",&nuSQUIDSDecohAtm::Set_CPPhase)
        .def("Set_SquareMassDifference",&nuSQUIDSDecohAtm::Set_SquareMassDifference)
        .def("Set_ProgressBar",&nuSQUIDSDecohAtm::Set_ProgressBar)
        .def("Set_MixingParametersToDefault",&nuSQUIDSDecohAtm::Set_MixingParametersToDefault)
        // .def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP)
        .def("Set_rel_error",(void(nuSQUIDSDecohAtm::*)(double))&nuSQUIDSDecohAtm::Set_rel_error)
        .def("Set_rel_error",(void(nuSQUIDSDecohAtm::*)(double, unsigned int))&nuSQUIDSDecohAtm::Set_rel_error)
        .def("Set_abs_error",(void(nuSQUIDSDecohAtm::*)(double))&nuSQUIDSDecohAtm::Set_abs_error)
        .def("Set_abs_error",(void(nuSQUIDSDecohAtm::*)(double, unsigned int))&nuSQUIDSDecohAtm::Set_abs_error)
        .def("GetNumE",&nuSQUIDSDecohAtm::GetNumE)
        .def("GetNumCos",&nuSQUIDSDecohAtm::GetNumCos)
        .def("GetNumNeu",&nuSQUIDSDecohAtm::GetNumNeu)
        .def("GetNumRho",&nuSQUIDSDecohAtm::GetNumRho)
        //.def("EvalMass",(double(nuSQUIDS::*)(unsigned int,double,unsigned int) const)&nuSQUIDS::EvalMass)
        .def("GetnuSQuIDS",(nuSQUIDSDecoh&(nuSQUIDSDecohAtm::*)(unsigned int))&nuSQUIDSDecohAtm::GetnuSQuIDS,boost::python::return_internal_reference<>())
        // .def("Set_initial_state",(void(nuSQUIDSDecohAtm::*)(const marray<double,3>&, Basis))&nuSQUIDSDecohAtm::Set_initial_state,nuSQUIDSDecohAtm_Set_initial_state())
        // .def("Set_initial_state",(void(nuSQUIDSDecohAtm::*)(const marray<double,4>&, Basis))&nuSQUIDSDecohAtm::Set_initial_state,nuSQUIDSDecohAtm_Set_initial_state())
        .def("Set_initial_state",(void(nuSQUIDSDecohAtm::*)(const marray<double,3>&, Basis))&nuSQUIDSDecohAtm::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSDecohAtm>())
        .def("Set_initial_state",(void(nuSQUIDSDecohAtm::*)(const marray<double,4>&, Basis))&nuSQUIDSDecohAtm::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSDecohAtm>())
        .def("GetERange",&nuSQUIDSDecohAtm::GetERange)
        .def("EnableDecoherence",&nuSQUIDSDecohAtm::EnableDecoherence)
        .def("GetCosthRange",&nuSQUIDSDecohAtm::GetCosthRange)
        .def("Set_DecoherenceGammaMatrix",(void(nuSQUIDSDecohAtm::*)(const marray<double,2>&))&nuSQUIDSDecohAtm::Set_DecoherenceGammaMatrix)
        .def("Set_DecoherenceGammaMatrixDiagonal",(void(nuSQUIDSDecohAtm::*)(const marray<double,1>&))&nuSQUIDSDecohAtm::Set_DecoherenceGammaMatrixDiagonal)
        .def("Get_DecoherenceGammaMatrix",&nuSQUIDSDecohAtm::Get_DecoherenceGammaMatrix)
        .def("Set_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohAtm::Set_DecoherenceGammaEnergyDependence)
        .def("Get_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohAtm::Get_DecoherenceGammaEnergyDependence)
        .def("Set_DecoherenceGammaEnergyScale",&nuSQUIDSDecohAtm::Set_DecoherenceGammaEnergyScale)
        .def("Get_DecoherenceGammaEnergyScale",&nuSQUIDSDecohAtm::Get_DecoherenceGammaEnergyScale)
    ;



}