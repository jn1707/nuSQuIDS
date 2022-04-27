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

    nusquids_decoh_class_object->def("Set_UseLightconeFluctuations",&nuSQUIDSDecoh::Set_UseLightconeFluctuations);
    nusquids_decoh_class_object->def("Get_UseLightconeFluctuations",&nuSQUIDSDecoh::Get_UseLightconeFluctuations);

    nusquids_decoh_class_object->def("Set_mLengthDependenceIndex",&nuSQUIDSDecoh::Set_mLengthDependenceIndex);
    nusquids_decoh_class_object->def("Get_mLengthDependenceIndex",&nuSQUIDSDecoh::Get_mLengthDependenceIndex);

    nusquids_decoh_class_object->def("Set_dL0",&nuSQUIDSDecoh::Set_dL0);
    nusquids_decoh_class_object->def("Get_dL0",&nuSQUIDSDecoh::Get_dL0);

    nusquids_decoh_class_object->def("Set_L0LengthScale",&nuSQUIDSDecoh::Set_L0LengthScale);
    nusquids_decoh_class_object->def("Get_L0LengthScale",&nuSQUIDSDecoh::Get_L0LengthScale);

    nusquids_decoh_class_object->def("Set_DampingPower",&nuSQUIDSDecoh::Set_DampingPower);
    nusquids_decoh_class_object->def("Get_DampingPower",&nuSQUIDSDecoh::Get_DampingPower);

    nusquids_decoh_class_object->def("Set_EtaParam",&nuSQUIDSDecoh::Set_EtaParam);
    nusquids_decoh_class_object->def("Get_EtaParam",&nuSQUIDSDecoh::Get_EtaParam);
    

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
        .def("Set_UseLightconeFluctuations",&nuSQUIDSDecohAtm::Set_UseLightconeFluctuations)
        .def("Get_UseLightconeFluctuations",&nuSQUIDSDecohAtm::Get_UseLightconeFluctuations)
        .def("Set_mLengthDependenceIndex",&nuSQUIDSDecohAtm::Set_mLengthDependenceIndex)
        .def("Get_mLengthDependenceIndex",&nuSQUIDSDecohAtm::Get_mLengthDependenceIndex)
        .def("Set_dL0",&nuSQUIDSDecohAtm::Set_dL0)
        .def("Get_dL0",&nuSQUIDSDecohAtm::Get_dL0)
        .def("Set_L0LengthScale",&nuSQUIDSDecohAtm::Set_L0LengthScale)
        .def("Get_L0LengthScale",&nuSQUIDSDecohAtm::Get_L0LengthScale)
        .def("Set_DampingPower",&nuSQUIDSDecohAtm::Set_DampingPower)
        .def("Get_DampingPower",&nuSQUIDSDecohAtm::Get_DampingPower)
        .def("Set_EtaParam",&nuSQUIDSDecohAtm::Set_EtaParam)
        .def("Get_EtaParam",&nuSQUIDSDecohAtm::Get_EtaParam)
    ;





    //
    // nuSQUIDSDecohLayers
    //

    // auto nusquids_decoh_layers_register = RegisterBasicLayerNuSQuIDSPythonBindings<nuSQUIDSDecohLayers>("nuSQUIDSDecohLayers");


    //TODO See TODO for nuSQUIDSDecohAtm (same applied here)

    // auto nusquids_decoh_atm_register = RegisterBasicAtmNuSQuIDSPythonBindings<nuSQUIDSDecohAtm>("nuSQUIDSDecohAtm");
    // auto nusquids_decoh_atm_class_object = nusquids_decoh_atm_register.GetClassObject();

    class_<nuSQUIDSDecohLayers, boost::noncopyable, std::shared_ptr<nuSQUIDSDecohLayers> >("nuSQUIDSDecohLayers", no_init)
        // Standard nuSQUIDSLayers functions
        .def(init<marray<double,2>,marray<double,2>,marray<double,2>,marray<double,1>,unsigned int,NeutrinoType>(args("lengths", "densities", "ye", "energies","numneu","NT")))
        .def("EvolveState",&nuSQUIDSDecohLayers::EvolveState)
        .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,1>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSLayers_Set_initial_state<nuSQUIDSDecohLayers>())
        .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,2>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSLayers_Set_initial_state<nuSQUIDSDecohLayers>())
        .def("EvalWithState",
            (
                marray<double,1>(nuSQUIDSDecohLayers::*)(
                unsigned int, const marray<double,1>&, const marray<double,1>&,
                const marray<double,2>&,
                unsigned int, double, double,
                const marray<double,1>&,
                const marray<double,1>&,
                const marray<double,1>&
            )
            )&nuSQUIDSDecohLayers::ArrEvalWithStateTRangeLPFilter,
                args("flavor", "time", "energy", "state", "rho", "avg_cutoff", "avg_scale",
                "t_range", "lowpass_cutoff", "lowpass_scale"),
                "Returns the flavor composition with a given (interpolated) state.\n\n"
                "Same as the other array version of the function, but allows a different\n"
                "averaging time and low-pass filter settings for each evaluation.\n"
                "All other arguments must be passed.\n\n"
            )
        .def("GetStates", (marray<double,2>(nuSQUIDSDecohLayers::*)(unsigned int))&nuSQUIDSDecohLayers::GetStatesArr, nuSQUIDSLayers_GetStates_overload<nuSQUIDSDecohLayers>(args("rho"), "Get evolved states of all nodes."))
        .def("Set_MixingAngle",&nuSQUIDSDecohLayers::Set_MixingAngle)
        .def("Set_CPPhase",&nuSQUIDSDecohLayers::Set_CPPhase)
        .def("Set_SquareMassDifference",&nuSQUIDSDecohLayers::Set_SquareMassDifference)
        .def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_rel_error)
        .def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_rel_error)
        .def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_abs_error)
        .def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_abs_error)
        .def("Set_EvolLowPassCutoff",&nuSQUIDSDecohLayers::Set_EvolLowPassCutoff)
        .def("Set_EvolLowPassScale",&nuSQUIDSDecohLayers::Set_EvolLowPassScale)
        .def("Set_AllowConstantDensityOscillationOnlyEvolution",&nuSQUIDSDecohLayers::Set_AllowConstantDensityOscillationOnlyEvolution)
        .def("Set_EvalThreads",&nuSQUIDSDecohLayers::Set_EvalThreads)

        // Decoherence functions
        .def("Set_DecoherenceGammaMatrix",(void(nuSQUIDSDecohLayers::*)(const marray<double,2>&))&nuSQUIDSDecohLayers::Set_DecoherenceGammaMatrix)
        .def("Set_DecoherenceGammaMatrixDiagonal",(void(nuSQUIDSDecohLayers::*)(const marray<double,1>&))&nuSQUIDSDecohLayers::Set_DecoherenceGammaMatrixDiagonal)
        .def("Get_DecoherenceGammaMatrix",&nuSQUIDSDecohLayers::Get_DecoherenceGammaMatrix)
        .def("Set_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohLayers::Set_DecoherenceGammaEnergyDependence)
        .def("Get_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohLayers::Get_DecoherenceGammaEnergyDependence)
        .def("Set_DecoherenceGammaEnergyScale",&nuSQUIDSDecohLayers::Set_DecoherenceGammaEnergyScale)
        .def("Get_DecoherenceGammaEnergyScale",&nuSQUIDSDecohLayers::Get_DecoherenceGammaEnergyScale)
        .def("Set_UseLightconeFluctuations",&nuSQUIDSDecohLayers::Set_UseLightconeFluctuations)
        .def("Get_UseLightconeFluctuations",&nuSQUIDSDecohLayers::Get_UseLightconeFluctuations)
        .def("Set_mLengthDependenceIndex",&nuSQUIDSDecohLayers::Set_mLengthDependenceIndex)
        .def("Get_mLengthDependenceIndex",&nuSQUIDSDecohLayers::Get_mLengthDependenceIndex)
        .def("Set_dL0",&nuSQUIDSDecohLayers::Set_dL0)
        .def("Get_dL0",&nuSQUIDSDecohLayers::Get_dL0)
        .def("Set_L0LengthScale",&nuSQUIDSDecohLayers::Set_L0LengthScale)
        .def("Get_L0LengthScale",&nuSQUIDSDecohLayers::Get_L0LengthScale)
        .def("Set_DampingPower",&nuSQUIDSDecohLayers::Set_DampingPower)
        .def("Get_DampingPower",&nuSQUIDSDecohLayers::Get_DampingPower)
        .def("Set_EtaParam",&nuSQUIDSDecohLayers::Set_EtaParam)
        .def("Get_EtaParam",&nuSQUIDSDecohLayers::Get_EtaParam)


    //     .def("EvalFlavor",(double(nuSQUIDSDecohLayers::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSDecohLayers::EvalFlavor, nuSQUIDSDecohLayers_EvalFlavor_overload<nuSQUIDSDecohLayers>(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"), "nuSQuIDSAtm evaluate flux.."))
    //     .def("Set_EvalThreads",&nuSQUIDSDecohLayers::Set_EvalThreads)
    //     .def("Get_EvalThreads",&nuSQUIDSDecohLayers::Get_EvalThreads)
    //     .def("Set_EarthModel",&nuSQUIDSDecohLayers::Set_EarthModel)
    //     .def("Set_Body",&nuSQUIDSDecohLayers::Set_Body)
    //     // .def("WriteStateHDF5",&nuSQUIDSDecohLayers::WriteStateHDF5) // Haven't implemented this
    //     // .def("ReadStateHDF5",&nuSQUIDSDecohLayers::ReadStateHDF5) // Haven't implemented this
    //     .def("Set_MixingAngle",&nuSQUIDSDecohLayers::Set_MixingAngle)
    //     .def("Set_CPPhase",&nuSQUIDSDecohLayers::Set_CPPhase)
    //     .def("Set_SquareMassDifference",&nuSQUIDSDecohLayers::Set_SquareMassDifference)
    //     .def("Set_ProgressBar",&nuSQUIDSDecohLayers::Set_ProgressBar)
    //     .def("Set_MixingParametersToDefault",&nuSQUIDSDecohLayers::Set_MixingParametersToDefault)
    //     // .def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP)
    //     .def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_rel_error)
    //     .def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_rel_error)
    //     .def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_abs_error)
    //     .def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_abs_error)
    //     .def("GetNumE",&nuSQUIDSDecohLayers::GetNumE)
    //     .def("GetNumCos",&nuSQUIDSDecohLayers::GetNumCos)
    //     .def("GetNumNeu",&nuSQUIDSDecohLayers::GetNumNeu)
    //     .def("GetNumRho",&nuSQUIDSDecohLayers::GetNumRho)
    //     //.def("EvalMass",(double(nuSQUIDS::*)(unsigned int,double,unsigned int) const)&nuSQUIDS::EvalMass)
    //     .def("GetnuSQuIDS",(nuSQUIDSDecoh&(nuSQUIDSDecohLayers::*)(unsigned int))&nuSQUIDSDecohLayers::GetnuSQuIDS,boost::python::return_internal_reference<>())
    //     // .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,3>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSDecohLayers_Set_initial_state())
    //     // .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,4>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSDecohLayers_Set_initial_state())
    //     .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,3>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSDecohLayers>())
    //     .def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,4>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSDecohLayers>())
    //     .def("GetERange",&nuSQUIDSDecohLayers::GetERange)
    //     .def("EnableDecoherence",&nuSQUIDSDecohLayers::EnableDecoherence)
    //     .def("GetCosthRange",&nuSQUIDSDecohLayers::GetCosthRange)
    //     .def("Set_DecoherenceGammaMatrix",(void(nuSQUIDSDecohLayers::*)(const marray<double,2>&))&nuSQUIDSDecohLayers::Set_DecoherenceGammaMatrix)
    //     .def("Set_DecoherenceGammaMatrixDiagonal",(void(nuSQUIDSDecohLayers::*)(const marray<double,1>&))&nuSQUIDSDecohLayers::Set_DecoherenceGammaMatrixDiagonal)
    //     .def("Get_DecoherenceGammaMatrix",&nuSQUIDSDecohLayers::Get_DecoherenceGammaMatrix)
    //     .def("Set_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohLayers::Set_DecoherenceGammaEnergyDependence)
    //     .def("Get_DecoherenceGammaEnergyDependence",&nuSQUIDSDecohLayers::Get_DecoherenceGammaEnergyDependence)
    //     .def("Set_DecoherenceGammaEnergyScale",&nuSQUIDSDecohLayers::Set_DecoherenceGammaEnergyScale)
    //     .def("Get_DecoherenceGammaEnergyScale",&nuSQUIDSDecohLayers::Get_DecoherenceGammaEnergyScale)
    ;



    //   //class_object->def(init<marray<double,2>,marray<double,2>,marray<double,2>,marray<double,1>,unsigned int,NeutrinoType,bool>(args("lengths", "densities", "ye", "energies","numneu","NT","iinteraction")));
    //   //class_object->def(init<marray<double,2>,marray<double,2>,marray<double,2>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<NeutrinoCrossSections>>(args("lengths", "densities", "ye", "energies","numneu","NT","iinteraction","ncs")));
    //   //class_object->def(init<std::string>(args("filename")));
    //   class_object->def("Set_TauRegeneration",&nuSQUIDSDecohLayers::Set_TauRegeneration);
    //   class_object->def("EvalFlavorAtNode", (double(nuSQUIDSDecohLayers::*)(unsigned int, unsigned int, unsigned int))&nuSQUIDSDecohLayers::EvalFlavorAtNode,
    //     nuSQUIDSLayers_EvalFlavorAtNode_overload<nuSQUIDSDecohLayers>(args("flavor", "node_idx", "rho"), "Evaluate flavor state at a node."));
    //   class_object->def("EvalFlavorAtNodes", (marray<double,1>(nuSQUIDSDecohLayers::*)(unsigned int, unsigned int))&nuSQUIDSDecohLayers::EvalFlavorAtNodes,
    //     nuSQUIDSLayers_EvalFlavorAtNodes_overload<nuSQUIDSDecohLayers>(args("flavor", "rho"), "Evaluate flavor state at all nodes and return as array."));
    //   class_object->def("GetStates", (marray<double,2>(nuSQUIDSDecohLayers::*)(unsigned int))&nuSQUIDSDecohLayers::GetStatesArr,
    //     nuSQUIDSLayers_GetStates_overload<nuSQUIDSDecohLayers>(args("rho"), "Get evolved states of all nodes."));
    //   class_object->def("EvalWithState",
    //     (
    //       double(nuSQUIDSDecohLayers::*)(
    //         unsigned int, double, double,
    //         const marray<double,1>&, unsigned int,
    //         double, double, double, double, double
    //       )
    //     )&nuSQUIDSDecohLayers::EvalWithState,
    //     nuSQUIDSLayers_EvalWithState_overload<nuSQUIDSDecohLayers>(
    //       args("flavor", "time", "energy", "state", "rho", "avg_cutoff", "avg_scale",
    //         "t_range", "lowpass_cutoff", "lowpass_scale"),
    //       "Returns the flavor composition with a given (interpolated) state.\n\n"
    //       "Arguments\n"
    //       "---------\n"
    //       "  flavor (int): Flavor to be evaluated.\n"
    //       "  time (double): Total evolution time. Should match the total distance of \n"
    //       "      the state evolution.\n"
    //       "  enu (double): Neutrino energy [eV].\n"
    //       "  state (1D ndarray of double): Interaction picture state to calculate the\n"
    //       "      trace with.\n"
    //       "  rho (int): Index of neutrino type. If neutrino type is `both`, neutrinos\n"
    //       "      are 0 and antineutrinos are 1. If neutrino type is not `both`, the\n"
    //       "      index is always 0.\n"
    //       "  avg_cutoff (double): Scale to use for averaging fast oscillations. This is\n"
    //       "      the number of oscillations at which the evaluation of sine and cosine\n"
    //       "      terms is cut off.\n"
    //       "  avg_scale (double): Distance from `avg_cutoff` over which to apply linear\n"
    //       "      ramp to smooth the transition between averaged and non-averaged regions.\n"
    //       "  t_range (double): Time or distance over which to average oscillations.\n"
    //       "      This averaging is done numerically in SQuIDS. In compatible with\n"
    //       "      `avg_cutoff`.\n"
    //       "  lowpass_cutoff (double): Frequency cut-off of the low-pass filter. Sine\n"
    //       "      and cosine evaluations at a higher frequency are evaluated as zero.\n"
    //       "      This is distinct from the effect of `avg_cutoff` because the number of\n"
    //       "      completed oscillations doesn't matter, i.e. there is no\n"
    //       "      time/distance dependence of this filter.\n"
    //       "  lowpass_scale (double): Distance in frequency over which a linear ramp is\n"
    //       "      applied in the low-pass filter. Frequencies below\n"
    //       "      `(lowpass_cutoff - lowpass_scale)` are allowed to pass fully, and the\n"
    //       "      linear ramp scales down to zero at `lowpass_cutoff`. If this scale is\n"
    //       "      zero, the filter is a step-function.\n"
    //     )
    //   );
    //   class_object->def("EvalWithState",
    //     (
    //       marray<double,1>(nuSQUIDSDecohLayers::*)(
    //         unsigned int, const marray<double,1>&, const marray<double,1>&,
    //         const marray<double,2>&, unsigned int, double, double, double, double, double
    //       )
    //     )&nuSQUIDSDecohLayers::ArrEvalWithState,
    //     nuSQUIDSLayers_ArrEvalWithState_overload<nuSQUIDSDecohLayers>(
    //       args("flavor", "time", "energy", "state", "rho", "avg_cutoff", "avg_scale",
    //         "t_range", "lowpass_cutoff", "lowpass_scale"),
    //       "Returns the flavor composition with a given (interpolated) state.\n\n"
    //       "This is the array version of the function, taking in numpy arrays of time,\n"
    //       "energy, and states.\n\n"
    //       "Arguments\n"
    //       "---------\n"
    //       "  flavor (int): Flavor to be evaluated.\n"
    //       "  time (1D ndarray of double): Total evolution time. Should match the total\n"
    //       "      distance of the state evolution.\n"
    //       "  enu (1D ndarray of double): Neutrino energy [eV].\n"
    //       "  state (2D ndarray of double): Interaction picture states to calculate the\n"
    //       "      trace with. The first dimension of the array has to match the number\n"
    //       "      of times and energies.\n"
    //       "  rho (int): Index of neutrino type. If neutrino type is `both`, neutrinos\n"
    //       "      are 0 and antineutrinos are 1. If neutrino type is not `both`, the\n"
    //       "      index is always 0.\n"
    //       "  avg_cutoff (double): Scale to use for averaging fast oscillations. This is\n"
    //       "      the number of oscillations at which the evaluation of sine and cosine\n"
    //       "      terms is cut off.\n"
    //       "  avg_scale (double): Distance from `avg_cutoff` over which to apply linear\n"
    //       "      ramp to smooth the transition between averaged and non-averaged regions.\n"
    //       "  t_range (double): Time or distance over which to average oscillations. This\n"
    //       "      averaging is done numerically in SQuIDS. In compatible with\n"
    //       "       `avg_cutoff`.\n"
    //       "  lowpass_cutoff (double): Frequency cut-off of the low-pass filter. Sine\n"
    //       "      and cosine evaluations at a higher frequency are evaluated as zero.\n"
    //       "      This is distinct from the effect of `avg_cutoff` because the number\n"
    //       "      of completed oscillations doesn't matter, i.e. there is no\n"
    //       "      time/distance dependenceof this filter.\n"
    //       "  lowpass_scale (double): Distance in frequency over which a linear ramp is\n"
    //       "      applied in the low-pass filter. Frequencies below\n"
    //       "      `(lowpass_cutoff - lowpass_scale)` are allowed to pass fully, and the\n"
    //       "      linear ramp scales down to zero at `lowpass_cutoff`. If this scale is\n"
    //       "      zero, the filter is a step-function.\n"
    //     )
    //   );
    //   class_object->def("EvalWithState",
    //     (
    //       marray<double,1>(nuSQUIDSDecohLayers::*)(
    //         unsigned int, const marray<double,1>&, const marray<double,1>&,
    //         const marray<double,2>&, unsigned int, double, double,
    //         const marray<double,1>&, double, double
    //       )
    //     )&nuSQUIDSDecohLayers::ArrEvalWithStateTRange,
    //     nuSQUIDSLayers_ArrEvalWithStateTRange_overload<nuSQUIDSDecohLayers>(
    //       args("flavor", "time", "energy", "state", "rho", "avg_cutoff", "avg_scale",
    //         "t_range", "lowpass_cutoff", "lowpass_scale"),
    //       "Returns the flavor composition with a given (interpolated) state.\n\n"
    //       "Same as the other array version of the function, but allows a different\n"
    //       "averaging time for each evaluation. All other arguments must be passed.\n\n"
    //     )
    //   );
    //   class_object->def("EvalWithState",
    //     (
    //       marray<double,1>(nuSQUIDSDecohLayers::*)(
    //         unsigned int, const marray<double,1>&, const marray<double,1>&,
    //         const marray<double,2>&,
    //         unsigned int, double, double,
    //         const marray<double,1>&,
    //         const marray<double,1>&,
    //         const marray<double,1>&
    //       )
    //     )&nuSQUIDSDecohLayers::ArrEvalWithStateTRangeLPFilter,
    //       args("flavor", "time", "energy", "state", "rho", "avg_cutoff", "avg_scale",
    //         "t_range", "lowpass_cutoff", "lowpass_scale"),
    //       "Returns the flavor composition with a given (interpolated) state.\n\n"
    //       "Same as the other array version of the function, but allows a different\n"
    //       "averaging time and low-pass filter settings for each evaluation.\n"
    //       "All other arguments must be passed.\n\n"
    //   );
    //   class_object->def("Set_EvalThreads",&nuSQUIDSDecohLayers::Set_EvalThreads);
    //   class_object->def("Get_EvalThreads",&nuSQUIDSDecohLayers::Get_EvalThreads);
    //   class_object->def("Set_EvolLowPassCutoff",&nuSQUIDSDecohLayers::Set_EvolLowPassCutoff);
    //   class_object->def("Set_EvolLowPassScale",&nuSQUIDSDecohLayers::Set_EvolLowPassScale);
    //   // class_object->def("WriteStateHDF5",&nuSQUIDSDecohLayers::WriteStateHDF5);
    //   // class_object->def("ReadStateHDF5",&nuSQUIDSDecohLayers::ReadStateHDF5);
    //   class_object->def("Set_MixingAngle",&nuSQUIDSDecohLayers::Set_MixingAngle);
    //   class_object->def("Get_MixingAngle",&nuSQUIDSDecohLayers::Get_MixingAngle);
    //   class_object->def("Set_CPPhase",&nuSQUIDSDecohLayers::Set_CPPhase);
    //   class_object->def("Get_CPPhase",&nuSQUIDSDecohLayers::Get_CPPhase);
    //   class_object->def("Set_SquareMassDifference",&nuSQUIDSDecohLayers::Set_SquareMassDifference);
    //   class_object->def("Get_SquareMassDifference",&nuSQUIDSDecohLayers::Get_SquareMassDifference);
    //   class_object->def("Set_h",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_h);
    //   class_object->def("Set_h",(void(nuSQUIDSDecohLayers::*)(double,unsigned int))&nuSQUIDSDecohLayers::Set_h);
    //   class_object->def("Set_h_max",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_h_max);
    //   class_object->def("Set_h_max",(void(nuSQUIDSDecohLayers::*)(double,unsigned int))&nuSQUIDSDecohLayers::Set_h_max);
    //   class_object->def("Set_h_min",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_h_min);
    //   class_object->def("Set_h_min",(void(nuSQUIDSDecohLayers::*)(double,unsigned int))&nuSQUIDSDecohLayers::Set_h_min);
    //   //class_object->def("Set_ProgressBar",&nuSQUIDSDecohLayers::Set_ProgressBar);
    //   class_object->def("Set_MixingParametersToDefault",&nuSQUIDSDecohLayers::Set_MixingParametersToDefault);
    //   class_object->def("Set_GSL_step",wrap_nusqlayer_Set_GSL_STEP<BaseType>);
    //   class_object->def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_rel_error);
    //   class_object->def("Set_rel_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_rel_error);
    //   class_object->def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double))&nuSQUIDSDecohLayers::Set_abs_error);
    //   class_object->def("Set_abs_error",(void(nuSQUIDSDecohLayers::*)(double, unsigned int))&nuSQUIDSDecohLayers::Set_abs_error);
    //   class_object->def("GetNumNeu",&nuSQUIDSDecohLayers::GetNumNeu);
    //   class_object->def("GetNumRho",&nuSQUIDSDecohLayers::GetNumRho);
    //   class_object->def("GetnuSQuIDS",(std::vector<BaseType>&(nuSQUIDSDecohLayers::*)())&nuSQUIDSDecohLayers::GetnuSQuIDS,boost::python::return_internal_reference<>());
    //   class_object->def("GetnuSQuIDS",(BaseType&(nuSQUIDSDecohLayers::*)(unsigned int))&nuSQUIDSDecohLayers::GetnuSQuIDS,boost::python::return_internal_reference<>());
    //   // TODO Do we want to handle neutrinos and antineutrinos at the same time?
    //   class_object->def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,1>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSLayers_Set_initial_state<nuSQUIDSDecohLayers>());
    //   class_object->def("Set_initial_state",(void(nuSQUIDSDecohLayers::*)(const marray<double,2>&, Basis))&nuSQUIDSDecohLayers::Set_initial_state,nuSQUIDSLayers_Set_initial_state<nuSQUIDSDecohLayers>());
    //   class_object->def("Set_IncludeOscillations",&nuSQUIDSDecohLayers::Set_IncludeOscillations);
    //   class_object->def("Set_GlashowResonance",&nuSQUIDSDecohLayers::Set_GlashowResonance);
    //   class_object->def("Set_TauRegeneration",&nuSQUIDSDecohLayers::Set_TauRegeneration);
    //   class_object->def("Set_AllowConstantDensityOscillationOnlyEvolution",&nuSQUIDSDecohLayers::Set_AllowConstantDensityOscillationOnlyEvolution);
    //   class_object->def("Set_PositivyConstrain",&nuSQUIDSDecohLayers::Set_PositivityConstrain);
    //   class_object->def("Set_PositivyConstrainStep",&nuSQUIDSDecohLayers::Set_PositivityConstrainStep);
    //   class_object->def("Get_EvalThreads",&nuSQUIDSDecohLayers::Get_EvalThreads);
    //   class_object->def("Set_EvalThreads",&nuSQUIDSDecohLayers::Set_EvalThreads);
    //   class_object->def("SetNeutrinoCrossSections",&nuSQUIDSDecohLayers::SetNeutrinoCrossSections);
    //   class_object->def("GetNeutrinoCrossSections",&nuSQUIDSDecohLayers::GetNeutrinoCrossSections);
    // }




}