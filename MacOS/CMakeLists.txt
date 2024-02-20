set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("./src/")
include_directories("./include/")
include_directories("./include/CHEP/Evaluator")
include_directories("./include/CHEP/Geometry")
include_directories("./include/CHEP/Random")
include_directories("./include/CHEP/Units")
include_directories("./include/CHEP/Utility")
include_directories("./include/CHEP/Vector")
include_directories("./include/PTL/")
include_directories("./SIPM/")

include_directories("/usr/local/gsl-2.7.1/")
link_directories("/usr/local/lib/")

project(SiPM LANGUAGES CXX)

find_package(LibXml2 REQUIRED)

add_executable(SiPM main.cpp  ./src/GateXMLDocument.cc ./src/TreeFile.cc
    ./SIPM/SiPM.hh ./SIPM/SiPM.cc  ./include/TreeFile.hh ./include/GateXMLDocument.hh
    ./src/NumpyFile.cc ./include/NumpyFile.hh ./src/cnpy.cpp ./include/cnpy.h
    ./include/G4UnitsTable.hh  ./include/G4Tokenizer.hh ./include/G4UIcmdWithADoubleAndUnit.hh
    ./src/G4UIcmdWithADoubleAndUnit.cc  ./src/G4UnitsTable.cc ./src/G4Threading.cc ./include/G4Threading.hh
    ./src/G4ThreadLocalSingleton.cc ./include/G4ThreadLocalSingleton.hh  ./include/G4CacheDetails.hh
./src/G4CacheDetails.cc ./src/G4UIcmdWithADouble.cc ./include/G4UIcmdWithADouble.hh
./src/G4UIcommand.cc ./include/G4UIcommand.hh ./src/G4UImessenger.cc ./include/G4UImessenger.hh
 ./src/G4UIdirectory.cc ./include/G4UIdirectory.hh  ./src/G4UImanager.cc ./include/G4UImanager.hh
 ./src/G4VStateDependent.cc ./include/G4VStateDependent.hh  ./src/G4StateManager.cc ./include/G4StateManager.hh
 ./src/G4UIcommandTree.cc ./include/G4UIcommandTree.hh ./src/G4UIsession.cc ./include/G4UIsession.hh
 ./src/G4coutDestination.cc ./include/G4coutDestination.hh  ./src/G4UIbatch.cc ./include/G4UIbatch.hh
 ./src/G4UIcontrolMessenger.cc ./include/G4UIcontrolMessenger.hh  ./src/G4UnitsMessenger.cc ./include/G4UnitsMessenger.hh
./src/G4LocalThreadCoutMessenger.cc ./include/G4LocalThreadCoutMessenger.hh ./src/G4ProfilerMessenger.cc ./include/G4ProfilerMessenger.hh
./src/G4UIaliasList.cc ./include/G4UIaliasList.hh ./src/G4MTcoutDestination.cc ./include/G4MTcoutDestination.hh ./src/G4UIbridge.cc ./include/G4UIbridge.hh
./src/G4UIcmdWithAString.cc ./include/G4UIcmdWithAString.hh ./src/G4UIcmdWithABool.cc ./include/G4UIcmdWithABool.hh
./src/G4UIcmdWithAnInteger.cc ./include/G4UIcmdWithAnInteger.hh ./src/G4UIcmdWithoutParameter.cc ./include/G4UIcmdWithoutParameter.hh
./src/G4BuffercoutDestination.cc ./include/G4BuffercoutDestination.hh ./src/G4FilecoutDestination.cc ./include/G4FilecoutDestination.hh
./src/G4LockcoutDestination.cc ./include/G4LockcoutDestination.hh ./src/G4MasterForwardcoutDestination.cc ./include/G4MasterForwardcoutDestination.hh
 ./include/G4ios.hh  ./src/G4ios.cc  ./include/G4Types.hh ./src/G4UIparameter.cc ./include/G4UIparameter.hh
./src/G4Exception.cc ./include/G4Exception.hh ./src/G4Profiler.cc ./include/G4Profiler.hh
./src/Random.cc ./include/Randomize.hh


include/CLHEP/Evaluator/Evaluator.h
include/CLHEP/Evaluator/hash_map.icc
include/CLHEP/Evaluator/stack.icc
include/CLHEP/Evaluator/string.icc
include/CLHEP/Geometry/BasicVector3D.h
include/CLHEP/Geometry/Normal3D.h
include/CLHEP/Geometry/Plane3D.h
include/CLHEP/Geometry/Point3D.h
include/CLHEP/Geometry/Transform3D.h
include/CLHEP/Geometry/Vector3D.h
include/CLHEP/Random/ranluxpp/helpers.h
include/CLHEP/Random/ranluxpp/mulmod.h
include/CLHEP/Random/ranluxpp/ranlux_lcg.h
include/CLHEP/Random/DoubConv.h
include/CLHEP/Random/DualRand.h
include/CLHEP/Random/EngineFactory.h
include/CLHEP/Random/engineIDulong.h
include/CLHEP/Random/JamesRandom.h
include/CLHEP/Random/mixmax_skip_N8.icc
include/CLHEP/Random/mixmax_skip_N17.icc
include/CLHEP/Random/mixmax_skip_N240.icc
include/CLHEP/Random/MixMaxRng.h
include/CLHEP/Random/MTwistEngine.h
include/CLHEP/Random/NonRandomEngine.h
include/CLHEP/Random/RandBinomial.h
include/CLHEP/Random/RandBinomial.icc
include/CLHEP/Random/RandBit.h
include/CLHEP/Random/RandBreitWigner.h
include/CLHEP/Random/RandBreitWigner.icc
include/CLHEP/Random/RandChiSquare.h
include/CLHEP/Random/RandChiSquare.icc
include/CLHEP/Random/RandExponential.h
include/CLHEP/Random/RandExponential.icc
include/CLHEP/Random/RandFlat.h
include/CLHEP/Random/RandFlat.icc
include/CLHEP/Random/RandGamma.h
include/CLHEP/Random/RandGamma.icc
include/CLHEP/Random/RandGauss.h
include/CLHEP/Random/RandGauss.icc
include/CLHEP/Random/RandGaussQ.h
include/CLHEP/Random/RandGaussQ.icc
include/CLHEP/Random/RandGeneral.h
include/CLHEP/Random/RandGeneral.icc
include/CLHEP/Random/RandLandau.h
include/CLHEP/Random/RandLandau.icc
include/CLHEP/Random/RandPoisson.h
include/CLHEP/Random/RandPoisson.icc
include/CLHEP/Random/RandPoissonQ.h
include/CLHEP/Random/RandPoissonQ.icc
include/CLHEP/Random/RandStudentT.h
include/CLHEP/Random/RandStudentT.icc
include/CLHEP/Random/Random.h
include/CLHEP/Random/Random.icc
include/CLHEP/Random/RandomEngine.h
include/CLHEP/Random/RandomEngine.icc
include/CLHEP/Random/Randomize.h
include/CLHEP/Random/RanecuEngine.h
include/CLHEP/Random/Ranlux64Engine.h
include/CLHEP/Random/RanluxppEngine.h
include/CLHEP/Random/RanluxEngine.h
include/CLHEP/Random/RanshiEngine.h
include/CLHEP/Random/SeedTable.h
include/CLHEP/Random/Stat.h
include/CLHEP/Random/StaticRandomStates.h
include/CLHEP/Random/gaussQTables.cdat
include/CLHEP/Random/gaussTables.cdat
include/CLHEP/Random/poissonTables.cdat
include/CLHEP/Units/PhysicalConstants.h
include/CLHEP/Units/SystemOfUnits.h
include/CLHEP/Utility/defs.h
include/CLHEP/Utility/keywords.h
include/CLHEP/Utility/memory.h
include/CLHEP/Utility/noncopyable.h
include/CLHEP/Utility/use_atomic.h
include/CLHEP/Vector/AxisAngle.h
include/CLHEP/Vector/AxisAngle.icc
include/CLHEP/Vector/Boost.h
include/CLHEP/Vector/Boost.icc
include/CLHEP/Vector/BoostX.h
include/CLHEP/Vector/BoostX.icc
include/CLHEP/Vector/BoostY.h
include/CLHEP/Vector/BoostY.icc
include/CLHEP/Vector/BoostZ.h
include/CLHEP/Vector/BoostZ.icc
include/CLHEP/Vector/EulerAngles.h
include/CLHEP/Vector/EulerAngles.icc
include/CLHEP/Vector/LorentzRotation.h
include/CLHEP/Vector/LorentzRotation.icc
include/CLHEP/Vector/LorentzVector.h
include/CLHEP/Vector/LorentzVector.icc
include/CLHEP/Vector/Rotation.h
include/CLHEP/Vector/Rotation.icc
include/CLHEP/Vector/RotationInterfaces.h
include/CLHEP/Vector/RotationInterfaces.icc
include/CLHEP/Vector/RotationX.h
include/CLHEP/Vector/RotationX.icc
include/CLHEP/Vector/RotationY.h
include/CLHEP/Vector/RotationY.icc
include/CLHEP/Vector/RotationZ.h
include/CLHEP/Vector/RotationZ.icc
include/CLHEP/Vector/ThreeVector.h
include/CLHEP/Vector/ThreeVector.icc
include/CLHEP/Vector/TwoVector.h
include/CLHEP/Vector/TwoVector.icc

 src/AxisAngle.cc
  src/BasicVector3D.cc
  src/Boost.cc
  src/BoostX.cc
  src/BoostY.cc
  src/BoostZ.cc
  src/DoubConv.cc
  src/DualRand.cc
  src/EngineFactory.cc
  src/engineIDulong.cc
  src/erfQ.cc
  src/EulerAngles.cc
  src/Evaluator.cc
  src/flatToGaussian.cc
  src/gammln.cc
  src/JamesRandom.cc
  src/LorentzRotation.cc
  src/LorentzRotationC.cc
  src/LorentzRotationD.cc
  src/LorentzVectorB.cc
  src/LorentzVector.cc
  src/LorentzVectorC.cc
  src/LorentzVectorK.cc
  src/LorentzVectorL.cc
  src/LorentzVectorR.cc
  src/MixMaxRng.cc
  src/MTwistEngine.cc
  src/NonRandomEngine.cc
  src/Normal3D.cc
  src/Plane3D.cc
  src/Point3D.cc
  src/RandBinomial.cc
  src/RandBit.cc
  src/RandBreitWigner.cc
  src/RandChiSquare.cc
  src/RandExponential.cc
  src/RandFlat.cc
  src/RandGamma.cc
  src/RandGauss.cc
  src/RandGaussQ.cc
  src/RandGeneral.cc
  src/RandLandau.cc
  src/Random.cc
  src/RandomEngine.cc
  src/RandPoisson.cc
  src/RandPoissonQ.cc
  src/RandStudentT.cc
  src/RanecuEngine.cc
  src/Ranlux64Engine.cc
  src/RanluxppEngine.cc
  src/RanluxEngine.cc
  src/RanshiEngine.cc
  src/RotationA.cc
  src/Rotation.cc
  src/RotationC.cc
  src/RotationE.cc
  src/RotationInterfaces.cc
  src/RotationIO.cc
  src/RotationL.cc
  src/RotationP.cc
  src/RotationX.cc
  src/RotationY.cc
  src/RotationZ.cc
  src/setStdMath.cc
  src/setSystemOfUnits.cc
  src/SpaceVector.cc
  src/SpaceVectorD.cc
  src/SpaceVectorP.cc
  src/SpaceVectorR.cc
  src/StaticRandomStates.cc
  src/ThreeVector.cc
  src/ThreeVectorR.cc
  src/Transform3D.cc
  src/TwoVector.cc
  src/Vector3D.cc
  src/ZMinput.cc
)

target_compile_options(SiPM PRIVATE  --std=c++17 -pthread)


TARGET_LINK_LIBRARIES(SiPM z
       gsl gslcblas LibXml2::LibXml2)

install(TARGETS SiPM
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


