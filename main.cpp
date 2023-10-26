#include <iostream>

using namespace std;

#include"cnpy.h"
#include<complex>
#include<cstdlib>
#include<iostream>
#include<map>
#include<string>
#include "G4UIcmdWithADouble.hh"
#include"GateXMLDocument.hh"
#include "NumpyFile.hh"
#include <math.h>
#include <stdio.h>
#include <random>
#include "Randomize.hh"
#include <gsl/gsl_integration.h>
#include "SIPM/SiPM.hh"


int main(int argc, char *argv[])

{


    SiPM SiPM(argv[4]);




    SiPM.setSiPMFromXml (argv[1], argv[2]);//"hamamatsucross");
    SiPM.initializeSignal();
    SiPM.createPobsIntegrand();
    cnpy::NpyArray arr1 = cnpy::npy_load(argv[3]);//"XYZTime.npy");

    double* data1;
    data1 = arr1.data<double>();
    int y=arr1.shape[1];

    std::vector<SimuPulse> Simupulses;
    SimuPulse Pulse;


    for(int i = 0; i < y;i++){
        Pulse.X=data1[i];
        Pulse.Y=data1[i+1*y];
        Pulse.Z=data1[i+2*y];
        Pulse.time=data1[i+3*y];
        Pulse.ev=data1[i+4*y];

    Simupulses.insert(Simupulses.end(),Pulse);
    }

    SiPM.ProcessPulseList(Simupulses);

    return 0;
}

