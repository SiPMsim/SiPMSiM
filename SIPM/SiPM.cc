#include "SiPM.hh"
#include <typeinfo>
#include "G4UIcmdWithADouble.hh"
#include "G4UnitsTable.hh"
#include "NumpyFile.hh"
#include "GateXMLDocument.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "Randomize.hh"

#include <math.h>
#include <stdio.h>
#include <random>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <ostream>

uint64_t SiPMPulse::currentId = 0;



/*Create the files pulses.npy and signal.npy
Note that pulses.npy is not fully implemented yet*/
SiPM::SiPM(G4String result):
      m_histRes(0), m_startSignal(0), m_durationSignal(0),m_durationPulse(0),
      m_stepSignal(0), m_surface(NULL), m_deadTime(0.01), m_time_primary_value(0), m_t0(0.),
      m_computerSignal(NULL),m_SPTR(0)
{

  G4String res = result;
  m_file_Pulses.open(strcat(res.data(), "/pulses.npy"));

  //these are specific to GATE and not used here
  //m_file_Pulses.register_variable("uid", &m_pulseID);
  //m_file_Pulses.register_variable("runID", &m_runID);
  //m_file_Pulses.register_variable("eventID", &m_eventID);
  //m_file_Pulses.register_variable("trackID", &m_trackID);
  m_file_Pulses.register_variable("generalDetId", &m_generalDetId);

  m_file_Pulses.register_variable("time", &m_current_time_value);
  m_file_Pulses.register_variable("primary_time", &m_time_primary_value);
  m_file_Pulses.register_variable("value", &m_amplitude);
  m_file_Pulses.register_variable("type", &m_pulseType);
  m_file_Pulses.register_variable("amplitude", &m_current_pulse_amplitude);
  m_file_Pulses.register_variable("nbCrosstalks", &m_nbCrosstalks);
  m_file_Pulses.register_variable("crosstalk_rand", &m_crosstalk_rand);
  m_file_Pulses.register_variable("seed", &m_seed);

  G4String res1 = result;
  m_file.open(strcat(res1.data(),"/signal.npy"));
  m_file.register_variable("time", &m_current_time_value);
  m_file.register_variable("value", &m_current_pulse_value);
}



/*Deconstructor*/
SiPM::~SiPM()
{
  createSignal();
  saveAndPurge();
}


/* save the signal*/
void SiPM::saveAndPurge()
{
  m_file.writeHeader();
    for(auto &&signal_in_one_pixel: m_signalTable) {
      auto signal = &signal_in_one_pixel.signal;
      for (auto current_data_signal_iterator = signal->begin();
           current_data_signal_iterator != signal->end(); current_data_signal_iterator++) {
        auto current_data_signal = *current_data_signal_iterator;
        auto current_time_value = current_data_signal.time;
        m_current_time_value = current_time_value + this->getStartSignal();
        m_current_pulse_value = current_data_signal.value;
        m_file.fill();
      }
    }
    m_file.close();
}

/* process the impinging photons to the SiPM surface */
 void SiPM::ProcessPulseList(std::vector<SimuPulse> Simupulses)
{
   size_t n_pulses = Simupulses.size();
 if(n_pulses>0)
{
      G4int iterator =0;
      for(auto pulse_iterator = Simupulses.begin(); pulse_iterator != Simupulses.end(); pulse_iterator++)
      {
        iterator=iterator+1;
        SimuPulse pulse = *pulse_iterator;
        auto one_pulse = SiPMPulse();
        one_pulse.generalDetId= (G4int) m_generalDetId; // pulse->GetVolumeID().GetCopyNo(m_depth);
        ThreeVector v;
        v.ax=pulse.X;
        v.ay=pulse.Y;
        v.az=pulse.Z;
        one_pulse.microcell=wichMicroCell(v);

        if (one_pulse.microcell >= microCells.number){
            std::cout << v.ax << " " << v.ay << std::endl;
            std::cout <<  one_pulse.microcell << std::endl;
            std::cout << "Error, photon out of SiPM" << std::endl;
            }
        one_pulse.time= pulse.time;
        double sptr_noise = G4RandGauss::shoot(0, m_SPTR);
        one_pulse.time=one_pulse.time+sptr_noise;
        one_pulse.typePulse=primary;
        one_pulse.ev=pulse.ev;
        auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), one_pulse);
        this->m_SiPMPulsesTable.insert(it, one_pulse);
        m_time = pulse.time / 10e-12;
      }
    }
}



/* process the pulses by choosing the noise that will be applied as well as amplitude and if micro cell triggers*/
void SiPM::createSignal(){

    m_file_Pulses.writeHeader();
    //initialize SIPMs

    auto t_microCells_concerned = &((m_signalTable[0]).cellsList_ns);

    while (m_SiPMPulsesTable.size()>0)
    {
      auto& pulse1= m_SiPMPulsesTable.front();
      //specific to GATE not used here
      //m_pulse.trackID = pulse1.trackID;
      //m_pulse.typePulse=pulse1.typePulse;
      m_pulse.generalDetId=pulse1.generalDetId;
      m_generalDetId=pulse1.generalDetId;
      //m_pulse.runID = pulse1.runID;


      m_pulse.microcell=pulse1.microcell;
      m_pulse.ev=pulse1.ev;
      m_pulse.time= pulse1.time;
      m_SiPMPulsesTable.erase(m_SiPMPulsesTable.begin());
      m_nbCrosstalks = 0;
      m_crosstalk_rand = -1;


      p_signal_concerned = &(m_signalTable[m_pulse.generalDetId].signal);

      auto t_microCells_concerned = &((m_signalTable[m_pulse.generalDetId]).cellsList_ns);
      auto last_time=(t_microCells_concerned->at(m_pulse.microcell)).time;

      G4double ev= m_pulse.ev;
      G4double pde= 1;
      if (m_use_pde)
        pde= m_pde->compute(ev);
      //verfify if the microCells has begin to charge
      if ( m_pulse.time > last_time + m_deadTime) {
        //verifiy if triggering of the cell ( because of the recovery phase)
        G4double A = 1 - exp(-((m_pulse.time - (last_time + m_deadTime)) / m_tauRecovery_ns));

        if (m_pulse.typePulse != primary || G4UniformRand() < A*pde )
        {
            // a crosstalk does not do crosstalks
            if (m_pulse.typePulse != crosstalk)
            {
                generateCrosstalk(A);

            }
            if (true)
            {
                generateAfterpulses();
            }
            (t_microCells_concerned->at(m_pulse.microcell)).time = m_pulse.time;
            m_current_pulse_amplitude = writeSignal(last_time);


            if (m_pulse.typePulse == primary || m_pulse.typePulse == darkNoise )
                m_time_primary_value = m_pulse.time;

            //spefific to GATE not used here
            //m_runID = m_pulse.runID;
            //m_trackID = m_pulse.trackID;

            m_pulseType = m_pulse.typePulse;
            m_amplitude=A;
            m_current_time_value = m_pulse.time;
            m_pulseID = m_pulse.m_id;
            m_file_Pulses.fill();
        }
      }
    }

    m_file_Pulses.close();
}


G4double SiPM::writeSignal (G4double last_time)
{
    //compute pulse amplitude
    G4double amplitude = 0. ;
    if(  m_pulse.time - (last_time + m_deadTime ) < 0. )
      amplitude = 0.;
    else
      amplitude= m_signalAmplitude*(1-exp(-(m_pulse.time - (last_time + m_deadTime ) )/m_tauRecovery_ns)) ;
    double noise = G4RandGauss::shoot(0, m_signalAmplitudeSigma);
    amplitude = noise + amplitude ;

    //find location to write the signal
    GateSiPMDataSignal dump;
    dump.time= m_pulse.time;
    auto it1 = std::lower_bound(p_signal_concerned->begin(), p_signal_concerned->end(), dump);
    dump.time= m_pulse.time +m_durationPulse;
    auto it2 = std::lower_bound(p_signal_concerned->begin(), p_signal_concerned->end(), dump);

    //iterate in signal and write
    for (auto current_data_signal_iterator = it1; current_data_signal_iterator != it2; current_data_signal_iterator++)
    {
      auto compute_pulse_tmp = m_computerSignal->compute(current_data_signal_iterator->time - m_pulse.time) *  amplitude;
      (*current_data_signal_iterator).value += compute_pulse_tmp;
      (*current_data_signal_iterator).runID = m_pulse.runID ;
      (*current_data_signal_iterator).trackID = m_pulse.trackID ;
    }

    return amplitude;
}


// Generation of after pulses functions

double afterPulse_function(double t, void  *params){
  SiPM *sipm = (SiPM*) params;
  if (t- sipm->getT0() > 0) {
    G4double ret=sipm->getCap()*pow(t,sipm->geta())*exp(-t/sipm->getTauBuilk())*(1-exp(-(t-sipm->getT0())/sipm->getTauRecovery()));
    return double (ret);
  }
  else
    return 0.;
}

double afterCrosstalk_function(G4double t,void  *params){
  SiPM *sipm = (SiPM*) params;
  return sipm->getCct()*pow(t,sipm->getb())*exp(- t/sipm->getTauBuilk());
}


// This function allow to creat a probablity distribution of afterpulses it uses gsl
void SiPM::createPobsIntegrand(){
  G4double res1;
  G4double res2=0;
  G4double t=m_stepSignal;
  G4double resap, resct;
  double  error;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);

  gsl_function Ap;
  Ap.function = &afterPulse_function;
  Ap.params=this;

  gsl_function Act;
  Act.function = &afterCrosstalk_function;
  Act.params=this;

  SiPMProbs topush;
  topush.value=0.;
  topush.number=0.;
  m_afterPulse.push_back(topush);

  do
  {
    res1=res2;
    if (t>m_deadTime) {
      gsl_integration_qags(&Ap,m_deadTime,t,0,1e-8,100000, w,&resap,&error);
      gsl_integration_qags(&Act,m_histRes,t,0,1e-8,100000, w,&resct, &error);
    }
    else {
      resct=0;
      resap=0;
    }

    topush.value=resap;
    topush.number=t;
    m_afterPulse.push_back(topush);

    topush.value=resct;
    m_afterCrosstalk.push_back(topush);


    res2=resap+resct;
    t=t+m_stepSignal;
  }
  while (res2-res1 > 10e-8 || t <= (m_t0 + 2* m_stepSignal + m_deadTime ));
  gsl_integration_workspace_free (w);
}



// find the possition of crosstalk using the method of distance distribution
auto SiPM::crosstalk_pos_dist (){

  SiPMProbs pos; // a probaliblity
  G4double interpolated; // value resulting from interpolation
  G4double theta;
  G4int x;
  G4int y;
  G4int i;

  pos.value = G4UniformRand();  // fnum  ]0,1[
  auto microCell_concernedPos=microCellPos_k_to_i_j(m_pulse.microcell);

  //interpolation for the pixel distance

  if (pos.value > m_CrosstalkDispertion.at(0).value)
  {
    auto bornInf = std::lower_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
    bornInf--;
    auto bornSup = std::upper_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
    interpolated = bornInf->number + (bornSup->number - bornInf->number)*(pos.value - bornInf->value)/(bornSup->value - bornInf->value);
    interpolated=interpolated+1;
  }
  else
  {
    interpolated = 2;
  }
  // find cell for crosstalk and empiling the new pulse
  //        theta= rand()/(G4double)(RAND_MAX) * 2 * M_PI;
  theta = G4RandFlat::shoot(0.  ,2 * M_PI );
  x = interpolated * cos(theta) + microCell_concernedPos.first;
  y = interpolated * sin(theta) + microCell_concernedPos.second;

  return  std::make_tuple(x,y);
}



//find position of crosstalk based on map distribution
auto SiPM::crosstalk_pos_map (){

  struct res {    // Declare a local structure
    int i1, i2;

  };
  SiPMCrossProbs pos; // a probaliblity
  G4double theta;
  G4int x;
  G4int y;
  G4int i;
  std::vector<SiPMCrossProbs>::iterator born;

  pos.value = G4UniformRand();  // fnum  ]0,1[
  auto microCell_concernedPos=microCellPos_k_to_i_j(m_pulse.microcell);

  if (pos.value < m_Crossalk_map.back().value)
  {
    born= std::upper_bound(m_Crossalk_map.begin(), m_Crossalk_map.end(), pos);
    x =  microCell_concernedPos.first + born->X - (m_crosstalk_mapdim-1)/2 -1;
    y =  microCell_concernedPos.second + born->Y - (m_crosstalk_mapdim-1)/2 -1;

  }
  else
  {
    auto born1 = m_Crossalk_map.back();
    x =  microCell_concernedPos.first + born1.X - (m_crosstalk_mapdim-1)/2 -1;
    y =  microCell_concernedPos.second + born1.Y - (m_crosstalk_mapdim-1)/2 -1;
  }
  return std::make_tuple(x,y);
}


void SiPM::generateAfterpulses ( )
{
  SiPMProbs pos, posct; // a probaliblity
  G4double interpolated; // value resulting from interpolation
  G4double theta;
  G4int x;
  G4int y;
  G4int i;
  G4int h=0;
      

  pos.value = G4UniformRand();

  //After-pulses /////

  while (pos.value <= m_afterPulse.back().value){

    auto itap= std::lower_bound(m_afterPulse.begin(), m_afterPulse.end(), pos);
    auto time=itap->number;
    SiPMPulse new_pulse;
    new_pulse.microcell=m_pulse.microcell; 
    new_pulse.time= m_pulse.time + time ;
    new_pulse.typePulse=afterPulse;

    //Specific to GATE not used here
    //new_pulse.runID = m_pulse.runID;
    //new_pulse.trackID = m_pulse.trackID;
    new_pulse.generalDetId=m_generalDetId;
    auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), new_pulse);
    m_SiPMPulsesTable.insert(it, new_pulse);

    pos.value= pos.value + G4UniformRand();
  }


  //Delayed Crosstalks /////

  std::vector<std::pair<G4int, G4int> > fired_cells;
  posct.value = G4UniformRand();
  while (posct.value <= m_afterCrosstalk.back().value && h<=1000){
    //    G4cout << "generateAfterpulses creating aftercrosstalk  posct.value = " << posct.value << " m_afterCrosstalk.back().value = " <<  m_afterCrosstalk.back().value << G4endl;

    auto itap= std::lower_bound(m_afterCrosstalk.begin(), m_afterCrosstalk.end(), posct);
    auto time=itap->number;
    h++;

    // calculate microCells to fire for delayed crosstalk
    auto microCell_concernedPos=microCellPos_k_to_i_j(m_pulse.microcell);

    //        pos.value=rand()/(G4double)(RAND_MAX);

    if (m_Crossalk_map_bool)
    {
      std::tie(x, y)  = crosstalk_pos_map();
    }
    else
    {
      std::tie(x, y) = crosstalk_pos_dist();
    }

    auto position= std::make_pair(x,y);

    // if the microcell exists (inside the sipm) then we add a pulse
    if ( (x >= 0) & (y >= 0) && (x < microCells.numberDim1) && (y < microCells.numberDim2) && !(std::find(fired_cells.begin(), fired_cells.end(), position) != fired_cells.end()) )
    {
      SiPMPulse new_pulse;
      new_pulse.microcell=microCellPos_i_j_to_k(std::make_pair(x,y));
      new_pulse.time= m_pulse.time + time;
      new_pulse.runID = m_pulse.runID;
      new_pulse.trackID = m_pulse.trackID;
      new_pulse.typePulse=afterCrosstalk;
      new_pulse.generalDetId=m_generalDetId;
      fired_cells.push_back(position);

      auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), new_pulse);
      m_SiPMPulsesTable.insert(it, new_pulse);
      posct.value= posct.value + G4UniformRand();
    }
  }
}


void SiPM::generateCrosstalk (double A0)
{
  G4int i=0,y,x,h=0;
  SiPMProbs pos; // a probaliblity
  // calculate the number of crosstalks
  m_crosstalk_rand = pos.value = G4UniformRand();  // fnum  ]0,1[
  auto itnbCrosstalk= std::upper_bound(m_CrosstalkTable.begin(), m_CrosstalkTable.end(), pos);
  auto nbCrosstalk=itnbCrosstalk->number;
  m_nbCrosstalks = nbCrosstalk;

  // calculate microCells to fire for crosstalk

  std::vector<std::pair<G4int, G4int> > fired_cells;

  while (i < nbCrosstalk && h<=1000)
  {
    h++;
    //    pos.value=rand()/(G4double)(RAND_MAX);
    if (m_Crossalk_map_bool)
    {
      std::tie(x, y)  = crosstalk_pos_map();
    }
    else
    {
      std::tie(x, y) = crosstalk_pos_dist();
    }
    auto position= std::make_pair(x,y);
    // if the microcell exists (inside the sipm) then we add a pulse
    if ( (x >= 0) && (y >= 0) && (x < microCells.numberDim1) && (y < microCells.numberDim2) && !(std::find(fired_cells.begin(), fired_cells.end(), position) != fired_cells.end()))
    {
      SiPMPulse new_pulse;
      new_pulse.microcell=microCellPos_i_j_to_k( std::make_pair(x,y));
      double sptr_noise = G4RandGauss::shoot(0, m_SPTR*sqrt(i+2));
      while ((m_SPTR*2.355*sqrt(i+1) + sptr_noise) < 0){
            sptr_noise = G4RandGauss::shoot(0, m_SPTR*sqrt(i+2));}
      new_pulse.time= m_pulse.time+abs(m_SPTR*2.355*sqrt(i+1)) + sptr_noise ;
      new_pulse.runID = m_pulse.runID;
      new_pulse.trackID = m_pulse.trackID;
      new_pulse.typePulse=crosstalk;
      new_pulse.generalDetId=m_generalDetId;
      //G4cout << " les microcells  " << new_pulse.microcell << " " << pulse.microcell << G4endl;
      (m_SiPMPulsesTable).insert((m_SiPMPulsesTable).begin(), new_pulse);  
      fired_cells.push_back(position);
      i++;
    }
  }
}


void SiPM::initializeSignal()
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine());
  G4long seed;

  if (m_seed_bool)
  {
    seed =  m_seed;
  }
  else
  {
    seed = time(NULL);
  }
  G4Random::setTheSeed(seed);
  G4RandGauss::setTheSeed(seed);
  G4RandGeneral::setTheSeed(seed);

  //creat nested volume
  SignalInOnePixel signal_in_one_pixel;
  m_signalTable.push_back(signal_in_one_pixel);
  setMircoCellsProperties();
  m_darkNoiseProb =  (m_stepSignal/1e9) * m_darkNoise;
  G4double durationSignal = this->getDurationSignal();
  GateSiPMCell tocreate; // micro_cells list ( for creation)

  for(size_t i = 0; i < m_signalTable.size(); i++)
  {
    G4double prob;

    //Get the good sipm concerned
    p_signal_concerned = &((m_signalTable[i]).signal);

    m_signalTable[i].generalDetId=i;

    G4double startSignal = -m_durationPulse -3 * m_tauRecovery_ns;

    //creat microcells in SiPM
    auto t_microCells_concerned = &((m_signalTable[i]).cellsList_ns);
    if (t_microCells_concerned->size()<=1)
    {
      tocreate.time=startSignal;
      *t_microCells_concerned= *new std::vector<GateSiPMCell> (microCells.number,tocreate);
    }
    if (p_signal_concerned->size() == 0)
    {
      uint64_t reserve = (uint64_t) (durationSignal-startSignal) / this->getStepSignal();
      //      if(nVerboseLevel > 1)
      //        G4cout << "[GateSiPM::ProcessPulseList]" << "Allocate signal, I reserve " << reserve << " elements" << Gateendl;
      p_signal_concerned->reserve(reserve);

      for (G4double current_time_value = startSignal; current_time_value < this->getDurationSignal(); current_time_value += this->getStepSignal()   )
      {

            GateSiPMDataSignal d;
            d.time = current_time_value;

            d.value = G4RandGauss::shoot(0, m_whiteNoiseSigma);
            //                d.value = 0.;
            d.runID = 0;
            p_signal_concerned->push_back(d);
      }


      G4double timing=0;
      G4double current_time_value;
      G4double cpt=0;
      while (timing < this->getDurationSignal())
      {
            prob = G4UniformRand()/m_darkNoiseProb * 2;  // fnum  ]0,1[
            timing= timing + prob;

            if (timing <= this->getDurationSignal())
            {
                cpt=cpt+1;
                current_time_value=timing;
                SiPMPulse one_pulse;
                one_pulse.microcell = G4RandFlat::shootInt((long)0,  microCells.number - 1);
                one_pulse.time = current_time_value + G4RandFlat::shoot(0., this->getStepSignal());
                one_pulse.runID = -1;
                one_pulse.trackID = -1;
                one_pulse.typePulse=darkNoise;
                one_pulse.generalDetId=i;
                auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), one_pulse);
                if(m_SiPMPulsesTable.size() > 0) //iteraror can be deferenced
                {
                    one_pulse.runID  = it->runID;
                }

                m_SiPMPulsesTable.insert(it, one_pulse);


            }
      }
    }
  }

}


// Attaches a G4MaterialPropertiesTable to the optical surface.

void SiPM::DescribeMyself(size_t indent)
{
  std::cout <<  "GateSiPM " << std::endl;
  std::cout << "TauRise: " << getTauRise()<< std::endl;
  std::cout << "TauFall: " << getTauFall() << std::endl;
  std::cout << "Start signal: " << getStartSignal() << std::endl;
  std::cout << "Duration signal: " << getDurationSignal() << std::endl;
  std::cout << "Step signal: " << getStepSignal() << std::endl;


}

void  SiPM::setMircoCellsProperties()
{
  microCells.numberDim1= m_mircoCellsPitch[m_surface[0]]/ microCells.lengthDim1;
  microCells.numberDim2= m_mircoCellsPitch[m_surface[1]]/ microCells.lengthDim2;
  microCells.number= (G4int) microCells.numberDim1 * microCells.numberDim2;
}


void SiPM::setCellsDimentions(G4double X, G4double Y, G4double Z)
{
  m_mircoCellsPitch=new G4double[3] {X, Y, Z};
}

void SiPM::setMircoCellsDim(G4double microCellsDim1,G4double microCellsDim2 )
{
  microCells.lengthDim1= microCellsDim1;
  microCells.lengthDim2= microCellsDim2;
}

G4int SiPM::wichMicroCell(ThreeVector localPos) const
{
  G4double dump [3]={localPos.ax, localPos.ay,localPos.az};
  G4int i=  (dump[m_surface[0]] + (m_mircoCellsPitch[m_surface[0]]/2)) / (microCells.lengthDim1);
  G4int j=  (dump[m_surface[1]] + (m_mircoCellsPitch[m_surface[1]] /2) ) / (microCells.lengthDim2);

  return i*microCells.numberDim2 + j;

}

G4int SiPM::microCellPos_i_j_to_k (const std::pair<G4int, G4int> &pos ){
  return pos.first*microCells.numberDim2 + pos.second;
}

std::pair<G4int, G4int> SiPM::microCellPos_k_to_i_j(G4int pos)
{
  G4int i= pos/microCells.numberDim2;
  G4int j=pos % microCells.numberDim2;
  return std::make_pair(i, j);
}

G4double SiPM::getTauRise() const
{
  return m_tauRise_ns;
}

void SiPM::setTauRise(G4double tauRise)
{
  m_tauRise_ns = tauRise / 1e-9;
}

G4double SiPM::getTauFall() const
{
  return m_tauFall_ns;
}

void SiPM::setTauFall(G4double tauFall)
{
  m_tauFall_ns = tauFall / 1e-9;
}

G4double SiPM::getStartSignal() const
{
  return m_startSignal;
}

void SiPM::setStartSignal(G4double startSignal)
{
  m_startSignal = startSignal;
}

G4double SiPM::getDurationSignal() const
{
  return m_durationSignal;
}

void SiPM::setDurationSignal(G4double endSignal)
{
  m_durationSignal = endSignal;
}


G4double SiPM::getDurationPulse() const
{
  return m_durationPulse;
}

void SiPM::setDurationPulse(G4double duration)
{
  m_durationPulse = duration;
}


G4double SiPM::getSPTR() const
{
  return m_SPTR;
}

void SiPM::setSPTR(G4double SPTR)
{
  m_SPTR= SPTR;
}


void SiPM::setStepSignal(G4double stepSignal)
{
  m_stepSignal = stepSignal;
}


const G4int* SiPM::getSurface() const
{
  return m_surface;
}

G4double* SiPM::getCellsDimentions() const
{
  return m_mircoCellsPitch;
}

G4double SiPM::getStepSignal() const
{
  return m_stepSignal;
}
G4double SiPM::getTauRecovery() const {
  return m_tauRecovery_ns;
}

void SiPM::setTauRecovery(G4double tauRecovery){
  m_tauRecovery_ns=tauRecovery;
}


void SiPM::setSurface(G4String surface){

  if ( surface=="ZY"){
    m_surface=new G4int [2] {2,1};
  }
  else if (surface=="XZ"){
    m_surface=new G4int [2] {0,2};
  }
  else if (surface=="ZX"){
    m_surface=new G4int [2] {2,0};
  }
  else if (surface=="YZ"){
    m_surface=new G4int [2] {1,2};
  }
  else if (surface=="XY"){
    m_surface=new G4int [2] {0,1};
  }
  else  {
    m_surface=new G4int [2] {1,0};
  }
}


void SiPM::setSiPMFromXml  (G4String location, G4String type){
  GateXMLDocument* doc = new GateXMLDocument( location + "/SiPM.xml");
  if (doc->Ok())
  {
    doc->Enter();
    if (doc->Find("sipm",type))
    {
      doc->Enter();
      doc->First();
      if (doc->Find("propertiestable"))
      {   G4double value;
        G4double unit;
        doc->Enter();
        while (doc->Next())
        {
          if (doc->GetName() == "property")
          {

            G4String property = doc->GetProperty("name");
            G4String valuestr = doc->GetProperty("value");
            value = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());

            if (property == "surface"){
              this->setSurface(valuestr);
            }
            else {
                G4String unitstr = "1 " + doc->GetProperty("unit");
                if (property == "tauFall"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());

                  if(!m_computerSignal)
                    m_computerSignal = new ComputeSignalFormula();
                  static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauFall_ns = value*unit;
                }
                if (property == "tauRise"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  //this->setTauRise(value*unit);
                  if(!m_computerSignal)
                    m_computerSignal = new ComputeSignalFormula();
                  static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauRise_ns = value*unit;
                }
                if (property == "tau" & !m_computerSignal)
                {

                      m_computerSignal = new ComputeSignalFormula();
                      G4String unitstr = "1 " + doc->GetProperty("unit");
                      unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                      G4String valuestr  = doc->GetProperty("taurise");
                      G4double value1    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                      static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauRise_ns = value1*unit;

                      G4String valuestr2  = doc->GetProperty("taufall");
                      G4double value2 = G4UIcmdWithADouble::GetNewDoubleValue(valuestr2.c_str());
                      static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauFall_ns = value2*unit;

                }

                if (property == "circuit" & !m_computerSignal)
                {

                      m_computerSignal = new ComputeSignalFormula();
                      G4String valuestr  = doc->GetProperty("Rs");
                      G4double Rs    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                      valuestr  = doc->GetProperty("Rq");
                      G4double Rq    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                      valuestr  = doc->GetProperty("Cd");
                      G4double Cd    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                      valuestr  = doc->GetProperty("Cq");
                      G4double Cq    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                      valuestr  = doc->GetProperty("Cg");
                      G4double Cg    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());

                      static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauFall_ns = Rs * (Cq + Cg);
                      static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauRise_ns = Rq * (Cd + Cq);
                }
                if (property == "durationSignal"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  this->setDurationSignal(value*unit);
                }
                if (property == "stepSignal"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  this->setStepSignal(value*unit);
                }
                if (property == "histRes"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_histRes=(value*unit);
                }
                if (property == "durationPulse"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  this->setDurationPulse(value*unit);
                }
                if (property == "whiteNoiseSigma"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_whiteNoiseSigma=value*unit*10e5;
                }
                if (property == "seed"){
                  m_seed=value;
                  m_seed_bool=true;
                }
                if (property == "tauRecovery"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  this->setTauRecovery(value*unit);
                }
                if (property == "SPTR"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  this->setSPTR(value*unit);
                }
                if (property == "signalDeconvolvedAmplitude"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_signalAmplitude=value*unit*10e5;
                }
                if (property == "signalDeconvolvedAmplitudeSigma"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_signalAmplitudeSigma=value*unit*10e5;
                }
                if (property == "darkNoise"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_darkNoise=value*unit*10e8;
                }
                if (property == "tauBuilk"){
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_tauBuilk=(value*unit);
                }
                if (property == "Cap")
                  m_Cap=value;
                if (property == "t0")
                {
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_t0=value*unit;
                }
                if (property == "deadTime")
                {
                  unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
                  m_deadTime=value*unit;
                }
                if (property == "Cct")
                  m_Cct=value;
                if (property == "a")
                  m_a=value;
                if (property == "b")
                  m_b=value;
              }
          }
          else if (doc->GetName() == "propertyvector")
          {
            G4String property = doc->GetProperty("name");
            if (property == "CROSSTALK")
            {
              // read vector
              G4double accumulate=0;
              doc->Enter();
              G4int iterator=0;
              SiPMProbs topush;
              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("value");
                G4double value     = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                accumulate+=value;
                topush.value=accumulate;
                topush.number=iterator;
                m_CrosstalkTable.push_back(topush);
                iterator ++;
              }
            }

            else if (property == "PDE")
            {
              // read vector
              m_use_pde=true;
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();

              m_pde = new PDE();
              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("energy");
                G4double time = (G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str()))*unit;
                //G4cout << time << G4endl;
                valuestr  = doc->GetProperty("value");
                G4double value = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                static_cast<PDE*>(m_pde)->add_pair(time, value);
               //G4cout << value << " " << time<< G4endl;
              }
                sortVectors(m_pde->ev, std::less<G4double>(), m_pde->ev, m_pde->pde );
            }


            else if (property == "PULSE" & !m_computerSignal)
            {
              // read vector
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();

              m_computerSignal = new ComputeSignalVector();

              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("time");
                G4double time = (G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str()))*unit;
                valuestr  = doc->GetProperty("value");
                G4double value = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                static_cast<ComputeSignalVector*>(m_computerSignal)->add_pair(time, value);               //G4cout << value << " " << time<< G4endl;
              }

            }
            else if (property == "DIMENTIONS")
            {
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();
              doc->Find("ve");
              G4String valuestr  = doc->GetProperty("value");
              G4double value1    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              valuestr  = doc->GetProperty("value");
              G4double value2 = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              this->setMircoCellsDim(value1*unit,value2*unit);
            }
            else if (property == "SiPM DIMENTIONS")
            {
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();
              doc->Find("ve");
              G4String valuestr  = doc->GetProperty("value");
              G4double value1    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              valuestr  = doc->GetProperty("value");
              G4double value2 = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              valuestr  = doc->GetProperty("value");
              G4double value3 = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              this->setCellsDimentions(value1*unit,value2*unit, value3*unit);
            }
            else if (property == "CROSSTALK_MAP")
            {
              G4String unitstr = doc->GetProperty("dim");
               m_crosstalk_mapdim    = G4UIcmdWithADouble::GetNewDoubleValue(unitstr.c_str());
               m_Crossalk_map_bool = true;

              // read vector
              G4double accumulate=0;
              doc->Enter();
              SiPMCrossProbs topush;
              while (doc->Find("ve"))
              {
                G4String xstr  = doc->GetProperty("X");
                G4double X     = G4UIcmdWithADouble::GetNewDoubleValue(xstr.c_str());
                G4String ystr  = doc->GetProperty("Y");
                G4double Y     = G4UIcmdWithADouble::GetNewDoubleValue(ystr.c_str());
                G4String valuestr  = doc->GetProperty("value");
                G4double value     = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                accumulate+=value;
                topush.value=accumulate;
                topush.X=X;
                topush.Y=Y;
                m_Crossalk_map.push_back(topush);
              }
            }
            else if (property == "CROSSTALK_DISPERTION")
            {
              // read vector
              G4double accumulate=0;
              doc->Enter();
              G4int iterator=0;
              SiPMProbs topush;
              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("value");
                G4double value     = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                accumulate+=value;
                topush.value=accumulate;
                topush.number=iterator;
                m_CrosstalkDispertion.push_back(topush);
                iterator ++;
              }
            }
            doc->Leave();
          }
        }
        doc->Leave();
      }
    }
  }
}

G4double SiPM::interp1(G4double x_new)
{
  G4double dx, dy, m, b;
  size_t x_max_idx = m_Pulse_time.size() - 1;

  auto it = std::lower_bound(m_Pulse_time.begin(), m_Pulse_time.end(), x_new);
  std::size_t idx = std::distance(m_Pulse_time.begin(), it);

  if (m_Pulse_time[idx] > x_new)
  {
    dx = idx > 0 ? (m_Pulse_time[idx] - m_Pulse_time[idx - 1]) : (m_Pulse_time[idx + 1] - m_Pulse_time[idx]);
    dy = idx > 0 ? (m_Pulse_value[idx] - m_Pulse_value[idx - 1]) : (m_Pulse_value[idx + 1] - m_Pulse_value[idx]);
  }
  else
  {
    dx = idx < x_max_idx ? (m_Pulse_time[idx + 1] - m_Pulse_time[idx]) : (m_Pulse_time[idx] - m_Pulse_time[idx - 1]);
    dy = idx < x_max_idx ? (m_Pulse_value[idx + 1] - m_Pulse_value[idx]) : (m_Pulse_value[idx] - m_Pulse_value[idx - 1]);
  }
  m = dy / dx;
  b = m_Pulse_value[idx] - m_Pulse_time[idx] * m;
  return x_new * m + b;

}


