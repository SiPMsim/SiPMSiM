#ifndef SIPM_H
#define SIPM_H
#include <vector>
#include <gsl/gsl_integration.h>
#include "G4OpticalSurface.hh"
#include "NumpyFile.hh"
#include "GateXMLDocument.hh"
#include <ThreeVector.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iomanip>
#include <numeric>



template <typename T, typename Compare>
void getSortPermutation(
    std::vector<unsigned>& out,
    const std::vector<T>& v,
    Compare compare = std::less<T>())
{
    out.resize(v.size());
    std::iota(out.begin(), out.end(), 0);

    std::sort(out.begin(), out.end(),
              [&](unsigned i, unsigned j){ return compare(v[i], v[j]); });
}

template <typename T>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t)
{
    assert(order.size() == t.size());
    std::vector<T> st(t.size());
    for(unsigned i=0; i<t.size(); i++)
    {
        st[i] = t[order[i]];
    }
    t = st;
}

template <typename T, typename... S>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t,
    std::vector<S>&... s)
{
    applyPermutation(order, t);
    applyPermutation(order, s...);
}

template<typename T, typename Compare, typename... SS>
void sortVectors(
    const std::vector<T>& t,
    Compare comp,
    std::vector<SS>&... ss)
{
    std::vector<unsigned> order;
    getSortPermutation(order, t, comp);
    applyPermutation(order, ss...);
}

// make less verbose for the usual ascending order
template<typename T, typename... SS>
void sortVectorsAscending(
    const std::vector<T>& t,
    std::vector<SS>&... ss)
{
    sortVectors(t, std::less<T>(), ss...);
}

struct ThreeVector 	{
        G4double  	ax;
        G4double  	ay;
        G4double  	az;
    };

class ComputeSignal{

public:
    virtual G4double compute(G4double t) = 0;

};

class ComputeSignalFormula: public ComputeSignal{

    public:
        ComputeSignalFormula()
        {
            m_tauRise_ns = 0.;
            m_tauFall_ns = 0.;
        }


    public:
        G4double compute(G4double t) override {
            return (1-exp(-t/m_tauRise_ns))*exp(-t/m_tauFall_ns);
        }


    public:
        G4double m_tauFall_ns;
        G4double m_tauRise_ns;

};





class ComputeSignalVector: public ComputeSignal{
    public:
        std::vector<G4double> m_x;
        std::vector<G4double> m_y;

public:
    ComputeSignalVector()
    {

    }

    void add_pair(G4double x, G4double y)
    {
        m_x.push_back(x);
        m_y.push_back(y);
    }

    G4double compute(G4double x_new) override {
        G4double dx, dy, m, b;
        size_t x_max_idx = m_x.size() - 1;

      //  size_t idx = nearestNeighbourIndex(x, x_new);
        auto it = std::lower_bound(m_x.begin(), m_x.end(), x_new);

        if (it==m_x.end())
            return 0;

        std::size_t idx = std::distance(m_x.begin(), it);

        if (m_x[idx] > x_new)
        {
          dx = idx > 0 ? (m_x[idx] - m_x[idx - 1]) : (m_x[idx + 1] - m_x[idx]);
          dy = idx > 0 ? (m_y[idx] - m_y[idx - 1]) : (m_y[idx + 1] - m_y[idx]);
        }
        else
        {
          dx = idx < x_max_idx ? (m_x[idx + 1] - m_x[idx]) : (m_x[idx] - m_x[idx - 1]);
          dy = idx < x_max_idx ? (m_y[idx + 1] - m_y[idx]) : (m_y[idx] - m_y[idx - 1]);
        }
        m = dy / dx;
        b = m_y[idx] - m_x[idx] * m;
        return x_new * m + b;
    }




};


class GateArrayParamsFinder;

/* structure for the signal */
struct GateSiPMDataSignal
{
  G4double time;
  G4double value;
  G4int runID;
  G4int trackID;


  bool operator<(const GateSiPMDataSignal & rhs) const
     {
         return (time < rhs.time);
     }

};

typedef std::vector<GateSiPMDataSignal> GateSiPMSignal;

struct GateSiPMCell
{
    G4double time;
};

class SignalInOnePixel
{
public:
  G4int generalDetId;
  GateSiPMSignal signal;
  std::vector<GateSiPMCell> cellsList_ns;
};


/* structure of microcells information*/
struct GateSiPMMicrocells
{
  G4int number;
  G4int numberDim1;
  G4int numberDim2;
  G4double lengthDim1;
  G4double lengthDim2;
};

/*Structure of probabilities*/

struct SiPMProbs
{
    G4double value;
    G4double number;

    bool operator<(const SiPMProbs & rhs) const
       {
           return (value < rhs.value);
       }
};

struct SiPMCrossProbs
{
    G4double X;
    G4double Y;
    G4double value;

    bool operator<(const SiPMCrossProbs & rhs) const
       {
           return (value < rhs.value);
       }
};

enum typeOfPulse {primary = 1, darkNoise, crosstalk, afterPulse, afterCrosstalk, };


class SimuPulse
{
public:
  G4double X;
  G4double Y;
  G4double Z;
  G4double time;
  G4double ev;

};



class SiPMPulse
{
public:
  SiPMPulse()
  {

    m_id = currentId;
    currentId ++;
  }

  G4int microcell;
  G4double time;
  G4int runID;
  G4int trackID;
  typeOfPulse typePulse;
  G4int SiPM_number;
  G4int micocell_number;
  G4int generalDetId;
  G4double ev;

  uint64_t m_id;

  static uint64_t currentId;

  bool operator<(const SiPMPulse & rhs) const
     {
         return (time < rhs.time);
     }
};




class PDE
{
  public:
      std::vector<G4double> ev;
      std::vector<G4double> pde;




  void add_pair(G4double x, G4double y)
  {
     ev.push_back(x);
     pde.push_back(y);
  }

  G4double compute( G4double x_new){
      G4double dx, dy, m, b;
      auto it = std::lower_bound(ev.begin(), ev.end(), x_new);
      //for (G4double i: ev)

      std::size_t idx = std::distance(ev.begin(), it);

      if (idx==0)
      {

        dx = (ev[idx + 1] - ev[idx]);
        dy = (pde[idx + 1] - pde[idx]);
      }
      else if (ev[idx-1] < x_new)
      {   idx=idx-1;
          dx = ev[idx] - ev[idx-1] ;
          dy = pde[idx] - pde[idx-1] ;
      }
      else
      { idx=idx-1;
        dx = ev[idx + 1] - ev[idx] ;
        dy = pde[idx + 1] - pde[idx] ;
      }

      m = dy / dx;
      b = pde[idx] - ev[idx] * m;

      if ((x_new * m + b) <0)
          return 0;

      return x_new * m + b;
  }

 };


class SiPM
{

public:


  G4double m_current_time_value;
  G4double m_current_pulse_value;
  G4double m_current_pulse_amplitude;
  G4int m_trackID;
  G4int m_runID;
  G4int m_stepID;
  G4int m_parentID;
  G4int m_eventID;
  G4double m_time;
  G4double m_energy;
  G4int  m_PDGEncoding;
  G4String m_topVolume;
  G4String m_bottomVolume;
  G4int m_nbCrosstalks;
  G4double m_crosstalk_rand;
  G4double m_histRes;
  SiPMPulse m_pulse;
  GateSiPMSignal* p_signal_concerned;
  PDE *m_pde;
  G4bool m_use_pde= false;
  G4bool m_seed_bool=false;
  G4double m_tauRise_ns;
  G4double m_tauFall_ns;
  G4double m_startSignal;
  G4double m_deadTime;
  G4double m_durationSignal;
  G4double m_durationPulse;
  G4double m_stepSignal;
  G4double m_microCells;
  G4String m_volume;
  G4int *m_surface;   // indicate the surface of detection
  G4double *m_mircoCellsPitch; //dimension of cells
  GateSiPMMicrocells microCells; // defin microcells proprieties
  G4double m_tauRecovery_ns;
  G4int *m_posCellConcerned;
  G4double m_whiteNoiseSigma;
  G4double m_signalAmplitude;
  G4double m_signalAmplitudeSigma;
  G4double m_darkNoise;
  G4double m_darkNoiseProb;
  G4double m_tauBuilk;
  G4double m_Cap;
  G4double m_Cct;
  G4double m_a;
  G4double m_b;
  G4double m_seed;

  G4double m_t0;
  G4double m_SPTR;
  uint64_t m_pulseID;

  // time and value of pulse signal
  std::vector<G4double> m_Pulse_time;
  std::vector<G4double> m_Pulse_value;

  std::vector<SignalInOnePixel> m_signalTable;
  std::vector<SiPMPulse> m_SiPMPulsesTable;

  std::vector<SiPMProbs> m_CrosstalkTable;
  std::vector<SiPMProbs> m_CrosstalkDispertion;

  std::vector<SiPMProbs> m_afterPulse;
  std::vector<SiPMProbs> m_afterCrosstalk;

  NumpyFile m_file;
  NumpyFile m_file_Pulses;

  G4int m_level1_debug, m_level2_debug;

  G4int m_generalDetId;

  GateArrayParamsFinder *m_arrayFinder;
  size_t m_nbX, m_nbY, m_nbZ; //!< Parameters of the matrix of detection
  size_t m_i, m_j, m_k;  //!< position \i,\j,\k in the matrix

  std::vector<G4int> m_numberOfComponentForLevel;
  G4int  m_numberOfHigherLevels;

  std::vector<SiPMCrossProbs> m_Crossalk_map;
  G4int m_crosstalk_mapdim;
  G4bool m_Crossalk_map_bool=false;

  G4double m_time_primary_value;
  G4double m_amplitude;
  G4int m_pulseType;

  ComputeSignal *m_computerSignal;


  SiPM(G4String result);
  virtual ~SiPM();

  void ProcessPulseList(std::vector<SimuPulse> Simupulses);

  void DescribeMyself(size_t indent);

  G4double interp1(G4double x_new );

  void generateAfterpulses ();

  void createPobsIntegrand();

  void initializeSignal();

  // initialize parameters of SiPM

  G4double writeSignal ( G4double last_time);

  // function to compute crosstalk ( add new pulses in the list)
  G4int* crosstalkspos();

  G4int microCellPos_i_j_to_k (const std::pair<G4int, G4int> &pos);

  void generateCrosstalk (double A0);
  /* microCell pos Geometrically */
  std::pair<G4int, G4int> microCellPos_k_to_i_j(G4int pos);

  /* Create the signal corresponding to pulses (called at deconstructor)*/
  void createSignal();

  /* Save and purge the signal with a time lower that the new pulse Pulselist*/
  void saveAndPurge();

  /* microCell position in the list*/
  G4int wichMicroCell(ThreeVector localPos) const;

  // Compute number of cells per side and their dimentions. note that this function is called
  // in setSiPMFromXml ot in setVolume depending on how the mac file is filled
  void setMircoCellsProperties();

  G4double* getCellsDimentions() const;
  void setCellsDimentions(G4double X, G4double Y, G4double Z);

  /*This function allows to complete the SiPM properties with the XML file named SiPM.xml
  Note that the GateXMLDocument class is used. You can find more details about in source/general/
  GateXMLDocument is only instancied to read tables with energies and values. Do it if you need to read something different
  */
  void setSiPMFromXml ( G4String location, G4String type);

  auto crosstalk_pos_dist ();
  auto crosstalk_pos_map ();


 /* Get the surface of detection of the SiPM ( XY, YZ or XZ) */
  const G4int* getSurface() const;
  void setSurface(G4String surface);

 /* Set the pitch of microcells */
  void setMircoCellsDim(G4double microCellsDim1_micron, G4double microCellsDim2_micron );

  G4double getTauFall() const;
  void setTauFall(G4double tauFall);

  G4double getTauRecovery() const;
  void setTauRecovery(G4double tauRecovery);

  G4double getSPTR() const;
  void setSPTR(G4double SPTR);


  G4double getTauRise() const;
  void setTauRise(G4double tauRise);

  G4double getStartSignal() const;
  void setStartSignal(G4double startSignal);

  G4double getDurationSignal() const;
  void setDurationSignal(G4double endSignal);

  G4double getDurationPulse() const;
  void setDurationPulse(G4double duration);
  void setDurationPulSignalInOnePixelse(G4double endSignal);

  G4double getStepSignal() const;
  void setStepSignal(G4double stepSignal);


  G4double getTauRecovery()  {
      return this->m_tauRecovery_ns;
  }

  G4double getTauBuilk() const {
      return this->m_tauBuilk;
  }

  G4double getCap() const {
      return this->m_Cap;
  }

  G4double getCct() const {
      return this->m_Cct;
  }

  G4double geta(){
      return this->m_a;
  }

  G4double getb() const {
      return this->m_b;
  }

  G4double getT0() const {
      return m_t0;
  }

};


#endif // SIPM_H
