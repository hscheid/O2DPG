
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include "Generators/GeneratorPythia8.h"
#include "TRandom3.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "CCDB/BasicCCDBManager.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"

#include <map>
#include <unordered_set>

class GeneratorPythia8ThermalRadiation : public o2::eventgen::GeneratorPythia8
{
public:
  /// Constructor
  GeneratorPythia8ThermalRadiation()
  {
    // Load theory file
    loadTheoryFileFromCcdb(); // set member histo
    // init member variable
    minRap = -1.;
    maxRap = 1.;
    mDaughters = {0.00051099895, 0.00051099895};
    nGammaPerEvent = theoryPrediction_th2d->Integral("width");
  }

  ///  Destructor
  ~GeneratorPythia8ThermalRadiation() = default;

protected:
  //__________________________________________________________________
  Bool_t generateEvent() override
  {

    /// reset event
    mPythia.event.clear();
    // not sure if I should not somehow add a vertex, but it should eb the vertex of the events I use as UE. Not usre how to do in cocktail generator
    xProd = 0.0;
    yProd = 0.0;
    zProd = 0.0;
    // How many virtual photons are there in the event?
    nGamma = rand.Poisson(nGammaPerEvent);
    // loop to generate thermal photon/ dielectron pair
    for (size_t iGamma = 0; iGamma < nGamma; iGamma++) {
      generateThermalPair();
      sign = 1;
      if (rand.Uniform(0., 1.) < 0.5) sign *= -1;
      Pythia8::Particle lAddedParticle1 = createParticle(sign, d1->Px(), d1->Py(), d1->Pz(), d1->E(), d1->M(), 0.0, 0.0, 0.0);
      Pythia8::Particle lAddedParticle2 = createParticle(sign * -1, d2->Px(), d2->Py(), d2->Pz(), d2->E(), d2->M(), 0.0, 0.0, 0.0);
      mPythia.event.append(lAddedParticle1);
      mPythia.event.append(lAddedParticle2);
    }
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    /// go to next Pythia event
    mPythia.next();

    return true;
  }

  int loadTheoryFileFromCcdb()
  {
    o2::ccdb::CcdbApi ccdb_api;
    // TFile* theoryFile;
    // TH2D* theoryPrediction_th2d;
    std::string theoryName = "thermalRadiation_Rapp"; // this is used to pick the calculation in the file
    std::string pathToFile = "Users/h/hscheid/theoryFile/";
    std::string theoryFileName = "lmee_theoryFile.root";

    ccdb_api.init("https://alice-ccdb.cern.ch");
    std::map<string, string> metadataRCT;
    std::string tmpDir = ".";

    // get the theory file from CCDB
    if (!ccdb_api.retrieveBlob(pathToFile.data(), tmpDir, metadataRCT, -1, false, theoryFileName.data())) {
      LOG(fatal) << "Input file for theory prediction not found on CCDB, please check the pathInputFile and nameInputFile!";
    }
    /// use input correction file from local filesystem
    std::unique_ptr<TFile> theoryFile(TFile::Open(("./" + theoryFileName).c_str(), "READ"));

    if (!theoryFile.get()) {
      LOG(fatal) << "Something wrong with the input file. Fix it!";
    }

    theoryPrediction_th2d = dynamic_cast<TH2D*>(theoryFile->Get(theoryName.c_str()));
    theoryPrediction_th2d->SetDirectory(0);
    theoryPrediction_th2d->Draw("col");
  }

  int generateThermalPair()
  {
    // set random variables
    phiVirtPhoton = rand.Uniform(0, o2::constants::math::TwoPI); // flat in phi
    yVirtPhoton = rand.Uniform(minRap, maxRap);                  // flat in rapidity
    massVirtPhoton = hRalf->GetRandom();                         // from calculations of Ralf
    ptVirtPhoton = rand.Uniform(5.);                             // fix... need calculation from Ralf I guess.

    // Calculate eta from rapidity/mass/pt
    etaVirtPhoton = y2eta(ptVirtPhoton, massVirtPhoton, yVirtPhoton);
    // Set vector components for the mother particle
    mother.SetPtEtaPhiM(ptVirtPhoton, etaVirtPhoton, phiVirtPhoton, massVirtPhoton);
    // Set TGenPhaseSpace with mother particle to ee decay
    decayer.SetDecay(mother, 2, mDaughters);
    // Generate decay
    weight = decayer.Generate();
    // Get the two decay results
    d1 = decayer.GetDecay(0);
    d2 = decayer.GetDecay(1);
    return 0;
  }

  Pythia8::Particle createParticle(int pdg, double px, double py, double pz, double E, double m, double xProd, double yProd, double zProd)
  {
    // std::cout << "createParticle() mass " << m << " pdgCode " << pdg << std::endl;
    Pythia8::Particle myparticle;
    myparticle.id(pdg);
    myparticle.status(11); // not sure about this... check with experts.
    myparticle.px(px);
    myparticle.py(py);
    myparticle.pz(pz);
    myparticle.e(E);
    myparticle.m(m);
    myparticle.xProd(xProd);
    myparticle.yProd(yProd);
    myparticle.zProd(zProd);

    return myparticle;
  }

  double y2eta(Double_t pt, Double_t mass, Double_t y)
  {
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(mt / pt * TMath::SinH(y));
  }

private:
  o2::ccdb::CcdbApi ccdb_api;
  TFile theoryFile;
  TH2D* theoryPrediction_th2d;
  TLorentzVector mother;
  TLorentzVector* d1;
  TLorentzVector* d2;
  int sign;
  double mDaughters[2];
  int nGammaPerEvent;
  int nGamma;
  TRandom3 rand;
  TGenPhaseSpace decayer;
  double massVirtPhoton;
  double ptVirtPhoton;
  double yVirtPhoton;
  double etaVirtPhoton;
  double phiVirtPhoton;
  double minRap;
  double maxRap;
  int weight;
};

FairGenerator* generator_thermalRadiation()
{
  auto generator = new GeneratorPythia8ThermalRadiation();
  return reinterpret_cast<FairGenerator*>(generator);
}
