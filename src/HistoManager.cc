//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("testem3")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  
  // Tuple
  analysisManager->SetNtupleMerging(true);
  //
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms
  
  

  // Define histograms start values
  
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","22"};
  G4String title;

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;
  
  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    if (k < kMaxAbsor) title = "Edep in absorber " + id[k];
    if (k > kMaxAbsor) title = "Edep longit. profile (MeV/event) in absorber "
                               + id[k-kMaxAbsor];
    if (k == 2*kMaxAbsor+1) title = "energy flow (MeV/event)";
    if (k == 2*kMaxAbsor+2) title = "lateral energy leak (MeV/event)";
    G4int ih = analysisManager->CreateH1(id[k], title, nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, true);
    G4cout << "Histogram H1: " << ih << G4endl;
  }

  G4int 
ih = analysisManager->
    CreateH1("GammaPlane", "Gamma XY at Plane0", 1000, -100*mm, 100*mm, "mm");
  analysisManager->SetH1Activation(ih, true); //0
  G4cout << "Histogram H2: " << ih << G4endl;

  G4int ih2;
  // Plane 0
  ih2 = analysisManager->
    CreateH2("GammaPlane", "Gamma XY at Plane0", 1000, -5*mm, 5*mm, 1000, -5*mm, 5*mm, "mm", "mm");
  analysisManager->SetH2Activation(0, true); //0
  //G4cout << "Histogram H2: " << 0 << G4endl;
  
  
  //***************************//
  //******** Create nTuples *******//
  analysisManager->SetFirstNtupleId(1);
  //id 1
  analysisManager->CreateNtuple("Electrons_Screen", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  //
  //id 2
  analysisManager->CreateNtuple("Positrons_Screen", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  //
  //id 3
  analysisManager->CreateNtuple("Gamma_Screen", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  //
  //id 4
  analysisManager->CreateNtuple("All_Screen", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  //id 5
  analysisManager->CreateNtuple("Electrons_Screen_source", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  //id 6
  analysisManager->CreateNtuple("Electrons_Screen_beforePP", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
    //id 7
  analysisManager->CreateNtuple("Electrons_Screen_afterPP", "Hits");
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("px"); // column Id = 3
  analysisManager->CreateNtupleDColumn("py"); // column Id = 4
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 5
  analysisManager->CreateNtupleDColumn("energy"); // column Id = 6
  analysisManager->CreateNtupleDColumn("idpart"); // column Id = 7
  analysisManager->FinishNtuple();
  
  
  
}
