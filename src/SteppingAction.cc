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
/// \file electromagnetic/TestEm3/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4Positron.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PolarizationHelper.hh"
#include "PrimaryGeneratorAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:G4UserSteppingAction(),fDetector(det),fEventAct(evt) 
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //track informations
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  
  
  ////// Pre-Step
  G4Track* aTrack = aStep->GetTrack(); 
  G4ParticleDefinition* particleType = aTrack->GetDefinition();
  
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  //G4ThreeVector preStepposition = preStepPoint->GetPosition();
  
  G4StepPoint* endStepPoint = aStep->GetPostStepPoint();
  
  //Example of how to get the beam at the boundary
   // Beam at boundary, example
   
   
   
  G4ThreeVector preStepposition  = preStepPoint->GetPosition();
  G4ThreeVector preStepdirection = preStepPoint->GetMomentumDirection();
  G4double preStepkinEnergy = preStepPoint->GetKineticEnergy();
  G4ThreeVector preStepMomentum = preStepPoint->GetMomentum();
    
    
  if (preStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "Screen" &&
      endStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "World") {
      
    //G4cout <<  "Pre: " << preStepPoint->GetTouchableHandle()->GetVolume()->GetName() << " Post: " << endStepPoint->GetTouchableHandle()->GetVolume()->GetName() << G4endl;    
    //G4cout<<"a "<< particleType->GetParticleName() <<" left the Box with " << preStepkinEnergy/MeV << G4endl;
    
    G4AnalysisManager::Instance()->FillH2(0,preStepposition.y()/mm, preStepposition.z()/mm);
    
    // Fill with the particles crossing the screen
    G4double idTuple = 0;
    
    if (particleType == G4Gamma::GammaDefinition())
      idTuple = 3;
    else if (particleType == G4Electron::ElectronDefinition())
      idTuple = 1;
    else if (particleType == G4Positron::PositronDefinition())
      idTuple = 2;    
    
    if(idTuple!=0){
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 0, preStepposition.y()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 1, preStepposition.z()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 2, preStepposition.x()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 3, preStepMomentum.x()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 4, preStepMomentum.y()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 5, preStepMomentum.z()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 6, preStepkinEnergy/MeV);
     G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 7, 0);
     G4AnalysisManager::Instance()->AddNtupleRow(idTuple);
    }
    
    // Record all particles
    idTuple = 4;
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 0, preStepposition.y()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 1, preStepposition.z()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 2, preStepposition.x()/mm);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 3, preStepMomentum.x()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 4, preStepMomentum.y()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 5, preStepMomentum.z()/MeV);
    G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 6, preStepkinEnergy/MeV);
     G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 7, 0);
    G4AnalysisManager::Instance()->AddNtupleRow(idTuple);
      
  } else if (preStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "ScreenSource" &&
      endStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "World") {
    if (particleType == G4Electron::ElectronDefinition()){	
    	// Record Electrons only
	G4double idTuple = 5;
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 0, preStepposition.y()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 1, preStepposition.z()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 2, preStepposition.x()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 3, preStepMomentum.x()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 4, preStepMomentum.y()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 5, preStepMomentum.z()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 6, preStepkinEnergy/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 7, 0);
	G4AnalysisManager::Instance()->AddNtupleRow(idTuple);
    }
  } else if (preStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "ScreenBeforePP" &&
      endStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "World") {
    if (particleType == G4Electron::ElectronDefinition()){	
    	// Record Electrons only
	G4double idTuple = 6;
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 0, preStepposition.y()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 1, preStepposition.z()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 2, preStepposition.x()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 3, preStepMomentum.x()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 4, preStepMomentum.y()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 5, preStepMomentum.z()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 6, preStepkinEnergy/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 7, 0);
	G4AnalysisManager::Instance()->AddNtupleRow(idTuple);
    }
  } else if (preStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "ScreenAfterPP" &&
      endStepPoint->GetTouchableHandle()->GetVolume()->GetName() == "World") {
    if (particleType == G4Electron::ElectronDefinition()){	
    	// Record Electrons only
	G4double idTuple = 7;
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 0, preStepposition.y()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 1, preStepposition.z()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 2, preStepposition.x()/mm);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 3, preStepMomentum.x()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 4, preStepMomentum.y()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 5, preStepMomentum.z()/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 6, preStepkinEnergy/MeV);
	G4AnalysisManager::Instance()->FillNtupleDColumn(idTuple, 7, 0);
	G4AnalysisManager::Instance()->AddNtupleRow(idTuple);
    }
  }
    
    
  //if World, return
  //
  G4VPhysicalVolume* volume = prePoint->GetTouchableHandle()->GetVolume();    
  //if sum of absorbers do not fill exactly a layer: check material, not volume.
  const G4Material* mat = volume->GetLogicalVolume()->GetMaterial();
  if (mat == fDetector->GetWorldMaterial()) return; 

  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 
 
  //here we are in an absorber. Locate it
  //
  G4int absorNum  = prePoint->GetTouchableHandle()->GetCopyNumber(0);
  G4int layerNum  = prePoint->GetTouchableHandle()->GetCopyNumber(1);
  
  //get Run
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
                       
  // collect energy deposit taking into account track weight
  G4double edep = aStep->GetTotalEnergyDeposit()*aStep->GetTrack()->GetWeight();
  
  // collect step length of charged particles
  G4double stepl = 0.;
  if (particle->GetPDGCharge() != 0.) {
    stepl = aStep->GetStepLength();
    run->AddChargedStep();
  } else { run->AddNeutralStep(); }
  
  //  G4cout << "Nabs= " << absorNum << "   edep(keV)= " << edep << G4endl;
  
  // sum up per event
  fEventAct->SumEnergy(absorNum,edep,stepl);
  
  //longitudinal profile of edep per absorber
  if (edep>0.) {
    G4AnalysisManager::Instance()->FillH1(kMaxAbsor+absorNum, 
                                          G4double(layerNum+1), edep);
  }


/*
  //energy flow
  //
  // unique identificator of layer+absorber
  G4int Idnow = (fDetector->GetNbOfAbsor())*layerNum + absorNum;
  G4int plane;
  //
  //leaving the absorber ?
  if (endPoint->GetStepStatus() == fGeomBoundary) {

    if (preStepPoint->GetTouchable()->GetVolume()->GetName() == "Galatic") {
    //G4cout << "Step starts on geometry boundary PLANE0" << G4endl;
    if(particleType == G4Gamma::GammaDefinition()){
      //G4cout << "Gamma!!" << G4endl;
      G4AnalysisManager::Instance()->FillH1(23,position.y()/mm);
      G4AnalysisManager::Instance()->FillH2(0,position.y()/mm, position.z()/mm);
    } 
 }

    //G4cout << preStepPoint->GetTouchable()->GetVolume()->GetName() << G4endl;
    G4ThreeVector position  = endPoint->GetPosition();
    G4ThreeVector direction = endPoint->GetMomentumDirection();
    G4double sizeYZ = 0.5*fDetector->GetCalorSizeYZ();       
    G4double Eflow = endPoint->GetKineticEnergy();
    if(particle == G4Positron::Positron()) Eflow += 2*electron_mass_c2;
    if((std::abs(position.y()) >= sizeYZ) || (std::abs(position.z()) >= sizeYZ))
                                  run->SumLateralEleak(Idnow, Eflow);
    else if (direction.x() >= 0.) run->SumEnergyFlow(plane=Idnow+1, Eflow);
    else                          run->SumEnergyFlow(plane=Idnow,  -Eflow);    
  }
  
  */
  
  
  
  
  
  
     

////  example of Birk attenuation
///G4double destep   = aStep->GetTotalEnergyDeposit();
///G4double response = BirksAttenuation(aStep);
///G4cout << " Destep: " << destep/keV << " keV"
///       << " response after Birks: " << response/keV << " keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
 //Example of Birk attenuation law in organic scintillators.
 //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
 //
 G4Material* material = aStep->GetTrack()->GetMaterial();
 G4double birk1       = material->GetIonisation()->GetBirksConstant();
 G4double destep      = aStep->GetTotalEnergyDeposit();
 G4double stepl       = aStep->GetStepLength();  
 G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 //
 G4double response = destep;
 if (birk1*destep*stepl*charge != 0.)
   {
     response = destep/(1. + birk1*destep/stepl);
   }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

