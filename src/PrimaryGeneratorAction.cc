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
/// \file electromagnetic/TestEm3/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//#include "random.hh"
#include "CLHEP/Random/RandGauss.h"

#include <math.h>       /* sin */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:G4VUserPrimaryGeneratorAction(),
 fParticleGun(0),
 fDetector(det),   
 fRndmBeam(0.),  
 fGunMessenger(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematic()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(3.*GeV);
  G4double position = -0.5*(fDetector->GetWorldSizeX());
  fParticleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  // Random Gaussian distribution
  // Make a Random distribuition of electron beam
  // Set the rms divergence and rms source size of the e-beam
  // Divergence in mrad
  // Source size in um
  
  std::random_device rd;
  std::mt19937 generator(rd());
  
  // PP experiment parameters
  G4double source_size = 23.25 * um;
  G4double divergence = 1.85 * mrad;
  
  //Laser grating parameters
  //G4double source_size = 2.65 * um;
  //G4double divergence = 2.45 * mrad;
  
  G4double peakenergy = 72 * MeV;
  G4double energyspread = 50 * MeV;
  
  std::normal_distribution<double> posy(0, source_size);
  std::normal_distribution<double> divy(0, divergence);
    
  std::normal_distribution<double> posz(0, source_size);
  std::normal_distribution<double> divz(0, divergence);
  
  std::normal_distribution<double> energydist(peakenergy, energyspread);
  
  G4double energy_elec = energydist(generator);
  
  
  G4double theta_z = divz(generator);
  G4double theta_y = divy(generator);
  
  // Randomize the position
  G4double uz = std::sin(theta_z),
           uy = std::sin(theta_y),
           ux = 1;
  
  G4double z0 = posz(generator),
           y0 = posy(generator),
           x0 = 0;
           
  // Uncomment for debugging
  //G4cout << "z0 = " << z0/um << " um y0 = " << y0/um << " um" << G4endl;
  //G4cout << "theta_z = " << theta_z/mrad << " mrad theta_y = " << theta_y/mrad << " mrad" << G4endl;
  //G4cout << " " << G4endl;
  
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleEnergy(energy_elec); // uncomment if the energy is fixed and not set in the mac file
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

