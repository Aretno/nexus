// ----------------------------------------------------------------------------
// nexus | Electroluminescence.cc
//
// This class describes the generation of the EL light.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Electroluminescence.h"

#include "IonizationElectron.h"
#include "BaseDriftField.h"

#include <G4MaterialPropertiesTable.hh>
#include <G4ParticleChange.hh>
#include <G4OpticalPhoton.hh>
#include <Randomize.hh>
#include <G4Poisson.hh>
#include <G4GenericMessenger.hh>
#include <G4Scintillation.hh>
#include "globals.hh"

#include <CLHEP/Units/PhysicalConstants.h>

using namespace nexus;
using namespace CLHEP;



Electroluminescence::Electroluminescence(const G4String& process_name,
					                               G4ProcessType type):
  G4VDiscreteProcess(process_name, type), theFastIntegralTable_(0),
  table_generation_(false), photons_per_point_(0)
{
  ParticleChange_ = new G4ParticleChange();
  pParticleChange = ParticleChange_;

  BuildThePhysicsTable();

   /// Messenger
  msg_ = new G4GenericMessenger(this, "/Physics/Electroluminescence/",
				"Control commands of the Electroluminescence physics process.");
  msg_->DeclareProperty("table_generation", table_generation_,
			"EL Table generation");
  msg_->DeclareProperty("photons_per_point", photons_per_point_,
			"Photon per point");

 }



Electroluminescence::~Electroluminescence()
{
  delete theFastIntegralTable_;
}



G4bool Electroluminescence::IsApplicable(const G4ParticleDefinition& pdef)
{
  return (pdef == *IonizationElectron::Definition());
}



G4VParticleChange*
Electroluminescence::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // Initialize particle change with current track values
  ParticleChange_->Initialize(track);

  // Get the current region and its associated drift field.
  // If no drift field is defined, kill the track and leave
  G4Region* region = track.GetVolume()->GetLogicalVolume()->GetRegion();
  BaseDriftField* field =
    dynamic_cast<BaseDriftField*>(region->GetUserInformation());
  if (!field) {
    ParticleChange_->ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  // Get the light yield from the field
  const G4double yield = field->LightYield();
  G4double step_length = field->GetTotalDriftLength();

  if (yield <= 0.)
    return G4VDiscreteProcess::PostStepDoIt(track, step);

  // Generate a random number of photons around mean 'yield'
  G4double mean = yield * step_length;

  G4int num_photons;
  G4int total_num_photons;

  if (yield < 10.) { // Poissonian regime
    total_num_photons = G4int(G4Poisson(mean));
  }
  else {             // Gaussian regime
    G4double sigma = sqrt(mean);
    total_num_photons = G4int(G4RandGauss::shoot(mean, sigma) + 0.5);
  }

  if (table_generation_)
    total_num_photons = photons_per_point_;

  ParticleChange_->SetNumberOfSecondaries(total_num_photons);

  // Track secondaries first to avoid a memory bloat
  if ((total_num_photons > 0) && (track.GetTrackStatus() == fAlive))
    ParticleChange_->ProposeTrackStatus(fSuspend);


  //////////////////////////////////////////////////////////////////
  G4StepPoint* pPostStepPoint = step.GetPostStepPoint();
  G4StepPoint* pPreStepPoint  = step.GetPreStepPoint();
  G4ThreeVector position = pPreStepPoint->GetPosition();
  G4double time = pPreStepPoint->GetGlobalTime();
  G4LorentzVector initial_position(position, time);

  G4ThreeVector position_end = pPostStepPoint->GetPosition();
  G4double time_end = pPostStepPoint->GetGlobalTime();
  G4LorentzVector final_position(position_end, time_end);

  // Energy is sampled from integral (like it is
  // done in G4Scintillation)
  G4Material* mat = pPostStepPoint->GetTouchable()->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4MaterialPropertiesTable* mpt = mat->GetMaterialPropertiesTable();
  const G4MaterialPropertyVector* spectrum = mpt->GetProperty("ELSPECTRUM");
  if (!spectrum) return G4VDiscreteProcess::PostStepDoIt(track, step);

  G4PhysicsOrderedFreeVector* spectrum_integral =
    (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable_)(mat->GetIndex());

  G4double sc_max = spectrum_integral->GetMaxValue();
  
  // Add time components to EL
  G4int N_timeconstants = 1;

   if(mpt->GetConstProperty("ELTIMECONSTANT3"))
     N_timeconstants = 3;
   else if(mpt->GetConstProperty("ELTIMECONSTANT2"))
     N_timeconstants = 2;
   else if(!(mpt->GetConstProperty("ELTIMECONSTANT1")))
   {
     // no components were specified
     return G4VDiscreteProcess::PostStepDoIt(track, step);
   }

  G4double yield1     = 0.;
  G4double yield2     = 0.;
  G4double yield3     = 0.;
  G4double sum_yields = 0.;  
  
  yield1 = mpt->GetConstProperty("ELYIELD1")
                ? mpt->GetConstProperty("ELYIELD1")
                : 1.;
  yield2 = mpt->ConstPropertyExists("ELYIELD2")
                ? mpt->GetConstProperty("ELYIELD2")
                : 0.;
  yield3 = mpt->ConstPropertyExists("ELYIELD3")
                ? mpt->GetConstProperty("ELYIELD3")
                : 0.;

  sum_yields = yield1 + yield2 + yield3;

  G4double scint_time                 = 0.;
  G4double rise_time                  = 0.;

for(G4int scnt = 0; scnt < N_timeconstants; ++scnt)
  {
    // if there is 1 time constant it is #1, etc.
    if(scnt == 0)
    {
      if(N_timeconstants == 1)
      {
        num_photons = total_num_photons;
      }
      else
      {
        num_photons = yield1 / sum_yields * total_num_photons;
      }
      scint_time = mpt->GetConstProperty("ELTIMECONSTANT1");
    }
    else if(scnt == 1)
    {
      // to be consistent with old version (due to double->int conversion)
      if(N_timeconstants == 2)
      {
        num_photons = total_num_photons - num_photons;
      }
      else
      {
        num_photons = yield2 / sum_yields * total_num_photons;
      }
      scint_time = mpt->GetConstProperty("ELTIMECONSTANT2");
    }
    else if(scnt == 2)
    {
      num_photons   = yield3 / sum_yields * total_num_photons;
      scint_time = mpt->GetConstProperty("ELTIMECONSTANT3");
    }

    for (G4int i=0; i<num_photons; i++) {
      // Generate a random direction for the photon
      // (EL is supposed isotropic)
      G4double cos_theta = 1. - 2.*G4UniformRand();
      G4double sin_theta = sqrt((1.-cos_theta)*(1.+cos_theta));

      G4double phi = twopi * G4UniformRand();
      G4double sin_phi = sin(phi);
      G4double cos_phi = cos(phi);

      G4double px = sin_theta * cos_phi;
      G4double py = sin_theta * sin_phi;
      G4double pz = cos_theta;

      G4ThreeVector momentum(px, py, pz);

      // Determine photon polarization accordingly
      G4double sx = cos_theta * cos_phi;
      G4double sy = cos_theta * sin_phi;
      G4double sz = -sin_theta;

      G4ThreeVector polarization(sx, sy, sz);
      G4ThreeVector perp = momentum.cross(polarization);

      phi = twopi * G4UniformRand();
      sin_phi = sin(phi);
      cos_phi = cos(phi);

      polarization = cos_phi * polarization + sin_phi * perp;
      polarization = polarization.unit();

      // Generate a new photon and set properties
      G4DynamicParticle* photon =
        new G4DynamicParticle(G4OpticalPhoton::Definition(), momentum);

      photon->
        SetPolarization(polarization.x(), polarization.y(), polarization.z());

      // Determine photon energy
      G4double sc_value = G4UniformRand()*sc_max;
      G4double sampled_energy = spectrum_integral->GetEnergy(sc_value);
      photon->SetKineticEnergy(sampled_energy);

      G4LorentzVector xyzt =
        field->GeneratePointAlongDriftLine(initial_position, final_position);

      // Time decay 
      G4double deltaTime = 0;
      if(rise_time == 0.0)
      {
        deltaTime -= scint_time * std::log(G4UniformRand());
      }
      else
      {
        deltaTime += SampleTime(rise_time, scint_time);
      }

      xyzt.setT(xyzt.t() + deltaTime);      

      // Create the track
      G4Track* secondary = new G4Track(photon, xyzt.t(), xyzt.v());
      secondary->SetParentID(track.GetTrackID());
      ParticleChange_->AddSecondary(secondary);

  }
  }
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}



void Electroluminescence::BuildThePhysicsTable()
{
  if (theFastIntegralTable_) return;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table

  if(!theFastIntegralTable_)
    theFastIntegralTable_ = new G4PhysicsTable(numOfMaterials);

  for (G4int i=0 ; i<numOfMaterials; i++) {

  	G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
  	  new G4PhysicsOrderedFreeVector();

  	// Retrieve vector of scintillation wavelength intensity for
  	// the material from the material's optical properties table.

  	G4Material* material = (*theMaterialTable)[i];

    G4MaterialPropertiesTable* mpt = material->GetMaterialPropertiesTable();

    if (mpt) {

  	  G4MaterialPropertyVector* theFastLightVector =
  	    mpt->GetProperty("ELSPECTRUM");

  	  if (theFastLightVector) {
        ComputeCumulativeDistribution(*theFastLightVector, *aPhysicsOrderedFreeVector);
		  }
  	}

  	// The scintillation integral(s) for a given material
  	// will be inserted in the table(s) according to the
  	// position of the material in the material table.

  	theFastIntegralTable_->insertAt(i,aPhysicsOrderedFreeVector);
  }
}



void Electroluminescence::ComputeCumulativeDistribution(
  const G4PhysicsOrderedFreeVector& pdf, G4PhysicsOrderedFreeVector& cdf)
{
  G4double sum = 0.;
  cdf.InsertValues(pdf.Energy(0), sum);

  for (unsigned int i=1; i<pdf.GetVectorLength(); ++i) {
    G4double area =
      0.5 * (pdf.Energy(i) - pdf.Energy(i-1)) * (pdf[i] + pdf[i-1]);
    sum = sum + area;
    cdf.InsertValues(pdf.Energy(i), sum);
  }
}

G4double Electroluminescence::SampleTime(G4double tau1, G4double tau2)
 {
   // tau1: rise time and tau2: decay time
   // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
   while(true)
   {
     G4double ran1 = G4UniformRand();
     G4double ran2 = G4UniformRand();
 
     // exponential distribution as envelope function: very efficient
     G4double d = (tau1 + tau2) / tau2;
     // make sure the envelope function is
     // always larger than the bi-exponential
     G4double t  = -1.0 * tau2 * std::log(1. - ran1);
     G4double gg = d * single_exp(t, tau2);
     if(ran2 <= bi_exp(t, tau1, tau2) / gg)
       return t;
   }
   return -1.0;
 }


G4double Electroluminescence::GetMeanFreePath(const G4Track&, G4double,
                                              G4ForceCondition* condition)
{
  *condition = StronglyForced;
  return DBL_MAX;
}
