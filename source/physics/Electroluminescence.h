// ----------------------------------------------------------------------------
// nexus | Electroluminescence.h
//
// This class describes the generation of the EL light.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef ELECTROLUMINESCENCE_H
#define ELECTROLUMINESCENCE_H

#include <G4VDiscreteProcess.hh>
#include <G4PhysicsOrderedFreeVector.hh>

class G4ParticleChange;
class G4GenericMessenger;


namespace nexus {

  class Electroluminescence: public G4VDiscreteProcess
  {
  public:
    /// Constructor
    Electroluminescence(const G4String& process_name = "Electroluminescence",
			                  G4ProcessType type=fUserDefined);
    /// Destructor
    ~Electroluminescence();

    /// Returns true if particle is an ionization electron
    G4bool IsApplicable(const G4ParticleDefinition&);

  public:
    /// This is the method that implements the EL light emission
    /// as a post-step process, that is, photons are generated as
    /// secondaries at the end of the step.
    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

  private:

    /// Returns infinity; i.e., the process does not limit the step,
    /// but sets the 'StronglyForced' condition for the DoIt to be
    /// invoked at every step.
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
    G4double SampleTime(G4double tau1, G4double tau2);
    void BuildThePhysicsTable();
    void ComputeCumulativeDistribution(const G4PhysicsOrderedFreeVector&,
                                       G4PhysicsOrderedFreeVector&);

    G4double single_exp(G4double t, G4double tau2);
    G4double bi_exp(G4double t, G4double tau1, G4double tau2);


  private:
    G4ParticleChange* ParticleChange_;

    G4PhysicsTable* theFastIntegralTable_;

    G4GenericMessenger* msg_;

    G4bool table_generation_;
    G4int photons_per_point_;
  };


inline G4double Electroluminescence::single_exp(G4double t, G4double tau2)
 {
   return std::exp(-1.0 * t / tau2) / tau2;
 }
 
 inline G4double Electroluminescence::bi_exp(G4double t, G4double tau1,
                                         G4double tau2)
 {
   return std::exp(-1.0 * t / tau2) * (1 - std::exp(-1.0 * t / tau1)) / tau2 /
          tau2 * (tau1 + tau2);
 }
 } // end namespace nexus

#endif
