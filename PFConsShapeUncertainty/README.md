#PFConstituents Uncertainties based on sub detector resolution(Similar to MET uncertainty tool)
1. charged pions (abs(pdgId)=211) - Tracker resolution uncertainty
2. photons (pdgId==22) - ECAL resolution uncertainty 
3. neutral PF hadrons (pdgId==130) - HCAL resolution uncertainty

#Applied (Uncorrelated) on the Signal Extraction Variable: PT Fraction carried by the three leading Particle flow constituents 
#of the Leading (Pencil) Jet
1. First three PF constituents can be π+ π- π0(reco Photon) or Kaon (pdg==130) 
2. P_{T}^{123} Fraction = (P_{T}^{1}+ P_{T}^{2} + P_{T}^{3})/Sum_{i}(P_{T}^{i})
  where 0 < i < No.of PF Constituents of the Leading Jet


