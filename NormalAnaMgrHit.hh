#ifndef NormalAnaMgrHit_hh
#define NormalAnaMgrHit_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include "TH1I.h"
#include <map>

class NormalAnaMgrHit: public IAnalysisElement,
                    public ToolBase{
public:

    NormalAnaMgrHit(const std::string& name);
    ~NormalAnaMgrHit();
    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);

    virtual void UserSteppingAction(const G4Step* step);

private:
    bool save_into_data_model();

private:

    // Evt Data
    TTree* m_evt_tree;
   Int_t m_eventID;
    Int_t m_nPhotons;
 //   Int_t m_totalPE;
    Int_t m_nPE[2000000];
    Float_t m_energy[2000000];
    Double_t m_hitTime[2000000];
    Int_t m_pmtID[2000000];
    Int_t m_peTrackID[2000000];
    Int_t m_isCerenkov[2000000];
    Int_t m_isReemission[2000000];
    Int_t m_isOriginalOP[2000000];
    Double_t m_OriginalOPTime[2000000];
  
    // PMT Info
    //
    Int_t m_npmts_byPMT;
    Int_t m_nPE_byPMT[2000000];
    Int_t m_PMTID_byPMT[2000000];
    std::map<int, int> m_cache_bypmt;
    // - 2015.10.10 Tao Lin
    //   The hit's local position in PMT will be saved.
    //   However, to save the space, Float is enough.
   
    Float_t m_localpos_x[2000000];
    Float_t m_localpos_y[2000000];
    Float_t m_localpos_z[2000000];
    
    // - 2016.04.17 Tao Lin
    //   hit's local direction in PMT
   
    Float_t m_localdir_x[2000000];
    Float_t m_localdir_y[2000000];
    Float_t m_localdir_z[2000000];
   
    // - 2017.03.01 Tao Lin
    //   hit's global position in PMT
    
    Float_t m_globalpos_x[2000000];
    Float_t m_globalpos_y[2000000];
    Float_t m_globalpos_z[2000000];

    Float_t m_boundarypos_x[2000000];
    Float_t m_boundarypos_y[2000000];
    Float_t m_boundarypos_z[2000000];
    
/*
    Float_t m_edep;
    Float_t m_edep_x;
    Float_t m_edep_y;
    Float_t m_edep_z;
*/
    bool m_flag_ntuple;

    TH1I* m_step_no;
};

#endif
