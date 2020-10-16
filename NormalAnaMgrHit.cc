
#include "NormalAnaMgrHit.hh"
//  for event
#include <sstream>
#include <cassert>
#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "NormalTrackInfo.hh"
#include "EvtNavigator/NavBuffer.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/SimHeader.h"

#include "DataCollSvc/IDataCollSvc.h"
static IDataCollSvc* m_datacollsvc;
static std::string key;
#include "JunoTimer/IJunoTimerSvc.h"
#include "JunoTimer/JunoTimer.h"
static IJunoTimerSvc* m_timersvc;
static JunoTimerPtr m_timer_step;
static JunoTimerPtr m_timer_beginrun;
static JunoTimerPtr m_timer_endrun;
static JunoTimerPtr m_timer_beginevent;
static JunoTimerPtr m_timer_endevent;
static JunoTimerPtr m_timer_begintrack;
static JunoTimerPtr m_timer_endtrack;



DECLARE_TOOL(NormalAnaMgrHit);

NormalAnaMgrHit::NormalAnaMgrHit(const std::string& name) 
    : ToolBase(name)
{
    declProp("EnableNtuple", m_flag_ntuple=true);
    m_evt_tree = 0;
    m_step_no = 0;
}

NormalAnaMgrHit::~NormalAnaMgrHit()
{

}

void
NormalAnaMgrHit::BeginOfRunAction(const G4Run* /*aRun*/) {
    if (not m_flag_ntuple) {
        return;
    }


//
   SniperPtr<IDataCollSvc> _datacollsvc(this->getParent(), "DataCollSvc");
    if (_datacollsvc.invalid()) {
        LogError << "Can't Locate DataCollSvc. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
    }
    m_datacollsvc = _datacollsvc.data();

   SniperPtr<IJunoTimerSvc> _timersvc(this->getParent(), "JunoTimerSvc");
    if (_timersvc.invalid()) {
        LogError << "Can't Locate JunoTimerSvc. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
    }
    m_timersvc = _timersvc.data();
    m_timer_step = m_timersvc->get("steppingtimer");
    m_timer_beginrun=m_timersvc->get("beginruntimer");
    m_timer_endrun=m_timersvc->get("endruntimer");
    m_timer_beginevent=m_timersvc->get("begineventtimer");
    m_timer_endevent=m_timersvc->get("endeventtimer");
    m_timer_begintrack=m_timersvc->get("begintracktimer");
    m_timer_endtrack=m_timersvc->get("endtracktimer"); 
//





    // check the RootWriter is Valid.
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }
    m_timer_beginrun->start();

    m_evt_tree = svc->bookTree("SIMEVT/evt", "evt");
    m_evt_tree->Branch("evtID", &m_eventID, "evtID/I");
    m_evt_tree->Branch("nPhotons", &m_nPhotons, "nPhotons/I");
   // m_evt_tree->Branch("totalPE", &m_totalPE, "totalPE/I");
   
    m_evt_tree->Branch("nPE", m_nPE, "nPE[nPhotons]/I");
    m_evt_tree->Branch("energy", m_energy, "energy[nPhotons]/F");
    m_evt_tree->Branch("hitTime", m_hitTime, "hitTime[nPhotons]/D");
    m_evt_tree->Branch("pmtID", m_pmtID, "pmtID[nPhotons]/I");
    m_evt_tree->Branch("PETrackID", m_peTrackID, "PETrackID[nPhotons]/I");
   
   /* m_evt_tree->Branch("edep", &m_edep, "edep/F");
    m_evt_tree->Branch("edepX", &m_edep_x, "edepX/F");
    m_evt_tree->Branch("edepY", &m_edep_y, "edepY/F");
    m_evt_tree->Branch("edepZ", &m_edep_z, "edepZ/F");
  */
    m_evt_tree->Branch("isCerenkov", m_isCerenkov, "isCerenkov[nPhotons]/I");
    m_evt_tree->Branch("isReemission", m_isReemission, "isReemission[nPhotons]/I");
    m_evt_tree->Branch("isOriginalOP", m_isOriginalOP, "isOriginalOP[nPhotons]/I");
    m_evt_tree->Branch("OriginalOPTime", m_OriginalOPTime, "OriginalOPTime[nPhotons]/D");
  
    // PMT
    
    m_evt_tree->Branch("nPMTs", &m_npmts_byPMT, "nPMTs/I");
    m_evt_tree->Branch("nPE_byPMT", m_nPE_byPMT, "nPE_byPMT[nPMTs]/I");
    m_evt_tree->Branch("PMTID_byPMT", m_PMTID_byPMT, "PMTID_byPMT[nPMTs]/I");
   
    // - 2015.10.10 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's position
    m_evt_tree->Branch("LocalPosX", m_localpos_x, "LocalPosX[nPhotons]/F");
    m_evt_tree->Branch("LocalPosY", m_localpos_y, "LocalPosY[nPhotons]/F");
    m_evt_tree->Branch("LocalPosZ", m_localpos_z, "LocalPosZ[nPhotons]/F");
   
    // - 2016.04.17 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's direction
    m_evt_tree->Branch("LocalDirX", m_localdir_x, "LocalDirX[nPhotons]/F");
    m_evt_tree->Branch("LocalDirY", m_localdir_y, "LocalDirY[nPhotons]/F");
    m_evt_tree->Branch("LocalDirZ", m_localdir_z, "LocalDirZ[nPhotons]/F");
    
    // - 2017.03.01 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's Global Position
    m_evt_tree->Branch("GlobalPosX", m_globalpos_x, "GlobalPosX[nPhotons]/F");
    m_evt_tree->Branch("GlobalPosY", m_globalpos_y, "GlobalPosY[nPhotons]/F");
    m_evt_tree->Branch("GlobalPosZ", m_globalpos_z, "GlobalPosZ[nPhotons]/F");

    m_evt_tree->Branch("BoundaryPosX", m_boundarypos_x, "BoundaryPosX[nPhotons]/F");
    m_evt_tree->Branch("BoundaryPosY", m_boundarypos_y, "BoundaryPosY[nPhotons]/F");
    m_evt_tree->Branch("BoundaryPosZ", m_boundarypos_z, "BoundaryPosZ[nPhotons]/F");
  
    m_step_no = new TH1I("stepno", "step number of optical photons", 1000, 0, 1000);
    svc->attach("SIMEVT", m_step_no);
  
    m_timer_beginrun->stop();  
    key = "t_beginrun";
    m_datacollsvc->collectData(key, m_timer_beginrun->elapsed()); 
}

void
NormalAnaMgrHit::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
NormalAnaMgrHit::BeginOfEventAction(const G4Event* evt) {
    // initialize the evt tree
    m_timer_beginevent->start();
  
    m_eventID = evt->GetEventID();
    m_nPhotons = 0;
 //   m_totalPE = 0;
    for(int i = 0; i < 2000000; i++) {
      m_nPE[i] = 0;
      m_energy[i] = 0;
      m_hitTime[i] = 0;
      m_pmtID[i] = 0;
      m_peTrackID[i] = 0;
      m_isCerenkov[i] = 0;
      m_isReemission[i] = 0;
      m_isOriginalOP[i] = 0;
      m_OriginalOPTime[i] = 0;

      m_localpos_x[i] = 0.;
      m_localpos_y[i] = 0.;
      m_localpos_z[i] = 0.;

      m_localdir_x[i] = 0.;
      m_localdir_y[i] = 0.;
      m_localdir_z[i] = 0.;

      m_boundarypos_x[i] = 0.;
      m_boundarypos_y[i] = 0.;
      m_boundarypos_z[i] = 0.;
    }
    
   /* m_edep = 0.;
    m_edep_x = 0.;
    m_edep_y = 0.;
    m_edep_z = 0.;
   */
    m_cache_bypmt.clear();
  
    m_timer_beginevent->stop();
    key = "t_beginevent";
    m_datacollsvc->collectData(key, m_timer_beginevent->elapsed());
}

void
NormalAnaMgrHit::EndOfEventAction(const G4Event* evt) {
    m_timer_endevent->start();

    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    junoHit_PMT_Collection* col = 0; 
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    if (!HCE or CollID<0) {
        LogError << "No hits collection found." << std::endl;
    } else {
        col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));
    }
    
    // fill evt data
  //  int totPE = 0;
    if (col) {
        int n_hit = col->entries();
        m_nPhotons = n_hit;
        // FIXME: Make sure not overflow
        if (n_hit > 2000000) { m_nPhotons = 2000000; }

        for (int i = 0; i < n_hit; ++i) {
           // totPE += (*col)[i]->GetCount(); 
            // if overflow, don't save anything into the array.
            // but still count the totalPE.
            if (i >= 2000000) { continue; }
            m_energy[i] = (*col)[i]->GetKineticEnergy();
            m_nPE[i] = (*col)[i]->GetCount();
            m_hitTime[i] = (*col)[i]->GetTime();
            m_pmtID[i] = (*col)[i]->GetPMTID();

            m_cache_bypmt[m_pmtID[i]] += m_nPE[i];
           

            if ((*col)[i]->IsFromCerenkov()) {
                LogDebug << "+++++ from cerenkov" << std::endl;
                m_isCerenkov[i] = 1;
            }
            if ((*col)[i]->IsReemission()) {
                LogDebug << "+++++ reemission" << std::endl;
                m_isReemission[i] = 1;
            }

            m_isOriginalOP[i] = (*col)[i]->IsOriginalOP();
            m_OriginalOPTime[i] = (*col)[i]->GetOriginalOPStartT();
            m_peTrackID[i] = (*col)[i]->GetProducerID();

            G4ThreeVector local_pos = (*col)[i]->GetPosition();
            m_localpos_x[i] = local_pos.x();
            m_localpos_y[i] = local_pos.y();
            m_localpos_z[i] = local_pos.z();

            G4ThreeVector local_dir = (*col)[i]->GetMomentum();
            m_localdir_x[i] = local_dir.x();
            m_localdir_y[i] = local_dir.y();
            m_localdir_z[i] = local_dir.z();

            G4ThreeVector global_pos = (*col)[i]->GetGlobalPosition();
            m_globalpos_x[i] = global_pos.x();
            m_globalpos_y[i] = global_pos.y();
            m_globalpos_z[i] = global_pos.z();

            G4ThreeVector boundary_pos = (*col)[i]->GetBoundaryPosition();
            m_boundarypos_x[i] = boundary_pos.x();
            m_boundarypos_y[i] = boundary_pos.y();
            m_boundarypos_z[i] = boundary_pos.z();
          
        }

    }

     m_npmts_byPMT = 0;
    for (std::map<int,int>::iterator it = m_cache_bypmt.begin();
            it != m_cache_bypmt.end(); ++it) {
        m_PMTID_byPMT[m_npmts_byPMT] = it->first;
        m_nPE_byPMT[m_npmts_byPMT] = it->second;
        ++m_npmts_byPMT;
    }
  

  //  m_totalPE = totPE;

   /* if (m_edep>0) {
        m_edep_x /= m_edep;
        m_edep_y /= m_edep;
        m_edep_z /= m_edep;
    }
   */
    if (m_flag_ntuple and m_evt_tree) {
        m_evt_tree -> Fill();
    }
    
    save_into_data_model();
    m_timer_endevent->stop();
    key = "t_endevent";
    m_datacollsvc->collectData(key, m_timer_endevent->elapsed());
}


void
NormalAnaMgrHit::PreUserTrackingAction(const G4Track* aTrack) {
 
}

void
NormalAnaMgrHit::PostUserTrackingAction(const G4Track* aTrack) {
}

void
NormalAnaMgrHit::UserSteppingAction(const G4Step* step) {
    m_timer_step->start();  
 
   G4Track* track = step->GetTrack();
    // if the step number of optical photon bigger than X, mark it as killed
    if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
        G4int stepno = track->GetCurrentStepNumber();

        if (track->GetTrackStatus() == fStopAndKill) {
            // if the opticalphoton is killed, save the step no
            m_step_no->Fill(stepno);
        }

        if (stepno >= 1000) 
          {
            G4String phyname;
            if (track->GetVolume()) { phyname = track->GetVolume()->GetName(); }
            const G4ThreeVector& tmppos = track->GetPosition();
            LogWarn << "opticalphoton [" << track->GetTrackID() << "]"
                    << "@[" << phyname << "]"
                    << " (" << tmppos.x() << ", "
                    << tmppos.y() << ", "
                    << tmppos.z() << ", "
                    << track->GetGlobalTime() << ") "
                    << " step number >= " << 1000
                    << std::endl;
            track->SetTrackStatus(fStopAndKill);
          }
        
        
    }
   
    m_timer_step->stop();
    key = "t_step";
    m_datacollsvc->collectData(key, m_timer_step->elapsed());

}

bool NormalAnaMgrHit::save_into_data_model() {
    SniperDataPtr<JM::NavBuffer>  navBuf(*getParent(), "/Event");
    if (navBuf.invalid()) {
        return false;
    }
    LogDebug << "navBuf: " << navBuf.data() << std::endl;
    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    LogDebug << "evt_nav: " << evt_nav << std::endl;
    if (not evt_nav) {
        return false;
    }
    JM::SimHeader* m_simheader = dynamic_cast<JM::SimHeader*>(evt_nav->getHeader("/Event/Sim"));
    LogDebug << "simhdr: " << m_simheader << std::endl;
    if (not m_simheader) {
        return false;
    }
    JM::SimEvent* m_simevent = dynamic_cast<JM::SimEvent*>(m_simheader->event());
    LogDebug << "simevt: " << m_simevent << std::endl;
    if (not m_simevent) {
        return false;
    }

    // DO NOTHING
    return true;
}
