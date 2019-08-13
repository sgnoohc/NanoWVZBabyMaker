#include "METModule.h"

//==============================================================================================================
//
// MET Module
//
//==============================================================================================================

void wvzModule::METModule::AddOutput()
{

    tx->createBranch<float>("met_orig_pt");
    tx->createBranch<float>("met_orig_phi");
    tx->createBranch<LV>("met_p4");
    tx->createBranch<float>("met_pt");
    tx->createBranch<float>("met_phi");
    tx->createBranch<float>("met_gen_pt");
    tx->createBranch<float>("met_gen_phi");

    tx->createBranch<int>("passesMETfiltersRun2");
}

void wvzModule::METModule::FillOutput()
{

    correctMET();

    tx->setBranch<float>("met_orig_pt", nt.MET_pt());
    tx->setBranch<float>("met_orig_phi", nt.MET_phi());
    if (ana.is_data)
    {
        tx->setBranch<LV>("met_p4", RooUtil::Calc::getLV(nt.MET_pt(), 0, nt.MET_phi(), 0));
        tx->setBranch<float>("met_pt", nt.MET_pt());
        tx->setBranch<float>("met_phi", nt.MET_phi());
    }
    else
    {
        tx->setBranch<LV>("met_p4", RooUtil::Calc::getLV(metobj_corrected.extras.met, 0, metobj_corrected.extras.phi, 0));
        tx->setBranch<float>("met_pt", metobj_corrected.extras.met);
        tx->setBranch<float>("met_phi", metobj_corrected.extras.phi);
    }
    tx->setBranch<float>("met_gen_pt", nt.GenMET_pt());
    tx->setBranch<float>("met_gen_phi", nt.GenMET_phi());

    tx->setBranch<int>("passesMETfiltersRun2", passesMETfilters(ana.is_data));
}

//______________________________________________________________________________________________
void wvzModule::METModule::correctMET()
{
    metobj.extras.met = metobj.extras.met_original = metobj.extras.met_raw
    = metobj.extras.met_METup = metobj.extras.met_METdn
    = metobj.extras.met_JERup = metobj.extras.met_JERdn
    = metobj.extras.met_PUup = metobj.extras.met_PUdn
    = metobj_corrected.extras.met = metobj_corrected.extras.met_original = metobj_corrected.extras.met_raw
    = metobj_corrected.extras.met_METup = metobj_corrected.extras.met_METdn
    = metobj_corrected.extras.met_JERup = metobj_corrected.extras.met_JERdn
    = metobj_corrected.extras.met_PUup = metobj_corrected.extras.met_PUdn
    = nt.MET_pt(); // This is ptmiss from the MET recipe alone (nominal).

    metobj.extras.phi = metobj.extras.phi_original = metobj.extras.phi_raw
    = metobj.extras.phi_METup = metobj.extras.phi_METdn
    = metobj.extras.phi_JECup = metobj.extras.phi_JECdn
    = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
    = metobj.extras.phi_PUup = metobj.extras.phi_PUdn
    = metobj_corrected.extras.phi = metobj_corrected.extras.phi_original = metobj_corrected.extras.phi_raw
    = metobj_corrected.extras.phi_METup = metobj_corrected.extras.phi_METdn
    = metobj_corrected.extras.phi_JECup = metobj_corrected.extras.phi_JECdn
    = metobj_corrected.extras.phi_JERup = metobj_corrected.extras.phi_JERdn
    = metobj_corrected.extras.phi_PUup = metobj_corrected.extras.phi_PUdn
    = nt.MET_phi(); // Nominal MET phi from MET recipe alone

    // metobj.extras.met_JECup = metobj_corrected.extras.met_JECup = babymaker->coreMET.met_up_pt; // MET JES up
    // metobj.extras.met_JECdn = metobj_corrected.extras.met_JECdn = babymaker->coreMET.met_dn_pt; // MET JES dn
    // metobj.extras.phi_JECup = metobj_corrected.extras.phi_JECup = babymaker->coreMET.met_up_phi; // MET phi JES up
    // metobj.extras.phi_JECdn = metobj_corrected.extras.phi_JECdn = babymaker->coreMET.met_dn_phi; // MET phi JES dn

    metobj.extras.met_JECup = metobj_corrected.extras.met_JECup = nt.MET_pt(); // MET JES up
    metobj.extras.met_JECdn = metobj_corrected.extras.met_JECdn = nt.MET_pt(); // MET JES dn
    metobj.extras.phi_JECup = metobj_corrected.extras.phi_JECup = nt.MET_phi(); // MET phi JES up
    metobj.extras.phi_JECdn = metobj_corrected.extras.phi_JECdn = nt.MET_phi(); // MET phi JES dn

    metcorrector.correctMET(nt.GenMET_pt(), nt.GenMET_phi(), &metobj_corrected, false);

}
