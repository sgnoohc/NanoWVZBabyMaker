#include "AnalysisVariableModule.h"

//==============================================================================================================
//
// AnalysisVariable Module
//
//==============================================================================================================

//##############################################################################################################
void wvzModule::AnalysisVariableModule::AddOutput()
{
    tx->createBranch<int>("nNlep");
    tx->createBranch<int>("lep_ZCand_idx1");
    tx->createBranch<int>("lep_ZCand_idx2");
    tx->createBranch<int>("lep_Nom_idx1");
    tx->createBranch<int>("lep_Nom_idx2");
    tx->createBranch<int>("lep_Nom_idx3");
    tx->createBranch<int>("lep_2ndZCand_idx1");
    tx->createBranch<int>("lep_2ndZCand_idx2");
    tx->createBranch<int>("lep_WCand_idx1");

    // TODO move below to another module
    tx->createBranch<float>("genWeight");
}

//##############################################################################################################
void wvzModule::AnalysisVariableModule::FillOutput()
{
    const int& nVlep = tx->getBranch<int>("nVlep");

    if (nVlep >= 3)
    {

        //===========================================================
        // Finding Z candidates
        //===========================================================

        int lep_ZCand_idx1 = -999;
        int lep_ZCand_idx2 = -999;
        const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
        const vector<float>& lep_pt = tx->getBranch<vector<float>>("lep_pt");
        const vector<int>& lep_id = tx->getBranch<vector<int>>("lep_id");
        const vector<int>& lep_isWVZVeto = tx->getBranch<vector<int>>("lep_isWVZVeto");
        const vector<int>& lep_isWVZNominal = tx->getBranch<vector<int>>("lep_isWVZNominal");

        double compare = 999; // Min value of |Mll - MZ|

        // Loop over the leptons and find the Z boson pair
        for (unsigned int jj = 0; jj < (lep_p4.size() - 1); jj++)
        {

            int id1 = lep_id.at(jj);

            if (not lep_isWVZVeto.at(jj))
                continue;

            for (unsigned int kk = jj + 1; kk < lep_p4.size(); kk++)
            {

                int id2 = lep_id.at(kk);

                if (not lep_isWVZVeto.at(kk))
                    continue;

                LV dilepZCand = lep_p4.at(jj) + lep_p4.at(kk);
                bool ifpass = true;
                if (fabs(dilepZCand.M() - 91.1875) > 10) ifpass = false; // within Z mass window
                if (abs(id1) != abs(id2)) ifpass = false; // same-flavor
                if (id1 == id2)  ifpass = false; // opposite-sign
                if (ifpass && fabs(dilepZCand.M() - 91.1875) > compare) ifpass = false; // if bigger than previous
                if (ifpass)
                {
                    compare = fabs(dilepZCand.M() - 91.1876);
                    lep_ZCand_idx1 = jj;
                    lep_ZCand_idx2 = kk;
                }
            }
        }

        tx->setBranch<int>("lep_ZCand_idx1", lep_ZCand_idx1);
        tx->setBranch<int>("lep_ZCand_idx2", lep_ZCand_idx2);

        //===========================================================
        // Finding "Nominal leptons" (if lep_ZCand_idx1 is found)
        //===========================================================

        if (lep_ZCand_idx1 >= 0 and lep_ZCand_idx2 >= 0)
        {

            // The nominal leptons are the leptons that are not tagged as Z

            // reset
            int lep_Nom_idx1 = -999;
            int lep_Nom_idx2 = -999;
            int lep_Nom_idx3 = -999;
            vector<int> lep_nom_idxs;

            // Loop over the leptons
            for (unsigned int jj = 0 ; jj < lep_pt.size(); jj++)
            {

                if (jj == (unsigned int) lep_ZCand_idx1)
                    continue;

                if (jj == (unsigned int) lep_ZCand_idx2)
                    continue;

                if (not lep_isWVZNominal.at(jj))
                    continue;

                lep_nom_idxs.push_back(jj);

            }

            if (lep_nom_idxs.size() > 0) lep_Nom_idx1 = lep_nom_idxs.at(0);
            if (lep_nom_idxs.size() > 1) lep_Nom_idx2 = lep_nom_idxs.at(1);
            if (lep_nom_idxs.size() > 2) lep_Nom_idx3 = lep_nom_idxs.at(2);

            tx->setBranch<int>("nNlep", lep_nom_idxs.size());
            tx->setBranch<int>("lep_Nom_idx1", lep_Nom_idx1);
            tx->setBranch<int>("lep_Nom_idx2", lep_Nom_idx2);
            tx->setBranch<int>("lep_Nom_idx3", lep_Nom_idx3);

            //===========================================================
            // Finding "2nd Z candidate"
            //===========================================================

            // These are the variables to set
            int lep_2ndZCand_idx1 = -999;
            int lep_2ndZCand_idx2 = -999;
            int lep_WCand_idx1 = -999;

            std::vector<int> good_veto_idx; // index of good nominal leptons

            // Loop over the leptons
            for (unsigned int jj = 0 ; jj < lep_pt.size(); jj++)
            {

                if (jj == (unsigned int) lep_ZCand_idx1)
                    continue;

                if (jj == (unsigned int) lep_ZCand_idx2)
                    continue;

                if (not lep_isWVZVeto.at(jj));

                good_veto_idx.push_back(jj);

            }

            // At this point we must only find 5 total
            if (good_veto_idx.size() == 3)
            {

                int lep_veto_idx1 = good_veto_idx.at(0);
                int lep_veto_idx2 = good_veto_idx.at(1);
                int lep_veto_idx3 = good_veto_idx.at(2);

                // Assumes lep_Nom_idx1,2,3 are set
                double dM12 = fabs((lep_p4.at(lep_veto_idx1) + lep_p4.at(lep_veto_idx2)).M() - 91.1876);
                double dM13 = fabs((lep_p4.at(lep_veto_idx1) + lep_p4.at(lep_veto_idx3)).M() - 91.1876);
                double dM23 = fabs((lep_p4.at(lep_veto_idx2) + lep_p4.at(lep_veto_idx3)).M() - 91.1876);
                double mindM = std::min(dM12, std::min(dM13, dM23));

                if (mindM == dM12)
                {
                    lep_2ndZCand_idx1 = lep_pt.at(lep_veto_idx1) > lep_pt.at(lep_veto_idx2) ? lep_veto_idx1 : lep_veto_idx2;
                    lep_2ndZCand_idx2 = lep_pt.at(lep_veto_idx1) > lep_pt.at(lep_veto_idx2) ? lep_veto_idx2 : lep_veto_idx1;
                    lep_WCand_idx1 = lep_veto_idx3;
                }
                else if (mindM == dM13)
                {
                    lep_2ndZCand_idx1 = lep_pt.at(lep_veto_idx1) > lep_pt.at(lep_veto_idx3) ? lep_veto_idx1 : lep_veto_idx3;
                    lep_2ndZCand_idx2 = lep_pt.at(lep_veto_idx1) > lep_pt.at(lep_veto_idx3) ? lep_veto_idx3 : lep_veto_idx1;
                    lep_WCand_idx1 = lep_veto_idx2;
                }
                else if (mindM == dM23)
                {
                    lep_2ndZCand_idx1 = lep_pt.at(lep_veto_idx2) > lep_pt.at(lep_veto_idx3) ? lep_veto_idx2 : lep_veto_idx3;
                    lep_2ndZCand_idx2 = lep_pt.at(lep_veto_idx2) > lep_pt.at(lep_veto_idx3) ? lep_veto_idx3 : lep_veto_idx2;
                    lep_WCand_idx1 = lep_veto_idx1;
                }
                else
                {
                    RooUtil::error("I should not be here!", "select2ndZCandAndWCandLeptons");
                }

                tx->setBranch<int>("lep_2ndZCand_idx1", lep_2ndZCand_idx1);
                tx->setBranch<int>("lep_2ndZCand_idx2", lep_2ndZCand_idx2);
                tx->setBranch<int>("lep_WCand_idx1", lep_WCand_idx1);

            }

        }

    }

    tx->setBranch<float>("genWeight", nt.genWeight());

}

