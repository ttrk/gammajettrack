#include "TFile.h"
#include "TH1.h"

#include <fstream>
#include <string>
#include <vector>

#include "systematics.h"

enum SYSUNC
{
    k_PES,
    k_PhoIso,
    k_PhoPurity,
    k_EleRej,
    k_phoeff,
    k_JES,
    k_JER,
    k_Tracking,
    k_LongRange,
    k_BkgSub,
    k_nonclosure,
    kN_SYSUNC
};

std::string sys_labels[kN_SYSUNC] = {
     "pes", "iso", "purity_up_plus", "ele_rej", "phoeffcorr",
     "jes_qg_down", "jer", "tracking_up_plus", "longrange", "bkgsub", "xi_nonclosure"
};

std::string sys_labels_ratio[kN_SYSUNC] = {
     "pes", "iso", "purity_up_plus", "ele_rej", "phoeffcorr",
     "jes_qg_down", "jer", "tracking_ratio", "longrange", "bkgsub", "xi_nonclosure"
};

int sysMethod[kN_SYSUNC] = {
    2, 2, 2, 2, 2,
    2, 2, 0, 0, 0, 0
};

std::string sys_titles[kN_SYSUNC] = {
    "Photon energy scale     ",
    "Photon isolation        ",
    "Photon purity           ",
    "Electron contamination  ",
    "Photon efficiency       ",
    "Jet energy scale        ",
    "Jet energy resolution   ",
    "Tracking efficiency     ",
    "Long range contribution ",
    "Background subtraction  ",
    "Non-closure in $\\xi$    ",
};

std::string sys_title_tot = "Total                  ";

enum SYSCOLUMNS {
    k_xijet_pbpb,
    k_xijet_pp,
    k_xijet_ratio,
    k_xigamma_pbpb,
    k_xigamma_pp,
    k_xigamma_ratio,
    kN_SYSCOLUMNS
};

std::string sample[kN_SYSCOLUMNS] = {
    "_pbpbdata", "_ppdata", "_ratio", "_pbpbdata", "_ppdata", "_ratio"
};

std::string jsreco[kN_SYSCOLUMNS] = {
    "_corrjsrecoreco", "_corrjsrecoreco", "", "_corrjsrecoreco", "_corrjsrecoreco", ""
};

std::string ffreco[kN_SYSCOLUMNS] = {
    "_recoreco", "_srecoreco", "", "_recoreco", "_srecoreco", ""
};

int print_systematics(const char* filelist, const char* label = "", int hiBinMin = 0, int hiBinMax = 20, float binMin = 0, float binMax = -1, bool printRatio = false);

int print_systematics(const char* filelist, const char* label, int hiBinMin, int hiBinMax, float binMin, float binMax, bool printRatio) {
    std::string line;

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist);
    if (!file_stream) {
        std::cout << "filelist is not found : " << filelist << std::endl;
        return 1;
    }
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    int nfiles = file_list.size();
    if (nfiles != kN_SYSCOLUMNS) {
        std::cout << "total files is not " << kN_SYSCOLUMNS << std::endl;
        std::cout << "exiting" << std::endl;
        return 1;
    }

    TFile* fsys[nfiles];
    for (int i=0; i<nfiles; ++i) {
        fsys[i] = new TFile(file_list[i].c_str(), "read");
    }

    bool is_js = (strcmp(label, "js") == 0);

    std::string *reco = new std::string[kN_SYSCOLUMNS];
    if (is_js) {
        reco = jsreco;
        sys_labels[k_nonclosure] = "js_nonclosure";
        sys_labels_ratio[k_nonclosure] = "js_nonclosure";
        sys_titles[k_nonclosure] = "Non-closure in $r$      ";
    }
    else {
        reco = ffreco;
    }

    TH1D* h_nom = 0;
    TH1D* h_ratio_abs = 0;
    TH1D* hTmp = 0;

    std::vector<double> sys_uncTot(kN_SYSCOLUMNS, 0);
    // print systematics in Latex format
    if (is_js) {
        std::cout << "\\begin{tabular}{lcccc}" << std::endl;
        std::cout << "\\hline" << std::endl;
        if (!printRatio) {
            if (binMin <= binMax) {
                std::cout << Form("Systematic             & \\multicolumn{2}{c}{\\rhor, %.1f $<$ r $<$ %.1f} \\\\", binMin, binMax) << std::endl;
            }
            else {
                std::cout << "Systematic             & \\multicolumn{2}{c}{\\rhor} \\\\" << std::endl;
            }
            std::cout << "uncertainty            & PbPb        & pp             \\\\" << std::endl;
        }
        else {
            if (binMin <= binMax) {
                std::cout << Form("Systematic             & \\multicolumn{3}{c}{\\rhor, %.1f $<$ r $<$ %.1f} \\\\", binMin, binMax) << std::endl;
            }
            else {
                std::cout << "Systematic             & \\multicolumn{3}{c}{\\rhor} \\\\" << std::endl;
            }
            std::cout << "uncertainty            & PbPb   &  pp   & PbPb/pp    \\\\" << std::endl;
        }
    }
    else {
        std::cout << "\\begin{tabular}{lcccc}" << std::endl;
        std::cout << "\\hline" << std::endl;
        if (!printRatio) {
            std::cout << "Systematic             & \\multicolumn{2}{c}{\\xijet} & \\multicolumn{2}{c}{\\xigamma} \\\\" << std::endl;
            std::cout << "uncertainty            & PbPb        & pp           & PbPb         & pp            \\\\" << std::endl;
        }
        else {
            std::cout << "Systematic             & \\multicolumn{3}{c}{\\xijet} & \\multicolumn{3}{c}{\\xigamma} \\\\" << std::endl;
            std::cout << "uncertainty            & PbPb   &  pp   & PbPb/pp   & PbPb &  pp  && PbPb/pp       \\\\" << std::endl;
        }
    }
    std::cout << "\\hline" << std::endl;
    std::cout << "\\hline" << std::endl;

    for (int iSys=0; iSys<SYSUNC::kN_SYSUNC; ++iSys) {

        //if (iSys == k_JER)  continue;
        if (is_js && iSys == k_LongRange)  continue;
        //if (iSys == k_nonclosure)  continue;

        std::cout << sys_titles[iSys];
        std::vector<float> sys_uncs(kN_SYSCOLUMNS);
        for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {

            if (is_js && iCol >= 3)  continue;

            if (!printRatio) {
                if (iCol == k_xijet_ratio) continue;
                if (iCol == k_xigamma_ratio) continue;
            }

            int hiBinMinTmp = hiBinMin;
            int hiBinMaxTmp = hiBinMax;
            bool isPPCol = (iCol == k_xijet_pp || iCol == k_xigamma_pp);
            if (isPPCol) {
                hiBinMinTmp = 100;
                hiBinMaxTmp = 200;
            }

            std::string hist_name_nom = Form("h%s_final%s%s_%d_%d_nominal", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp);
            h_nom = 0;
            h_nom = (TH1D*)fsys[iCol]->Get(hist_name_nom.c_str());

            std::string sysMethodStr = "";
            if (sysMethod[iSys] == 0) sysMethodStr = "";
            else if (sysMethod[iSys] == 1) sysMethodStr = "_fitBand";
            else if (sysMethod[iSys] == 2) sysMethodStr = "_fit";

            std::string sys_label = sys_labels[iSys];
            if (iCol == k_xijet_ratio || iCol == k_xigamma_ratio)  sys_label = sys_labels_ratio[iSys];

            std::string hist_name_sys = Form("h%s_final%s%s_%d_%d_%s_ratio_abs%s", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp,
                    sys_label.c_str(), sysMethodStr.c_str());
            if (is_js) {
                hist_name_sys = Form("h%s_final%s%s_%d_%d_%s_ratio_fit_diff", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp,
                        sys_label.c_str());
                if (sysMethod[iSys] == 0)
                    hist_name_sys = Form("h%s_final%s%s_%d_%d_%s_ratio_abs", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp,
                        sys_label.c_str());
            }

            h_ratio_abs = 0;
            h_ratio_abs = (TH1D*)fsys[iCol]->Get(hist_name_sys.c_str());

            if (sysMethod[iSys] == 2) h_ratio_abs->Divide(h_nom);

            if (iSys == k_JES) {
                th1_sqrt_sum_squares(h_ratio_abs, h_ratio_abs); // 0.02^2 + 0.02^2

                std::string hist_name_Tmp = Form("h%s_final%s%s_%d_%d_jes_qg_down_ratio_abs%s", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp,
                        sysMethodStr.c_str());
                if (is_js) {
                    hist_name_Tmp = Form("h%s_final%s%s_%d_%d_jes_qg_down_ratio_fit_diff", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp);
                    if (sysMethod[iSys] == 0)
                        hist_name_Tmp = Form("h%s_final%s%s_%d_%d_jes_qg_down_ratio_abs", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp);
                }
                hTmp = (TH1D*)fsys[iCol]->Get(hist_name_Tmp.c_str());

                if (sysMethod[iSys] == 2) hTmp->Divide(h_nom);
                th1_sqrt_sum_squares(h_ratio_abs, hTmp);     // 0.02^2 + 0.02^2 + q/g scale
            }
            int binFirst = 1;
            int binLast = 0;
            if (binMin <= binMax) {
                binFirst = h_ratio_abs->FindBin(binMin);
                binLast = h_ratio_abs->FindBin(binMax) - 1;
            }
            sys_uncs[iCol] = 100 * th1_average_content_FF((TH1D*)h_ratio_abs, binFirst, binLast);
            sys_uncTot[iCol] += sys_uncs[iCol]*sys_uncs[iCol];
            if (sys_uncs[iCol] >= 10)        std::cout << Form("& %.1f\\%%    ", sys_uncs[iCol]);
            else if (sys_uncs[iCol] >= 0.1)  std::cout << Form("& %.1f\\%%     ", sys_uncs[iCol]);
            else if ((iSys == k_PES || iSys == k_BkgSub || iSys == k_nonclosure) && (iCol == k_xijet_pp || iCol == k_xigamma_pp))
                std::cout << Form("& $--$      ");
            else                             std::cout <<      "& $<$0.1\\%    ";
        }
        std::cout << " \\\\" << std::endl;
    }
    std::cout << "\\hline" << std::endl;
    std::cout << sys_title_tot.c_str();
    for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {

        if (is_js && iCol >= 3)  continue;

        if (!printRatio) {
            if (iCol == k_xijet_ratio) continue;
            if (iCol == k_xigamma_ratio) continue;
        }

        double uncTotTmp = sqrt(sys_uncTot[iCol]);
        if (uncTotTmp >= 10)        std::cout << Form("& %.1f\\%%    ", uncTotTmp);
        else if (uncTotTmp >= 0.1)  std::cout << Form("& %.1f\\%%     ", uncTotTmp);
        //else if (uncTotTmp == 0)    std::cout << Form("& --\\%%      ");
        else                        std::cout <<      "& $<$0.1\\%    ";
    }
    std::cout << " \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

    // bin by bin cross-check with total systematics
    if (binMin <= binMax) {
        int binFirst = h_ratio_abs->FindBin(binMin);
        int binLast = h_ratio_abs->FindBin(binMax) - 1;

        if (binFirst == binLast) {
            for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {

                if (is_js && iCol >= 3)  continue;

                if (!printRatio) {
                    if (iCol == k_xijet_ratio) continue;
                    if (iCol == k_xigamma_ratio) continue;
                }

                int hiBinMinTmp = hiBinMin;
                int hiBinMaxTmp = hiBinMax;
                bool isPPCol = (iCol == k_xijet_pp || iCol == k_xigamma_pp);
                if (isPPCol) {
                    hiBinMinTmp = 100;
                    hiBinMaxTmp = 200;
                }

                double uncTotTmp = sqrt(sys_uncTot[iCol]);

                std::string hist_name_totalSys = Form("h%s_final%s%s_%d_%d_systematics", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp);
                hTmp = (TH1D*)fsys[iCol]->Get(hist_name_totalSys.c_str());
                double diffTotTmp = hTmp->GetBinContent(binFirst);

                std::string hist_name_nominal = Form("h%s_final%s%s_%d_%d_nominal", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMinTmp, hiBinMaxTmp);
                hTmp = (TH1D*)fsys[iCol]->Get(hist_name_nominal.c_str());
                double nominalTmp = hTmp->GetBinContent(binFirst);

                double uncTotHist = 100 * diffTotTmp / nominalTmp;

                std::cout << "iCol = " << iCol << ", binMin = " << binMin << ", binMax = " << binMax;
                std::cout << ", uncTotTmp = " << uncTotTmp << ", uncTotHist = " << uncTotHist;
                std::cout << ", uncTotTmp - uncTotHist = " << uncTotTmp - uncTotHist << ". ";
                if (std::fabs(uncTotTmp - uncTotHist) < 0.000001) {
                    std::cout << "Total systematics agree.";
                }
                else {
                    std::cout << "Warning : Total systematics do NOT match.";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int i=0; i<nfiles; ++i) {
        fsys[i]->Close();
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 8)
        return print_systematics(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]), atoi(argv[7]));
    else if (argc == 7)
        return print_systematics(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]));
    else if (argc == 5)
        return print_systematics(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
    else if (argc == 3)
        return print_systematics(argv[1], argv[2]);
    else if (argc == 2)
        return print_systematics(argv[1]);
    else
        return 1;
}
