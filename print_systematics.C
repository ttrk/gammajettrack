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
    k_JES,
    k_JER,
    k_Tracking,
    k_LongRange,
    k_BkgSub,
    k_xi_nonclosure,
    kN_SYSUNC
};

std::string sys_labels[kN_SYSUNC] = {
     "pes", "iso", "purity_up_plus", "ele_rej",
     "jes_up_plus", "jer", "tracking_up_plus", "longrange", "bkgsub", "xi_nonclosure"
};

std::string sys_titles[kN_SYSUNC] = {
    "Photon energy scale     ",
    "Photon isolation        ",
    "Photon purity           ",
    "Electron contamination  ",
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
    k_xigamma_pbpb,
    k_xigamma_pp,
    kN_SYSCOLUMNS
};

std::string sample[kN_SYSCOLUMNS] = {
    "pbpbdata", "ppdata", "pbpbdata", "ppdata"
};

std::vector<std::string> jsreco = {
    "recoreco", "srecoreco", "recoreco", "srecoreco"
};

std::vector<std::string> ffreco = {
    "recoreco", "srecoreco", "recoreco", "srecoreco"
};

int print_systematics(const char* filelist, const char* label = "", int hiBinMin = 0, int hiBinMax = 20, float xiBinMin = 0, float xiBinMax = -1);

int print_systematics(const char* filelist, const char* label, int hiBinMin, int hiBinMax, float xiBinMin, float xiBinMax) {
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

    std::vector<std::string> reco;
    if (strcmp(label, "js") == 0)
        reco = jsreco;
    else
        reco = ffreco;

    TH1D* h_ratio_abs = 0;
    TH1D* hTmp = 0;
    std::vector<double> sys_uncTot = {0, 0, 0, 0};
    // print systematics in Latex format
    std::cout << "\\begin{tabular}{lcccc}" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "Systematic             & \\multicolumn{2}{c}{\\xijet} & \\multicolumn{2}{c}{\\xigamma} \\\\" << std::endl;
    std::cout << "uncertainty            & PbPb        & pp           & PbPb         & pp            \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\hline" << std::endl;

    for (int iSys=0; iSys<SYSUNC::kN_SYSUNC; ++iSys) {

        std::cout << sys_titles[iSys];
        std::vector<float> sys_uncs(4);
        for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {
            std::string hist_name = Form("h%s_final_%s_%s_%d_%d_%s_ratio_abs", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMin, hiBinMax, sys_labels[iSys].c_str());
            h_ratio_abs = (TH1D*)fsys[iCol]->Get(hist_name.c_str());
            if (iSys == k_JES) {
                th1_sqrt_sum_squares(h_ratio_abs, h_ratio_abs); // 0.02^2 + 0.02^2

                std::string hist_name_Tmp = Form("h%s_final_%s_%s_%d_%d_jes_qg_down_ratio_abs", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMin, hiBinMax);
                hTmp = (TH1D*)fsys[iCol]->Get(hist_name_Tmp.c_str());
                th1_sqrt_sum_squares(h_ratio_abs, hTmp);     // 0.02^2 + 0.02^2 + q/g scale
            }
            int binFirst = 1;
            int binLast = 0;
            if (xiBinMin <= xiBinMax) {
                binFirst = h_ratio_abs->FindBin(xiBinMin);
                binLast = h_ratio_abs->FindBin(xiBinMax) - 1;
            }
            sys_uncs[iCol] = 100 * th1_average_content_FF((TH1D*)h_ratio_abs, binFirst, binLast);
            sys_uncTot[iCol] += sys_uncs[iCol]*sys_uncs[iCol];
            if (sys_uncs[iCol] >= 10)        std::cout << Form("& %.1f\\%%    ", sys_uncs[iCol]);
            else if (sys_uncs[iCol] >= 0.1)  std::cout << Form("& %.1f\\%%     ", sys_uncs[iCol]);
            else if ((iSys == k_BkgSub || iSys == k_xi_nonclosure) && (iCol == k_xijet_pp || iCol == k_xigamma_pp))
                std::cout << Form("& $--$      ");
            else                             std::cout <<      "& $<$0.1\\%    ";
        }
        std::cout << " \\\\" << std::endl;
    }
    std::cout << "\\hline" << std::endl;
    std::cout << sys_title_tot.c_str();
    for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {
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

    for (int i=0; i<nfiles; ++i) {
        fsys[i]->Close();
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 7)
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
