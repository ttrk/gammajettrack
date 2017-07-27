#include "TFile.h"
#include "TH1.h"

#include <fstream>
#include <string>
#include <vector>

#include "systematics.h"

enum SYSVAR
{
    kSYS_PES,
    kSYS_PhoIso,
    kSYS_PhoPurity,
    kSYS_EleRej,
    kSYS_JES,
    kSYS_JER,
    kSYS_Tracking,
    kN_SYS
};

std::string sys_labels[kN_SYS] = {
     "pes", "iso", "purity_up_plus", "ele_rej",
     "jes_up_plus", "jer", "tracking_up_plus"
};

std::string sys_titles[kN_SYS] = {
    "Photon energy scale    ",
    "Photon isolation       ",
    "Photon purity          ",
    "Electron contamination ",
    "Jet energy scale       ",
    "Jet energy resolution  ",
    "Tracking efficiency    "
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

int print_systematics(const char* filelist, const char* label = "", int hiBinMin = 0, int hiBinMax = 20);

int print_systematics(const char* filelist, const char* label, int hiBinMin, int hiBinMax) {
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
    std::vector<double> sys_uncTot = {0, 0, 0, 0};
    // print systematics in Latex format
    std::cout << "\\begin{tabular}{lcccc}" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "Systematic             & \\multicolumn{2}{c}{\\xijet} & \\multicolumn{2}{c}{\\xigamma} \\\\" << std::endl;
    std::cout << "uncertainty            & PbPb        & pp           & PbPb         & pp            \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\hline" << std::endl;

    for (int iSys=0; iSys<SYSVAR::kN_SYS; ++iSys) {

        std::cout << sys_titles[iSys];
        std::string sys_label = sys_labels[iSys];
        std::vector<float> sys_uncs(4);
        for (int iCol = 0; iCol < kN_SYSCOLUMNS; ++iCol) {
            std::string hist_name = Form("h%s_final_%s_%s_%d_%d_%s_ratio_abs", label, sample[iCol].c_str(), reco[iCol].c_str(), hiBinMin, hiBinMax, sys_labels[iSys].c_str());
            h_ratio_abs = (TH1D*)fsys[iCol]->Get(hist_name.c_str());
            sys_uncs[iCol] = 100 * th1_average_content_FF((TH1F*)h_ratio_abs);
            sys_uncTot[iCol] += sys_uncs[iCol]*sys_uncs[iCol];
            if (sys_uncs[iCol] >= 10)        std::cout << Form("& %.1f\\%%    ", sys_uncs[iCol]);
            else if (sys_uncs[iCol] >= 0.1)  std::cout << Form("& %.1f\\%%     ", sys_uncs[iCol]);
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
    if (argc == 5)
        return print_systematics(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
    else if (argc == 3)
        return print_systematics(argv[1], argv[2]);
    else if (argc == 2)
        return print_systematics(argv[1]);
    else
        return 1;
}
