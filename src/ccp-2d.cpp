#include "ccp-2d.hpp"

int ccp2d::runPlasma2D(int argc, char** argv)
{
    try {
        std::string caseFile = "cases/lymberopoulos1993.case";
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "-c" && i + 1 < argc) caseFile = argv[++i];
            else if (a == "-h" || a == "--help") {
                std::cout << "Usage: " << argv[0] << " [-c case-file]\n";
                return 0;
            }
        }

        Config cfg = readConfig(caseFile);
        std::cout << std::fixed << std::setprecision(6)
                  << "[Config]\n"
                  << "  outputDir    = " << cfg.outputDir << "\n"
                  << "  ionCSV       = " << cfg.ionCSV << "\n"
                  << "  excCSV       = " << cfg.excCSV << "\n"
                  << "  meshFile     = " << cfg.meshFile << "\n"
                  << "  bccFile      = " << cfg.bccFile << "\n"
                  << "  Hion/Hexc[eV]= " << cfg.Hion << " / " << cfg.Hexc << "\n"
                  << "  stopPeriod   = " << cfg.stopPeriod << "\n"
                  << "  dtPhys       = " << std::scientific << cfg.dtPhys << "\n"
                  << "  pseudoCFL    = " << std::fixed << cfg.pseudoCFL << "\n"
                  << "  nInnerDT     = " << cfg.nInnerDT << "\n"
                  << "  nTIDT        = " << cfg.nTIDT << "\n";

        nondimensionalize(cfg);
        printNondimensionalizationInfo(cfg, std::cout);

        Solver solver(cfg);
        solver.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "fatal: " << e.what() << std::endl;
        return 1;
    }
}

int main(int argc, char** argv)
{
    return ccp2d::runPlasma2D(argc, argv);
}
