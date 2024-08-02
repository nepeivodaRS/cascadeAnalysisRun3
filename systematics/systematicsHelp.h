#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>

void saveJsonToFile(const std::string& json, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << json;
        file.close();
        std::cout << "JSON saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

std::string generateJson(float bachBaryonCosPA, float bachBaryonDCAxyToPV, float casccospa, float cascradius,
                         float dcabachtopv, float dcacascdau, float dcanegtopv, float dcapostopv, float dcav0dau,
                         float dcav0topv, float lambdamasswin, float masswin, float masswintpc, float minpt,
                         int mintpccrrows, float nsigmatpcKa, float nsigmatpcPi, float nsigmatpcPr,
                         float nsigmatofKa, float nsigmatofPi, float nsigmatofPr, float proplifetime, float ptthrtof,
                         float rejcomp, float v0cospa, float v0radius, float hastof, int isMC,
                         float isSelectBachBaryon, int isXi, int evSelFlag, int processGen, int processRec, std::string fileList)
{
    std::ostringstream json;

    json << "{" << std::endl;
    json << "    \"internal-dpl-clock\": \"\"," << std::endl;
    json << "    \"internal-dpl-aod-reader\": {" << std::endl;
    json << "        \"time-limit\": \"0\"," << std::endl;
    json << "        \"orbit-offset-enumeration\": \"0\"," << std::endl;
    json << "        \"orbit-multiplier-enumeration\": \"0\"," << std::endl;
    json << "        \"start-value-enumeration\": \"0\"," << std::endl;
    json << "        \"end-value-enumeration\": \"-1\"," << std::endl;
    json << "        \"step-value-enumeration\": \"1\"," << std::endl;
    json << "        \"aod-file\": \"@" << fileList << "\"" << std::endl;
    json << "    }," << std::endl;
    json << "    \"lf-cascpostprocessing\": {" << std::endl;
    json << "        \"bachBaryonCosPA\": \"" << std::fixed << std::setprecision(1) << bachBaryonCosPA << "\"," << std::endl;
    json << "        \"bachBaryonDCAxyToPV\": \"" << std::fixed << std::setprecision(3) << bachBaryonDCAxyToPV << "\"," << std::endl;
    json << "        \"casccospa\": \"" << std::fixed << std::setprecision(4) << casccospa << "\"," << std::endl;
    json << "        \"cascradiusMin\": \"" << std::fixed << std::setprecision(2) << cascradius << "\"," << std::endl;
    json << "        \"cascradiusMax\": \"" << std::fixed << std::setprecision(2) << 1000 << "\"," << std::endl;
    json << "        \"dcabachtopv\": \"" << std::fixed << std::setprecision(2) << dcabachtopv << "\"," << std::endl;
    json << "        \"dcacascdau\": \"" << std::fixed << std::setprecision(2) << dcacascdau << "\"," << std::endl;
    json << "        \"dcanegtopv\": \"" << std::fixed << std::setprecision(2) << dcanegtopv << "\"," << std::endl;
    json << "        \"dcapostopv\": \"" << std::fixed << std::setprecision(2) << dcapostopv << "\"," << std::endl;
    json << "        \"dcav0dau\": \"" << std::fixed << std::setprecision(1) << dcav0dau << "\"," << std::endl;
    json << "        \"dcav0topv\": \"" << std::fixed << std::setprecision(2) << dcav0topv << "\"," << std::endl;
    json << "        \"evSelFlag\": \"" << std::fixed << std::setprecision(0) << evSelFlag << "\"," << std::endl;
    json << "        \"hastof\": \"" << std::fixed << std::setprecision(0) << hastof << "\"," << std::endl;
    json << "        \"isMC\": \"" << std::fixed << std::setprecision(0) << isMC << "\"," << std::endl;
    json << "        \"isSelectBachBaryon\": \"" << std::fixed << std::setprecision(0) << isSelectBachBaryon << "\"," << std::endl;
    json << "        \"isXi\": \"" << std::fixed << std::setprecision(0) << isXi << "\"," << std::endl;
    json << "        \"lambdamasswin\": \"" << std::fixed << std::setprecision(4) << lambdamasswin << "\"," << std::endl;
    json << "        \"masswin\": \"" << std::fixed << std::setprecision(3) << masswin << "\"," << std::endl;
    json << "        \"masswintpc\": \"" << std::fixed << std::setprecision(3) << masswintpc << "\"," << std::endl;
    json << "        \"minpt\": \"" << std::fixed << std::setprecision(2) << minpt << "\"," << std::endl;
    json << "        \"mintpccrrows\": \"" << std::fixed << std::setprecision(1) << mintpccrrows << "\"," << std::endl;
    json << "        \"nsigmatofKa\": \"" << std::fixed << std::setprecision(1) << nsigmatofKa << "\"," << std::endl;
    json << "        \"nsigmatofPi\": \"" << std::fixed << std::setprecision(1) << nsigmatofPi << "\"," << std::endl;
    json << "        \"nsigmatofPr\": \"" << std::fixed << std::setprecision(1) << nsigmatofPr << "\"," << std::endl;
    json << "        \"nsigmatpcKa\": \"" << std::fixed << std::setprecision(1) << nsigmatpcKa << "\"," << std::endl;
    json << "        \"nsigmatpcPi\": \"" << std::fixed << std::setprecision(1) << nsigmatpcPi << "\"," << std::endl;
    json << "        \"nsigmatpcPr\": \"" << std::fixed << std::setprecision(1) << nsigmatpcPr << "\"," << std::endl;
    json << "        \"processGen\": \"" << std::fixed << std::setprecision(0) << processGen << "\"," << std::endl;
    json << "        \"processRec\": \"" << std::fixed << std::setprecision(0) << processRec << "\"," << std::endl;
    json << "        \"proplifetime\": \"" << std::fixed << std::setprecision(1) << proplifetime << "\"," << std::endl;
    json << "        \"ptthrtof\": \"" << std::fixed << std::setprecision(1) << ptthrtof << "\"," << std::endl;
    json << "        \"rejcomp\": \"" << std::fixed << std::setprecision(3) << rejcomp << "\"," << std::endl;
    json << "        \"v0cospa\": \"" << std::fixed << std::setprecision(4) << v0cospa << "\"," << std::endl;
    json << "        \"v0radiusMin\": \"" << std::fixed << std::setprecision(2) << v0radius << "\"," << std::endl;
    json << "        \"v0radiusMax\": \"" << std::fixed << std::setprecision(2) << 1000 << "\"," << std::endl;
    json << "        \"dooobrej\": \"" << std::fixed << std::setprecision(0) << 0 << "\"" << std::endl;
    json << "    }," << std::endl;
    json << "    \"internal-dpl-injected-dummy-sink\": \"\"," << std::endl;
    json << "    \"internal-dpl-aod-global-analysis-file-sink\": \"\"" << std::endl;
    json << "}" << std::endl;

    return json.str();
}

double generateRandomX(double min, double max)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dis(min,max);

    return dis(gen);
}

int nsigmamc = 100;

// Default values (Xi)
float bachBaryonCosPA = 2;
float bachBaryonDCAxyToPV = 0.0199999996;
float casccospa = 0.9947;
float cascradius = 1.01;
float dcabachtopv = 0.04;
float dcacascdau = 0.25;
float dcanegtopv = 0.06;
float dcapostopv = 0.04;
float dcav0dau = 0.5;
float dcav0topv = 0.03;
float etadau = 0.8;
float evSelFlag = 2;
float hastof = 0;
float isSelectBachBaryon = 1;
float isXi = 1;
float lambdamasswin = 0.0116;
float masswin = 0.050;
float masswintpc = 0.012;
float minpt = 0;
int mintpccrrows = 52;
float nsigmatofKa = 100;
float nsigmatofPi = 100;
float nsigmatofPr = 100;
float nsigmatpcKa = 4.7;
float nsigmatpcPi = 4.8;
float nsigmatpcPr = 5.0;
float proplifetime = 4.6;
float ptthrtof = 0;
float rap = 0.5;
float rejcomp = 0;
float v0cospa = 0.9876;
float v0radius = 0.55;

//Limits for 1by1 study:
double casccospalim[2] = {0.969, 0.999}; // >
double cascradiuslim[2] = {0.4, 1.3}; // >
double dcabachtopvlim[2] = {0.0, 0.1}; // >
double dcacascdaulim[2] = {0.1, 1.0}; // < 
double dcanegtopvlim[2] = {0.00, 0.10}; // >
double dcapostopvlim[2] = {0.00, 0.10}; // >
double dcav0daulim[2] = {0.1, 1.0}; // <
double dcav0topvlim[2] = {0.0, 0.1}; // >
double lambdamasswinlim[2] = {0.0010, 0.0200}; // <
double nsigmatpcKlim[2] = {2, 6}; // <
double nsigmatpcPilim[2] = {2, 6}; // <
double nsigmatpcPrlim[2] = {2, 6}; // <
double proplifetimelim[2] = {2, 8}; // <
double v0cospalim[2] = {0.969, 0.999}; // >
double v0radiuslim[2] = {0.4, 2.5}; // >
double tpccrrowslim[2] = {50, 80}; // >
double rejcomplim[2] = {0.000, 0.010}; // >

//Limits for multitrial:
double casccospalimMultiTrial[2] = {0.97, 0.995};
double cascradiuslimMultiTrial[2] = {0.4, 1.3};
double dcabachtopvlimMultiTrial[2] = {0.02, 0.1};
double dcacascdaulimMultiTrial[2] = {0.2, 1.00};
double dcanegtopvlimMultiTrial[2] = {0.01, 0.10};
double dcapostopvlimMultiTrial[2] = {0.01, 0.10};
double dcav0daulimMultiTrial[2] = {0.1, 1.0};
double dcav0topvlimMultiTrial[2] = {0., 0.1};
double lambdamasswinlimMultiTrial[2] = {0.006, 0.014};
double nsigmatpcKlimMultiTrial[2] = {4, 6};
double nsigmatpcPilimMultiTrial[2] = {4, 6};
double nsigmatpcPrlimMultiTrial[2] = {4, 6};
double proplifetimelimMultiTrial[2] = {4, 8};
double v0cospalimMultiTrial[2] = {0.97, 0.99};
double v0radiuslimMultiTrial[2] = {0.4, 2.5};
double tpccrrowslimMultiTrial[2] = {50, 80};
double rejcomplimMultiTrial[2] = {0., 0.};

//Limits for PurityMaxim:
double casccospalimPurityMaxim[2] = {0.980, 0.995};
double cascradiuslimPurityMaxim[2] = {0.4, 1.1};
double dcabachtopvlimPurityMaxim[2] = {0.02, 0.06};
double dcacascdaulimPurityMaxim[2] = {0.2, 0.4};
double dcanegtopvlimPurityMaxim[2] = {0.02, 0.06};
double dcapostopvlimPurityMaxim[2] = {0.02, 0.06};
double dcav0daulimPurityMaxim[2] = {0.15, 0.60};
double dcav0topvlimPurityMaxim[2] = {0.00, 0.04};
double lambdamasswinlimPurityMaxim[2] = {0.007, 0.012};
double nsigmatpcKlimPurityMaxim[2] = {4, 5};
double nsigmatpcPilimPurityMaxim[2] = {4, 5};
double nsigmatpcPrlimPurityMaxim[2] = {4, 5};
double proplifetimelimPurityMaxim[2] = {4, 8};
double v0cospalimPurityMaxim[2] = {0.97, 0.99};
double v0radiuslimPurityMaxim[2] = {0.5, 1.5};
double tpccrrowslimPurityMaxim[2] = {50, 65};
double rejcomplimPurityMaxim[2] = {0., 0.};

std::string pathToDATA = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_data.txt";
std::string pathToMC = "/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_mc.txt";
