#include "../../systematicsHelp.h"
// The macro generates N different jsons with N different cut settings
// when given as input the limits for each variable

// Usage: root -l -b -q generatejson_systematics.cpp

void generatejson_systematics()
{

    // DEFAULT JSON
    std::string jsondef = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospa, cascradius, dcabachtopv,
                                       dcacascdau, dcanegtopv, dcapostopv, dcav0dau, dcav0topv, lambdamasswin, masswin,
                                       masswintpc, minpt, mintpccrrows, nsigmatpcKa, nsigmatpcPi, nsigmatpcPr,
                                       nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetime, ptthrtof, rejcomp, v0cospa,
                                       v0radius, hastof, 0, isSelectBachBaryon, isXi, evSelFlag, 0,
                                       1,
                                       pathToDATA);

    std::string jsondefMC = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospa, cascradius, dcabachtopv,
                                         dcacascdau, dcanegtopv, dcapostopv, dcav0dau, dcav0topv, lambdamasswin, masswin,
                                         masswintpc, minpt, mintpccrrows, nsigmamc, nsigmamc, nsigmamc,
                                         nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetime, ptthrtof, rejcomp, v0cospa,
                                         v0radius, hastof, 1, isSelectBachBaryon, isXi, evSelFlag, 1,
                                         1,
                                         pathToMC);

    std::string filenameDEF;
    std::string filenameDEFMC;
    filenameDEF = "data/config/configDEF.json";
    filenameDEFMC = "mc/config/configDEF.json";

    saveJsonToFile(jsondef, filenameDEF);
    saveJsonToFile(jsondefMC, filenameDEFMC);

    const int njson1 = 16;
    const int njson2 = 16;

    std::string filename[njson1*njson2];
    std::string filenameMC[njson1*njson2];

    for (int i = 0; i < njson1; i++)
    {
        casccospa = casccospalim[0] + i * (casccospalim[1] - casccospalim[0]) / (njson1 - 1) ;
        for (int j = 0; j < njson2; j++)
        {
            mintpccrrows = tpccrrowslim[0] + j * (tpccrrowslim[1] - tpccrrowslim[0]) / (njson2 - 1) ;

            std::string json = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospa, cascradius, dcabachtopv,
                                            dcacascdau, dcanegtopv, dcapostopv, dcav0dau, dcav0topv, lambdamasswin, masswin,
                                            masswintpc, minpt, mintpccrrows, nsigmatpcKa, nsigmatpcPi, nsigmatpcPr,
                                            nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetime, ptthrtof, rejcomp, v0cospa,
                                            v0radius, hastof, 0, isSelectBachBaryon, isXi, evSelFlag, 0,
                                            1,
                                            pathToDATA);

            std::string jsonMC = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospa, cascradius, dcabachtopv,
                                              dcacascdau, dcanegtopv, dcapostopv, dcav0dau, dcav0topv, lambdamasswin, masswin,
                                              masswintpc, minpt, mintpccrrows, nsigmamc, nsigmamc, nsigmamc,
                                              nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetime, ptthrtof, rejcomp, v0cospa,
                                              v0radius, hastof, 1, isSelectBachBaryon, isXi, evSelFlag, 1,
                                              1,
                                              pathToMC);

            filename[i*njson2 + j] = "data/config/systjson" + std::to_string(i*njson2 + j + 1) + ".json";
            filenameMC[i*njson2 + j] = "mc/config/systjson" + std::to_string(i*njson2 + j + 1) + ".json";
            saveJsonToFile(json, filename[i*njson2 + j]);
            saveJsonToFile(jsonMC, filenameMC[i*njson2 + j]);
        }
    }
}
