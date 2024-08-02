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

    const int njson = 11;

    std::string filename[njson];
    std::string filenameMC[njson];

    for (int i = 0; i < njson; i++)
    {
        // Assign to the variables random values within their limits
        // casccospa = generateRandomX(casccospalim[0], casccospalim[1]);
        // cascradius = generateRandomX(cascradiuslim[0], cascradiuslim[1]);
        // dcabachtopv = generateRandomX(dcabachtopvlim[0], dcabachtopvlim[1]);
        // dcacascdau = generateRandomX(dcacascdaulim[0], dcacascdaulim[1]);
        // dcanegtopv = generateRandomX(dcanegtopvlim[0], dcanegtopvlim[1]);
        // dcapostopv = generateRandomX(dcapostopvlim[0], dcapostopvlim[1]);
        // dcav0dau = generateRandomX(dcav0daulim[0], dcav0daulim[1]);
        // dcav0topv = generateRandomX(dcav0topvlim[0], dcav0topvlim[1]);
        // lambdamasswin = generateRandomX(lambdamasswinlim[0], lambdamasswinlim[1]);
        // nsigmatpcKa = generateRandomX(nsigmatpcKlim[0], nsigmatpcKlim[1]);
        // nsigmatpcPi = generateRandomX(nsigmatpcPilim[0], nsigmatpcPilim[1]);
        // nsigmatpcPr = generateRandomX(nsigmatpcPrlim[0], nsigmatpcPrlim[1]);
        // proplifetime = generateRandomX(proplifetimelim[0], proplifetimelim[1]);
        // rejcomp = generateRandomX(rejcomplim[0], rejcomplim[1]);
        // v0cospa = generateRandomX(v0cospalim[0], v0cospalim[1]);
        // v0radius = generateRandomX(v0radiuslim[0], v0radiuslim[1]);
        // mintpccrrows = generateRandomX(tpccrrowslim[0], tpccrrowslim[1]);

        rejcomp = rejcomplim[0] + i * (rejcomplim[1] - rejcomplim[0]) / (njson - 1) ;


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

        filename[i] = "data/config/systjson" + std::to_string(i+1) + ".json";
        filenameMC[i] = "mc/config/systjson" + std::to_string(i+1) + ".json";
        saveJsonToFile(json, filename[i]);
        saveJsonToFile(jsonMC, filenameMC[i]);
    }
}
