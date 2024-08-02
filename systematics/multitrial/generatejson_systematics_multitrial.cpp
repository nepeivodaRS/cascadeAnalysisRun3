#include "../systematicsHelp.h"
// The macro generates N different jsons with N different cut settings
// when given as input the limits for each variable

void generatejson_systematics_multitrial()
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

    // Loosest
    std::string jsonloosest = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospalimMultiTrial[0], cascradiuslimMultiTrial[0], dcabachtopvlimMultiTrial[0],
                                       dcacascdaulimMultiTrial[1], dcanegtopvlimMultiTrial[0], dcapostopvlimMultiTrial[0], dcav0daulimMultiTrial[1], dcav0topvlimMultiTrial[0], lambdamasswinlimMultiTrial[1], masswin,
                                       masswintpc, minpt, tpccrrowslimMultiTrial[0], nsigmatpcKlimMultiTrial[1], nsigmatpcPilimMultiTrial[1], nsigmatpcPrlimMultiTrial[1],
                                       nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetimelimMultiTrial[1], ptthrtof, rejcomplimMultiTrial[0], v0cospalimMultiTrial[0],
                                       v0radiuslimMultiTrial[0], hastof, 0, isSelectBachBaryon, isXi, evSelFlag, 0,
                                       1,
                                       pathToDATA);

    std::string jsonloosestMC = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospalimMultiTrial[0], cascradiuslimMultiTrial[0], dcabachtopvlimMultiTrial[0],
                                       dcacascdaulimMultiTrial[1], dcanegtopvlimMultiTrial[0], dcapostopvlimMultiTrial[0], dcav0daulimMultiTrial[1], dcav0topvlimMultiTrial[0], lambdamasswinlimMultiTrial[1], masswin,
                                       masswintpc, minpt, tpccrrowslimMultiTrial[0], nsigmamc, nsigmamc, nsigmamc,
                                       nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetimelimMultiTrial[1], ptthrtof, rejcomplimMultiTrial[0], v0cospalimMultiTrial[0],
                                       v0radiuslimMultiTrial[0], hastof, 1, isSelectBachBaryon, isXi, evSelFlag, 1,
                                       1,
                                       pathToMC);


    std::string filenameDEF;
    std::string filenameDEFMC;
    filenameDEF = "data/config/configDEF.json";
    filenameDEFMC = "mc/config/configDEF.json";

    std::string filenameLoosestDEF;
    std::string filenameLoosestDEFMC;
    filenameLoosestDEF = "data/config/configLOOSEST.json";
    filenameLoosestDEFMC = "mc/config/configLOOSEST.json";

    saveJsonToFile(jsondef, filenameDEF);
    saveJsonToFile(jsondefMC, filenameDEFMC);

    saveJsonToFile(jsonloosest, filenameLoosestDEF);
    saveJsonToFile(jsonloosestMC, filenameLoosestDEFMC);

    const int njson = 1000;

    std::string filename[njson];
    std::string filenameMC[njson];

    for (int i = 0; i < njson; i++)
    {
        // Assign to the variables random values within their limits
        casccospa = generateRandomX(casccospalimMultiTrial[0], casccospalimMultiTrial[1]);
        cascradius = generateRandomX(cascradiuslimMultiTrial[0], cascradiuslimMultiTrial[1]);
        dcabachtopv = generateRandomX(dcabachtopvlimMultiTrial[0], dcabachtopvlimMultiTrial[1]);
        dcacascdau = generateRandomX(dcacascdaulimMultiTrial[0], dcacascdaulimMultiTrial[1]);
        dcanegtopv = generateRandomX(dcanegtopvlimMultiTrial[0], dcanegtopvlimMultiTrial[1]);
        dcapostopv = generateRandomX(dcapostopvlimMultiTrial[0], dcapostopvlimMultiTrial[1]);
        dcav0dau = generateRandomX(dcav0daulimMultiTrial[0], dcav0daulimMultiTrial[1]);
        dcav0topv = generateRandomX(dcav0topvlimMultiTrial[0], dcav0topvlimMultiTrial[1]);
        lambdamasswin = generateRandomX(lambdamasswinlimMultiTrial[0], lambdamasswinlimMultiTrial[1]);
        nsigmatpcKa = generateRandomX(nsigmatpcKlimMultiTrial[0], nsigmatpcKlimMultiTrial[1]);
        nsigmatpcPi = generateRandomX(nsigmatpcPilimMultiTrial[0], nsigmatpcPilimMultiTrial[1]);
        nsigmatpcPr = generateRandomX(nsigmatpcPrlimMultiTrial[0], nsigmatpcPrlimMultiTrial[1]);
        proplifetime = generateRandomX(proplifetimelimMultiTrial[0], proplifetimelimMultiTrial[1]);
        rejcomp = generateRandomX(rejcomplimMultiTrial[0], rejcomplimMultiTrial[1]);
        v0cospa = generateRandomX(v0cospalimMultiTrial[0], v0cospalimMultiTrial[1]);
        v0radius = generateRandomX(v0radiuslimMultiTrial[0], v0radiuslimMultiTrial[1]);
        mintpccrrows = generateRandomX(tpccrrowslimMultiTrial[0], tpccrrowslimMultiTrial[1]);

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
