#include "../systematicsHelp.h"
// The macro generates N different jsons with N different cut settings
// when given as input the limits for each variable

void generatejson_systematics_purity()
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
    std::string jsonloosest = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospalimPurityMaxim[0], cascradiuslimPurityMaxim[0], dcabachtopvlimPurityMaxim[0],
                                       dcacascdaulimPurityMaxim[1], dcanegtopvlimPurityMaxim[0], dcapostopvlimPurityMaxim[0], dcav0daulimPurityMaxim[1], dcav0topvlimPurityMaxim[0], lambdamasswinlimPurityMaxim[1], masswin,
                                       masswintpc, minpt, tpccrrowslimPurityMaxim[0], nsigmatpcKlimPurityMaxim[1], nsigmatpcPilimPurityMaxim[1], nsigmatpcPrlimPurityMaxim[1],
                                       nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetimelimPurityMaxim[1], ptthrtof, rejcomplimPurityMaxim[0], v0cospalimPurityMaxim[0],
                                       v0radiuslimPurityMaxim[0], hastof, 0, isSelectBachBaryon, isXi, evSelFlag, 0,
                                       1,
                                       pathToDATA);

    std::string jsonloosestMC = generateJson(bachBaryonCosPA, bachBaryonDCAxyToPV, casccospalimPurityMaxim[0], cascradiuslimPurityMaxim[0], dcabachtopvlimPurityMaxim[0],
                                       dcacascdaulimPurityMaxim[1], dcanegtopvlimPurityMaxim[0], dcapostopvlimPurityMaxim[0], dcav0daulimPurityMaxim[1], dcav0topvlimPurityMaxim[0], lambdamasswinlimPurityMaxim[1], masswin,
                                       masswintpc, minpt, tpccrrowslimPurityMaxim[0], nsigmamc, nsigmamc, nsigmamc,
                                       nsigmatofKa, nsigmatofPi, nsigmatofPr, proplifetimelimPurityMaxim[1], ptthrtof, rejcomplimPurityMaxim[0], v0cospalimPurityMaxim[0],
                                       v0radiuslimPurityMaxim[0], hastof, 1, isSelectBachBaryon, isXi, evSelFlag, 1,
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
        casccospa = generateRandomX(casccospalimPurityMaxim[0], casccospalimPurityMaxim[1]);
        cascradius = generateRandomX(cascradiuslimPurityMaxim[0], cascradiuslimPurityMaxim[1]);
        dcabachtopv = generateRandomX(dcabachtopvlimPurityMaxim[0], dcabachtopvlimPurityMaxim[1]);
        dcacascdau = generateRandomX(dcacascdaulimPurityMaxim[0], dcacascdaulimPurityMaxim[1]);
        dcanegtopv = generateRandomX(dcanegtopvlimPurityMaxim[0], dcanegtopvlimPurityMaxim[1]);
        dcapostopv = generateRandomX(dcapostopvlimPurityMaxim[0], dcapostopvlimPurityMaxim[1]);
        dcav0dau = generateRandomX(dcav0daulimPurityMaxim[0], dcav0daulimPurityMaxim[1]);
        dcav0topv = generateRandomX(dcav0topvlimPurityMaxim[0], dcav0topvlimPurityMaxim[1]);
        lambdamasswin = generateRandomX(lambdamasswinlimPurityMaxim[0], lambdamasswinlimPurityMaxim[1]);
        nsigmatpcKa = generateRandomX(nsigmatpcKlimPurityMaxim[0], nsigmatpcKlimPurityMaxim[1]);
        nsigmatpcPi = generateRandomX(nsigmatpcPilimPurityMaxim[0], nsigmatpcPilimPurityMaxim[1]);
        nsigmatpcPr = generateRandomX(nsigmatpcPrlimPurityMaxim[0], nsigmatpcPrlimPurityMaxim[1]);
        proplifetime = generateRandomX(proplifetimelimPurityMaxim[0], proplifetimelimPurityMaxim[1]);
        rejcomp = generateRandomX(rejcomplimPurityMaxim[0], rejcomplimPurityMaxim[1]);
        v0cospa = generateRandomX(v0cospalimPurityMaxim[0], v0cospalimPurityMaxim[1]);
        v0radius = generateRandomX(v0radiuslimPurityMaxim[0], v0radiuslimPurityMaxim[1]);
        mintpccrrows = generateRandomX(tpccrrowslimPurityMaxim[0], tpccrrowslimPurityMaxim[1]);

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
