// Read in all related dictionaries

IOdictionary APOLLOControlDict
(
    IOobject
    (
		"APOLLOControlDict",
		runTime.system(),
		runTime,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

const scalar& maxSweep = readScalar(APOLLOControlDict.lookup("maxSweep"));
const scalar& currentTolerance = readScalar(APOLLOControlDict.lookup("currentTolerance"));
const scalar& protonSweep = readScalar(APOLLOControlDict.lookup("protonSweep"));
const scalar& electronSweep = readScalar(APOLLOControlDict.lookup("electronSweep"));
const scalar& THV = readScalar(APOLLOControlDict.lookup("THV")) ;
const scalar& eta_blockl_guess = readScalar(APOLLOControlDict.lookup("etaBlockl"));
const scalar& eta_blockr_guess = readScalar(APOLLOControlDict.lookup("etaBlockr"));
const word& apolloSolverMethod = word(APOLLOControlDict.lookup("solverMethod"));
const word& apolloEqnSystem	= word(APOLLOControlDict.lookup("equationSystem"));
const scalar& implicitFactor = readScalar(APOLLOControlDict.lookup("implicitFactor"));
