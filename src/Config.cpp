#include "ccp-2d.hpp"

#include <filesystem>

namespace ccp2d {

std::string trim(std::string s)
{
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string upper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return s;
}

void makeDir(const std::string& path)
{
    std::filesystem::create_directories(path);
}

std::unordered_map<std::string, std::string> readMap(const std::string& name)
{
    std::ifstream in(name);
    if (!in) throw std::runtime_error("cannot open case file: " + name);
    std::unordered_map<std::string, std::string> m;
    std::string line;
    while (std::getline(in, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string k = trim(line.substr(0, eq));
        const std::string v = trim(line.substr(eq + 1));
        if (!k.empty()) m[k] = v;
    }
    return m;
}

namespace {

bool isAbsolutePath(const std::string& path)
{
    if (path.empty()) return false;
    if (path[0] == '/' || path[0] == '\\') return true;
    return path.size() > 1 && path[1] == ':';
}

std::string parentPath(const std::string& path)
{
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return ".";
    if (pos == 0) return path.substr(0, 1);
    return path.substr(0, pos);
}

std::string joinPath(const std::string& base, const std::string& rel)
{
    if (base.empty() || base == ".") return rel;
    const char last = base.back();
    if (last == '/' || last == '\\') return base + rel;
    return base + "/" + rel;
}

std::string resolveCasePath(const std::string& caseDir, const std::string& path)
{
    if (path.empty() || isAbsolutePath(path)) return path;
    return joinPath(caseDir, path);
}

bool parseBool(std::string s)
{
    s = upper(trim(std::move(s)));
    if (s == "YES" || s == "TRUE" || s == "ON" || s == "1") return true;
    if (s == "NO" || s == "FALSE" || s == "OFF" || s == "0") return false;
    throw std::runtime_error("invalid boolean value: " + s);
}

std::vector<scalar> parseScalarList(std::string s)
{
    for (char& ch : s) {
        if (ch == ',' || ch == ';' || ch == '[' || ch == ']' || ch == '(' || ch == ')') ch = ' ';
    }
    std::vector<scalar> out;
    std::istringstream is(s);
    scalar v = 0.0;
    while (is >> v) out.push_back(v);
    return out;
}


}

Config readConfig(const std::string& name)
{
    Config c;
    const auto m = readMap(name);
    const std::string caseDir = parentPath(name);
    auto has = [&](const std::string& k) { return m.find(k) != m.end(); };
    auto gs = [&](const std::string& k) { return m.at(k); };
    auto gi = [&](const std::string& k) { return std::stoi(m.at(k)); };
    auto gd = [&](const std::string& k) { return std::stod(m.at(k)); };


    if (has("ionCSV")) c.ionCSV = gs("ionCSV");
    if (has("excCSV")) c.excCSV = gs("excCSV");
    if (has("outputDir")) c.outputDir = gs("outputDir");
    if (has("meshFile")) c.meshFile = gs("meshFile");
    if (has("bccFile")) c.bccFile = gs("bccFile");
    if (!c.ionCSV.empty()) c.ionCSV = resolveCasePath(caseDir, c.ionCSV);
    if (!c.excCSV.empty()) c.excCSV = resolveCasePath(caseDir, c.excCSV);
    c.outputDir = resolveCasePath(caseDir, c.outputDir);
    c.meshFile = resolveCasePath(caseDir, c.meshFile);
    c.bccFile = resolveCasePath(caseDir, c.bccFile);

    if (has("cD")) c.cD = gd("cD");
    if (has("epsImplicit")) c.epsImplicit = gd("epsImplicit");
    if (has("chemRowSumDamping")) c.chemRowSumDamping = gd("chemRowSumDamping");
    if (has("omegaImplicit")) c.omegaImplicit = gd("omegaImplicit");
    if (has("dQmaxRel")) c.dQmaxRel = gd("dQmaxRel");
    if (has("qFloor")) c.qFloor = gd("qFloor");
    if (has("relFloor")) c.relFloor = gd("relFloor");
    if (has("TeMin")) c.TeMin = gd("TeMin");
    if (has("TeMax")) c.TeMax = gd("TeMax");

    if (has("PhiRF")) c.PhiRF = gd("PhiRF");
    if (has("PhiAR")) c.PhiAR = gd("PhiAR");
    if (has("fRF")) c.fRF = gd("fRF");
    if (has("phase")) c.phase = gd("phase");
    if (has("Ne0")) c.Ne0 = gd("Ne0");
    if (has("Ni0")) c.Ni0 = gd("Ni0");
    if (has("Te0")) c.Te0 = gd("Te0");
    if (has("Phi0")) c.Phi0 = gd("Phi0");

    if (has("LRef")) c.LRef = gd("LRef");
    if (has("TeRef")) c.TeRef = gd("TeRef");
    if (has("nRef")) c.nRef = gd("nRef");
    if (has("phiRef")) c.phiRef = gd("phiRef");

    if (has("kB")) c.kB = gd("kB");
    if (has("e")) c.e = gd("e");
    if (has("eps0")) c.eps0 = gd("eps0");
    if (has("eV2J")) c.eV2J = gd("eV2J");
    if (has("Hion")) c.Hion = gd("Hion");
    if (has("Hexc")) c.Hexc = gd("Hexc");
    if (has("useExcitationLoss")) c.useExcitationLoss = parseBool(gs("useExcitationLoss"));
    if (has("me")) c.me = gd("me");
    if (has("Gam")) c.Gam = gd("Gam");
    if (has("Ks")) c.Ks = gd("Ks");
    if (has("P")) c.P = gd("P");
    if (has("N0")) c.N0 = gd("N0");
    if (has("De0")) c.De0 = gd("De0");
    if (has("Di0")) c.Di0 = gd("Di0");
    if (has("muE0")) c.muE0 = gd("muE0");
    if (has("muI0")) c.muI0 = gd("muI0");

    if (has("stopPeriod")) c.stopPeriod = gd("stopPeriod");
    if (has("dtPhys")) c.dtPhys = gd("dtPhys");
    if (has("pseudoCFL")) c.pseudoCFL = gd("pseudoCFL");
    if (has("nInnerDT")) c.nInnerDT = gi("nInnerDT");
    if (has("innerResTol")) c.innerResTol = gd("innerResTol");
    if (has("innerRelResTol")) c.innerRelResTol = gd("innerRelResTol");
    if (has("writeInnerResidual")) c.writeInnerResidual = parseBool(gs("writeInnerResidual"));
    if (has("innerResWriteStep")) c.innerResWriteStep = gi("innerResWriteStep");
    if (has("nTIDT")) c.nTIDT = gi("nTIDT");
    if (has("dtOutputPhases")) c.dtOutputPhases = parseScalarList(gs("dtOutputPhases"));
    if (has("printStep")) c.printStep = gi("printStep");
    if (has("resWriteStep")) c.resWriteStep = gi("resWriteStep");
    if (has("writeStep")) c.writeStep = gi("writeStep");

    if (c.innerResTol < 0.0 || c.innerRelResTol < 0.0) {
        throw std::runtime_error("inner residual tolerances must be non-negative");
    }
    if (c.innerResWriteStep < 0) {
        throw std::runtime_error("innerResWriteStep must be non-negative");
    }
    if (c.ionCSV.empty()) {
        throw std::runtime_error("ionCSV is required");
    }
    if (!has("Hion")) {
        throw std::runtime_error("Hion is required");
    }
    if (c.useExcitationLoss && c.excCSV.empty()) {
        throw std::runtime_error("useExcitationLoss=yes requires excCSV");
    }
    if (c.useExcitationLoss && !has("Hexc")) {
        throw std::runtime_error("useExcitationLoss=yes requires Hexc");
    }
    if (c.stopPeriod <= 0.0) throw std::runtime_error("stopPeriod must be positive");
    if (c.dtPhys <= 0.0) throw std::runtime_error("dtPhys must be positive");
    if (c.pseudoCFL <= 0.0) throw std::runtime_error("pseudoCFL must be positive");
    if (c.nInnerDT <= 0) throw std::runtime_error("nInnerDT must be positive");
    if (c.nTIDT <= 0) throw std::runtime_error("nTIDT must be positive");
    return c;
}

void nondimensionalize(Config& c)
{
    c.N = c.N0 * c.P;
    c.De = c.De0 / c.N;
    c.Di = c.Di0 / c.N;
    c.muE = c.muE0 / c.N;
    c.muI = c.muI0 / c.N;
    c.Hion *= c.eV2J;
    c.Hexc *= c.eV2J;

    const scalar vth = std::sqrt(8.0 * c.kB * c.TeRef / (pi * c.me));
    c.tRef = c.LRef / vth;
    c.fRef = 1.0 / c.tRef;
    c.EeRef = 1.5 * c.kB * c.nRef * c.TeRef;
    c.DeRef = c.LRef * vth;
    c.DiRef = c.LRef * vth;
    c.muERef = c.DeRef / c.phiRef;
    c.muIRef = c.DiRef / c.phiRef;
    c.eRef = c.EeRef / (c.phiRef * c.nRef);
    c.klRef = 1.0 / (c.nRef * c.tRef);
    c.HRef = c.EeRef / c.nRef;
    c.epsRef = c.EeRef * c.LRef * c.LRef / (c.phiRef * c.phiRef);

    c.e /= c.eRef;
    c.eps0 /= c.epsRef;
    c.N /= c.nRef;
    c.Hion /= c.HRef;
    c.Hexc /= c.HRef;
    c.De /= c.DeRef;
    c.Di /= c.DiRef;
    c.muE /= c.muERef;
    c.muI /= c.muIRef;
    c.Ks /= vth;
    c.PhiRF /= c.phiRef;
    c.PhiAR /= c.phiRef;
    c.fRF /= c.fRef;
    c.period = 1.0 / c.fRF;
    c.dtPhys /= c.tRef;
    c.phase = c.phase * pi / 180.0;
    c.Ne0 /= c.nRef;
    c.Ni0 /= c.nRef;
    c.Te0 /= c.TeRef;
    c.TeMin /= c.TeRef;
    c.TeMax /= c.TeRef;
    if (std::isfinite(c.TeMax) && c.TeMin > c.TeMax) {
        throw std::runtime_error("TeMin must be <= TeMax after nondimensionalization");
    }
    c.Phi0 /= c.phiRef;
}

void printNondimensionalizationInfo(const Config& c, std::ostream& os)
{
    const auto oldFlags = os.flags();
    const auto oldPrecision = os.precision();

    os << std::scientific << std::setprecision(6);

    auto refRow = [&](const std::string& name, scalar value, const std::string& unit) {
        os << "  " << std::left << std::setw(10) << name
           << " = " << std::right << std::setw(14) << value
           << "  [" << unit << "]\n";
    };

    auto quantityRow = [&](const std::string& name,
                           scalar dimValue,
                           scalar refValue,
                           scalar nondimValue,
                           const std::string& unit) {
        os << "  " << std::left << std::setw(10) << name
           << "  dim = " << std::right << std::setw(14) << dimValue
           << "  ref = " << std::setw(14) << refValue
           << "  nondim = " << std::setw(14) << nondimValue
           << "  [" << unit << "]\n";
    };

    os << "[Nondimensionalization]\n";
    os << "  Reference values\n";
    refRow("LRef", c.LRef, "m");
    refRow("tRef", c.tRef, "s");
    refRow("fRef", c.fRef, "Hz");
    refRow("nRef", c.nRef, "m^-3");
    refRow("TeRef", c.TeRef, "K");
    refRow("phiRef", c.phiRef, "V");
    refRow("EeRef", c.EeRef, "J/m^3");
    refRow("DeRef", c.DeRef, "m^2/s");
    refRow("DiRef", c.DiRef, "m^2/s");
    refRow("muERef", c.muERef, "m^2/(V s)");
    refRow("muIRef", c.muIRef, "m^2/(V s)");
    refRow("eRef", c.eRef, "C");
    refRow("klRef", c.klRef, "m^3/s");
    refRow("HRef", c.HRef, "J");
    refRow("epsRef", c.epsRef, "F/m");

    os << "  Dimensional inputs after scaling check\n";
    quantityRow("Ne0", c.Ne0 * c.nRef, c.nRef, c.Ne0, "m^-3");
    quantityRow("Ni0", c.Ni0 * c.nRef, c.nRef, c.Ni0, "m^-3");
    quantityRow("Te0", c.Te0 * c.TeRef, c.TeRef, c.Te0, "K");
    if (c.TeMin > 0.0) quantityRow("TeMin", c.TeMin * c.TeRef, c.TeRef, c.TeMin, "K");
    if (std::isfinite(c.TeMax)) quantityRow("TeMax", c.TeMax * c.TeRef, c.TeRef, c.TeMax, "K");
    quantityRow("Phi0", c.Phi0 * c.phiRef, c.phiRef, c.Phi0, "V");
    quantityRow("PhiRF", c.PhiRF * c.phiRef, c.phiRef, c.PhiRF, "V");
    quantityRow("PhiAR", c.PhiAR * c.phiRef, c.phiRef, c.PhiAR, "V");
    quantityRow("fRF", c.fRF * c.fRef, c.fRef, c.fRF, "Hz");
    quantityRow("period", c.period * c.tRef, c.tRef, c.period, "s");
    quantityRow("dtPhys", c.dtPhys * c.tRef, c.tRef, c.dtPhys, "s");
    quantityRow("N", c.N * c.nRef, c.nRef, c.N, "m^-3");
    quantityRow("De", c.De * c.DeRef, c.DeRef, c.De, "m^2/s");
    quantityRow("Di", c.Di * c.DiRef, c.DiRef, c.Di, "m^2/s");
    quantityRow("muE", c.muE * c.muERef, c.muERef, c.muE, "m^2/(V s)");
    quantityRow("muI", c.muI * c.muIRef, c.muIRef, c.muI, "m^2/(V s)");
    quantityRow("Ks", c.Ks * c.LRef / c.tRef, c.LRef / c.tRef, c.Ks, "m/s");
    quantityRow("eps0", c.eps0 * c.epsRef, c.epsRef, c.eps0, "F/m");
    quantityRow("e", c.e * c.eRef, c.eRef, c.e, "C");
    quantityRow("Hion", c.Hion * c.HRef, c.HRef, c.Hion, "J");
    quantityRow("Hexc", c.Hexc * c.HRef, c.HRef, c.Hexc, "J");
    quantityRow("phase", c.phase * 180.0 / pi, 1.0, c.phase, "deg/rad");

    os.flags(oldFlags);
    os.precision(oldPrecision);
}

}
