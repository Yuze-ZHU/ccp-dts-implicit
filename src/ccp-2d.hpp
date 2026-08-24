#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Sacado.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ccp2d {

using scalar = double;
using label = int;
using ADScalar = Sacado::Fad::DFad<scalar>;
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;
using Mat4 = Eigen::Matrix<scalar, 4, 4, Eigen::RowMajor>;
using Mat3 = Eigen::Matrix<scalar, 3, 3, Eigen::RowMajor>;
using Vec3 = Eigen::Matrix<scalar, 3, 1>;

inline constexpr scalar pi = 3.141592653589793238462643383279502884;
inline constexpr scalar tiny = 1.0e-120;

std::string trim(std::string s);
std::string upper(std::string s);
void makeDir(const std::string& path);

struct Config {
    std::string ionCSV;
    std::string excCSV;
    std::string outputDir = "results_ad";
    std::string meshFile;
    std::string bccFile;

    label nCellX = 0;
    label nCellY = 0;
    label nNodeX = 0;
    label nNodeY = 0;
    static constexpr label numT = 1;
    scalar cD = 2.0;
    scalar epsImplicit = 0.7;
    scalar chemRowSumDamping = 1.0;
    scalar omegaImplicit = 0.7;
    scalar dQmaxRel = 0.3;
    scalar qFloor = 1.0e-3;
    scalar relFloor = 1.0e-3;
    scalar TeMin = 0.0;
    scalar TeMax = std::numeric_limits<scalar>::infinity();

    scalar PhiRF = 100.0;
    scalar PhiAR = 0.0;
    scalar fRF = 13.56e6;
    scalar phase = 90.0;
    scalar period = 0.0;

    scalar Ni0 = 1.0e14;
    scalar Ne0 = 1.0e14;
    scalar Te0 = 11604.0;
    scalar Phi0 = 0.0;

    scalar LRef = 2.54e-2;
    scalar TeRef = 11604.0;
    scalar nRef = 1.0e14;
    scalar phiRef = 100.0;
    scalar tRef = 0.0;
    scalar fRef = 0.0;
    scalar EeRef = 0.0;
    scalar DeRef = 0.0;
    scalar DiRef = 0.0;
    scalar muERef = 0.0;
    scalar muIRef = 0.0;
    scalar eRef = 0.0;
    scalar klRef = 0.0;
    scalar HRef = 0.0;
    scalar epsRef = 0.0;

    scalar kB = 1.380649e-23;
    scalar e = 1.6022e-19;
    scalar eps0 = 8.8542e-12;
    scalar eV2J = 1.602e-19;
    scalar Hion = 15.7;
    scalar Hexc = 11.56;
    bool useExcitationLoss = false;
    scalar me = 9.10938356e-31;
    scalar Gam = 0.01;
    scalar Ks = 0.0;
    scalar P = 1.0;
    scalar N0 = 3.22e22;
    scalar De0 = 3.86e24;
    scalar Di0 = 2.07e20;
    scalar muE0 = 9.66e23;
    scalar muI0 = 4.65e21;
    scalar N = 0.0;
    scalar De = 0.0;
    scalar Di = 0.0;
    scalar muE = 0.0;
    scalar muI = 0.0;

    scalar stopPeriod = 1.0;
    scalar dtPhys = 0.0;
    scalar pseudoCFL = 10000.0;
    label nInnerDT = 20;
    scalar innerResTol = 0.0;
    scalar innerRelResTol = 0.0;
    bool writeInnerResidual = false;
    label innerResWriteStep = 1;
    label nTIDT = 20;
    std::vector<scalar> dtOutputPhases;
    label printStep = 50;
    label resWriteStep = 10;
    label writeStep = 200;
};

std::unordered_map<std::string, std::string> readMap(const std::string& name);
Config readConfig(const std::string& name);
void nondimensionalize(Config& c);
void printNondimensionalizationInfo(const Config& c, std::ostream& os);

struct ChemPoint {
    scalar energy = 0.0;
    scalar rate = 0.0;
};

struct InnerResidualSample {
    label inner = 0;
    scalar dtau = 0.0;
    std::array<scalar, 5> g{};
    std::array<scalar, 5> g0{};
};

class ChemTable {
public:
    explicit ChemTable(const std::string& name);
    scalar interpolate(scalar eV) const;

private:
    std::vector<ChemPoint> rows_;
};

struct Flux {
    scalar ne = 0.0;
    scalar ni = 0.0;
    scalar ee = 0.0;
    scalar j = 0.0;
    scalar phi = 0.0;
    scalar ef = 0.0;
};

template <class T>
struct FluxResult {
    T ne, ni, ee, j, phi, ef;
};

enum class WallModel { None, Sakiyama, Lymberopoulos };

inline scalar getVal(const scalar& x) { return x; }
inline scalar getVal(const ADScalar& x) { return x.val(); }

template <class T>
T safeMax(const T& x, scalar floor)
{
    return getVal(x) < floor ? T(floor) : x;
}

template <class T>
FluxResult<T> internalFluxKernel(
    const T& NeL, const T& NiL, const T& EeL, const T& PhiL,
    const T& NeR, const T& NiR, const T& EeR, const T& PhiR,
    scalar dist, const Config& c)
{
    const T dNe = (NeR - NeL) / dist;
    const T dNi = (NiR - NiL) / dist;
    const T dPhi = (PhiR - PhiL) / dist;
    const T TeL = EeL / NeL;
    const T TeR = EeR / NeR;
    const T dTe = (TeR - TeL) / dist;
    const bool upwind = getVal(PhiR) >= getVal(PhiL);
    const T NeC = upwind ? NeL : NeR;
    const T NiC = upwind ? NiR : NiL;
    const T TeC = upwind ? TeL : TeR;
    const T PhiC = 0.5 * (PhiL + PhiR);

    FluxResult<T> f;
    f.ef = -dPhi;
    f.ne = -c.De * dNe + c.muE * NeC * dPhi;
    f.ni = -c.Di * dNi - c.muI * NiC * dPhi;
    f.phi = -f.ef * c.eps0;
    const T qE = -(5.0 / 3.0) * c.De * NeC * dTe + (5.0 / 3.0) * TeC * f.ne;
    f.j = c.e * f.ne;
    f.ee = qE - c.e * PhiC * f.ne;
    return f;
}

template <class T>
FluxResult<T> minusBCFluxKernel(
    const T& NeR, const T& NiR, const T& EeR, const T& PhiR,
    scalar PhiLbc, WallModel wallModel, scalar wallTe, scalar dist, const Config& c)
{
    using Sacado::Fad::sqrt;
    const T dPhi = (PhiR - PhiLbc) / dist;
    const T TeR = EeR / NeR;

    FluxResult<T> f;
    f.ef = -dPhi;
    f.ni = getVal(f.ef) <= 0.0 ? c.muI * NiR * f.ef : T(0.0);
    f.phi = -f.ef * c.eps0;
    if (wallModel == WallModel::Lymberopoulos) {
        const T dTe = (TeR - wallTe) / dist;
        f.ne = -c.Ks * NeR - c.Gam * f.ni;
        const T qE = -(5.0 / 3.0) * c.De * NeR * dTe + (5.0 / 3.0) * TeR * f.ne;
        f.j = c.e * f.ne;
        f.ee = qE - c.e * PhiR * f.ne;
        return f;
    }
    const T k = 0.25 * sqrt(TeR);
    f.ne = -k * NeR - c.Gam * f.ni;
    f.j = c.e * f.ne;
    f.ee = (5.0 / 3.0) * TeR * f.ne - c.e * PhiR * f.ne;
    return f;
}

template <class T>
FluxResult<T> plusBCFluxKernel(
    const T& NeL, const T& NiL, const T& EeL, const T& PhiL,
    scalar PhiRbc, WallModel wallModel, scalar wallTe, scalar dist, const Config& c)
{
    using Sacado::Fad::sqrt;
    const T dPhi = (PhiRbc - PhiL) / dist;
    const T TeL = EeL / NeL;

    FluxResult<T> f;
    f.ef = -dPhi;
    f.ni = getVal(f.ef) >= 0.0 ? c.muI * NiL * f.ef : T(0.0);
    f.phi = -f.ef * c.eps0;
    if (wallModel == WallModel::Lymberopoulos) {
        const T dTe = (wallTe - TeL) / dist;
        f.ne = c.Ks * NeL - c.Gam * f.ni;
        const T qE = -(5.0 / 3.0) * c.De * NeL * dTe + (5.0 / 3.0) * TeL * f.ne;
        f.j = c.e * f.ne;
        f.ee = qE - c.e * PhiL * f.ne;
        return f;
    }
    const T k = 0.25 * sqrt(TeL);
    f.ne = k * NeL - c.Gam * f.ni;
    f.j = c.e * f.ne;
    f.ee = (5.0 / 3.0) * TeL * f.ne - c.e * PhiL * f.ne;
    return f;
}

struct Cell {
    scalar x = 0.0;
    scalar y = 0.0;
    scalar vol = 0.0;
    scalar dt = 1.0e-6;

    Vec Ne, Ni, Ee, Te, Phi, Ex, Ey, kl, kIon, kExc;
    Vec NeOld, NiOld, EeOld, PhiOld;
    Vec NePhysOld, NiPhysOld, EePhysOld, PhiPhysOld;
    Vec NePhysOlder, NiPhysOlder, EePhysOlder, PhiPhysOlder;
    Vec rhsNe, rhsNi, rhsEe;
    Vec resNe, resNi, resEe, resPhi;
    Vec dQhNe, dQhNi, dQhEe;
    Vec dQNe, dQNi, dQEe;
    std::vector<Mat4> chemJ;
    std::vector<Mat4> jouleC, jouleL, jouleR, jouleB, jouleT;

    explicit Cell(label nt = 1)
        : Ne(Vec::Zero(nt)), Ni(Vec::Zero(nt)), Ee(Vec::Zero(nt)), Te(Vec::Zero(nt)),
          Phi(Vec::Zero(nt)), Ex(Vec::Zero(nt)), Ey(Vec::Zero(nt)), kl(Vec::Zero(nt)),
          kIon(Vec::Zero(nt)), kExc(Vec::Zero(nt)),
          NeOld(Vec::Zero(nt)), NiOld(Vec::Zero(nt)), EeOld(Vec::Zero(nt)),
          PhiOld(Vec::Zero(nt)),
          NePhysOld(Vec::Zero(nt)), NiPhysOld(Vec::Zero(nt)), EePhysOld(Vec::Zero(nt)),
          PhiPhysOld(Vec::Zero(nt)),
          NePhysOlder(Vec::Zero(nt)), NiPhysOlder(Vec::Zero(nt)), EePhysOlder(Vec::Zero(nt)),
          PhiPhysOlder(Vec::Zero(nt)),
          rhsNe(Vec::Zero(nt)), rhsNi(Vec::Zero(nt)), rhsEe(Vec::Zero(nt)),
          resNe(Vec::Zero(nt)), resNi(Vec::Zero(nt)), resEe(Vec::Zero(nt)),
          resPhi(Vec::Zero(nt)), dQhNe(Vec::Zero(nt)), dQhNi(Vec::Zero(nt)),
          dQhEe(Vec::Zero(nt)), dQNe(Vec::Zero(nt)), dQNi(Vec::Zero(nt)),
          dQEe(Vec::Zero(nt)), chemJ(nt, Mat4::Zero()),
          jouleC(nt, Mat4::Zero()), jouleL(nt, Mat4::Zero()),
          jouleR(nt, Mat4::Zero()), jouleB(nt, Mat4::Zero()),
          jouleT(nt, Mat4::Zero())
    {}
};

enum class FaceKind { Interior, MinusBC, PlusBC, Symmetry };

struct Node {
    scalar x = 0.0;
    scalar y = 0.0;
};

struct Face {
    FaceKind kind = FaceKind::Interior;
    label cL = -1;
    label cR = -1;
    bool bcSet = false;
    bool timeDependentPhi = false;
    scalar phiAmp = 0.0;
    scalar phiPhase = 0.0;
    WallModel wallModel = WallModel::None;
    scalar wallTe = 0.0;
    scalar area = 1.0;
    scalar dist = 1.0;
    Vec phiBC;
    std::vector<Flux> flux;
    std::vector<Mat4> JL;
    std::vector<Mat4> JR;

    Face() = default;
    explicit Face(label nt) : phiBC(Vec::Zero(nt)), flux(nt), JL(nt, Mat4::Zero()), JR(nt, Mat4::Zero()) {}
};

scalar dist2D(const Node& a, const Node& b);
scalar quadArea(const Node& a, const Node& b, const Node& c, const Node& d);

class Solver {
public:
    explicit Solver(Config cfg);
    void run();

private:
    Config cfg_;
    ChemTable chemIon_;
    ChemTable chemExc_;
    std::vector<Cell> cells_;
    std::vector<Face> xFaces_;
    std::vector<Face> yFaces_;
    std::vector<Node> nodes_;
    label nCells_ = 0;
    label nNodes_ = 0;
    scalar dtMin_ = 0.0;
    scalar dtPhys_ = 0.0;
    scalar dtPhysStep_ = 0.0;
    label step_ = 0;
    label innerIter_ = 0;
    label physicalStep_ = 0;
    scalar runTime_ = 0.0;
    scalar runCyc_ = 0.0;
    label nextAverageCycle_ = 0;
    label nextDTSample_ = 0;
    std::vector<scalar> dtSamplePhases_;
    bool dtFieldWritten_ = false;

    label nid(label i, label j) const { return i + cfg_.nNodeX * j; }
    label cid(label i, label j) const { return i + cfg_.nCellX * j; }
    label ci(label cID) const { return cID % cfg_.nCellX; }
    label cj(label cID) const { return cID / cfg_.nCellX; }
    label xf(label i, label j) const { return i + (cfg_.nCellX + 1) * j; }
    label yf(label i, label j) const { return i + cfg_.nCellX * j; }

    void initDTSampling();
    void storeOld();
    void storePhysicalOld();
    void preparePhysicalHistoryForStep();
    bool useBDF2DT() const;
    scalar physicalTimeDerivativeCoeff() const;

    bool readPlot3DMesh(const std::string& filename);
    void initGridAndFaces();
    void applyBoundaryFile();
    void updateTimeDependentBoundaries();
    Node cellCenter(label i, label j) const;
    Node faceCenterX(label i, label j) const;
    Node faceCenterY(label i, label j) const;

    void updateChemistryJacobian();
    void computeFaceAD(Face& f, label t);
    void computeFluxesAD();
    void addJRow(Mat4& M, scalar scale, const Mat4& J, label t);
    void computeJouleJacobianAD();

    void assembleRHS_DT();
    Mat3 buildBlockB(label cID, label t) const;
    void forwardSweepImplicit();
    void backwardSweepImplicit();
    void updateSolution();
    void computeDt();
    bool innerResidualConverged(const std::array<scalar, 4>& rr) const;
    bool innerRelativeResidualConverged(const std::array<scalar, 5>& g,
                                        const std::array<scalar, 5>& g0) const;

    scalar semiImplicitPoissonFaceCoeff(label c, label nb, label t) const;
    scalar semiImplicitPoissonBoundaryCoeff(label c, label t) const;
    scalar poissonTimeStep(const Cell& c) const;
    bool hasPoissonDirichletBoundary() const;
    void assemblePhysicalPoisson(label t, std::vector<Eigen::Triplet<scalar>>& triplets, Vec& b) const;
    void addSemiImplicitPoissonCorrection(label t, std::vector<Eigen::Triplet<scalar>>& triplets, Vec& b) const;
    void solvePoisson();
    void updateElectricField();
    scalar poissonResidual(label cID, label t) const;

    std::array<scalar, 4> residualLog10() const;
    std::array<scalar, 5> residualNorms() const;
    void writeResidual(const std::array<scalar, 4>& rr, scalar wall);
    bool shouldWriteInnerResidual() const;
    void writeInnerResidualZone(const std::vector<InnerResidualSample>& rows) const;
    void infoResiduals(const std::array<scalar, 4>& rr, scalar wall) const;
    scalar cellVal(const Cell& c, label t, label v) const;
    scalar minVar(label v) const;
    scalar maxVar(label v) const;
    bool allFinite() const;

    template <class F>
    void writeNodeBlock(std::ostream& out, label n, F valueOf) const
    {
        for (label k = 0; k < n; ++k) {
            out << valueOf(k) << ((k + 1) % 6 == 0 ? "\n" : " ");
        }
        if (n % 6 != 0) out << "\n";
    }

    template <class F>
    void writeCellBlock(std::ostream& out, F valueOf) const
    {
        label k = 0;
        for (label j = 0; j < cfg_.nCellY; ++j) {
            for (label i = 0; i < cfg_.nCellX; ++i) {
                const label cID = cid(i, j);
                out << valueOf(cells_[cID], cID)
                    << ((++k) % 6 == 0 ? "\n" : " ");
            }
        }
        if (k % 6 != 0) out << "\n";
    }

    void writeConnectivity(std::ostream& out) const;
    void writeMeshGeometry() const;
    void writeFieldTecplotMergedAllIT(const std::string& fname,
                                      std::vector<label> iTList,
                                      bool append,
                                      const std::string& zoneTag) const;
    void writeIterSolution(bool force = false) const;
    void writeFinalSolution(bool force = false) const;
    void writeRestart(bool force = false) const;
    void writeAverageNeNiEeDT(scalar wall, label targetCycle, scalar theta) const;
    void writeUnsteadyFlowFieldDT(label sampleID, scalar phase, scalar targetCycle,
                                  scalar theta, bool append) const;
    void writeDueDTAverages(scalar cycleBefore, scalar wall);
    void writeDueDTFields(scalar cycleBefore);
};

int runPlasma2D(int argc, char** argv);

}
