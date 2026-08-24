#include "ccp-2d.hpp"

namespace ccp2d {

Solver::Solver(Config cfg)
    : cfg_(std::move(cfg)),
      chemIon_(cfg_.ionCSV),
      chemExc_(cfg_.useExcitationLoss ? cfg_.excCSV : cfg_.ionCSV),
      dtPhys_(cfg_.dtPhys),
      dtPhysStep_(cfg_.dtPhys)
{
    initGridAndFaces();
    initDTSampling();
    updateTimeDependentBoundaries();
}

void Solver::initDTSampling()
{
    dtSamplePhases_ = cfg_.dtOutputPhases;
    if (dtSamplePhases_.empty()) {
        const label n = std::max<label>(cfg_.nTIDT, 1);
        dtSamplePhases_.reserve(static_cast<size_t>(n));
        for (label k = 0; k < n; ++k) {
            dtSamplePhases_.push_back(scalar(k) / scalar(n));
        }
    }
    for (scalar p : dtSamplePhases_) {
        if (p < -1.0e-12 || p > 1.0 + 1.0e-12) {
            throw std::runtime_error("dtOutputPhases values must be in [0, 1]");
        }
    }
    for (scalar& p : dtSamplePhases_) {
        p = std::max(scalar(0.0), std::min(scalar(1.0), p));
    }
    std::sort(dtSamplePhases_.begin(), dtSamplePhases_.end());
    dtSamplePhases_.erase(std::unique(dtSamplePhases_.begin(), dtSamplePhases_.end(),
                                      [](scalar a, scalar b) { return std::abs(a - b) < 1.0e-12; }),
                          dtSamplePhases_.end());
    nextAverageCycle_ = 1;
    nextDTSample_ = 0;
}

void Solver::storeOld()
{
    for (Cell& c : cells_) {
        c.NeOld = c.Ne;
        c.NiOld = c.Ni;
        c.EeOld = c.Ee;
        c.PhiOld = c.Phi;
    }
}

void Solver::storePhysicalOld()
{
    for (Cell& c : cells_) {
        c.NePhysOld = c.Ne;
        c.NiPhysOld = c.Ni;
        c.EePhysOld = c.Ee;
        c.PhiPhysOld = c.Phi;
        c.NePhysOlder = c.Ne;
        c.NiPhysOlder = c.Ni;
        c.EePhysOlder = c.Ee;
        c.PhiPhysOlder = c.Phi;
    }
}

void Solver::preparePhysicalHistoryForStep()
{
    for (Cell& c : cells_) {
        c.NePhysOlder = c.NePhysOld;
        c.NiPhysOlder = c.NiPhysOld;
        c.EePhysOlder = c.EePhysOld;
        c.PhiPhysOlder = c.PhiPhysOld;
        c.NePhysOld = c.Ne;
        c.NiPhysOld = c.Ni;
        c.EePhysOld = c.Ee;
        c.PhiPhysOld = c.Phi;
    }
}

bool Solver::useBDF2DT() const
{
    if (physicalStep_ <= 0) return false;
    const scalar rel = std::abs(dtPhysStep_ - dtPhys_) / std::max(dtPhys_, tiny);
    return rel < 1.0e-10;
}

scalar Solver::physicalTimeDerivativeCoeff() const
{
    return (useBDF2DT() ? 1.5 : 1.0) / std::max(dtPhysStep_, tiny);
}

void Solver::updateTimeDependentBoundaries()
{
    auto updateFace = [&](Face& f) {
        if (!f.timeDependentPhi) return;
        f.phiBC(0) = f.phiAmp * std::sin(2.0 * pi * cfg_.fRF * runTime_ + f.phiPhase);
    };
    for (Face& f : xFaces_) updateFace(f);
    for (Face& f : yFaces_) updateFace(f);
}

void Solver::run()
{
    auto t0 = std::chrono::high_resolution_clock::now();
    {
        std::ofstream os(cfg_.outputDir + "/stop.dat", std::ios::trunc);
        os << "0\n";
    }

    std::cout << "\nStarting dual-time BDF2 implicit iteration..."
              << "\n  Poisson    : semi-implicit"
              << "\n  cells      : " << nCells_
              << "\n  stopPeriod : " << cfg_.stopPeriod
              << "\n  dtPhys     : " << std::scientific << std::setprecision(6)
              << dtPhys_ * cfg_.tRef << " s"
              << std::fixed << std::setprecision(6)
              << "\n  pseudoCFL  : " << cfg_.pseudoCFL
              << "\n  nInnerDT   : " << cfg_.nInnerDT
              << "\n\n";

    updateTimeDependentBoundaries();
    computeDt();
    solvePoisson();
    computeFluxesAD();
    storeOld();
    storePhysicalOld();
    writeAverageNeNiEeDT(0.0, 0, 1.0);
    writeDueDTFields(runCyc_);

    const scalar stopTime = cfg_.stopPeriod * cfg_.period;
    const scalar stopTolerance = 1.0e-10 * std::max(stopTime, cfg_.period);
    for (step_ = 0; runTime_ + stopTolerance < stopTime; ++step_) {
        {
            std::ifstream is(cfg_.outputDir + "/stop.dat");
            int stop = 0;
            if (is && (is >> stop) && stop == 1) {
                std::cout << "Stop signal received.\n";
                break;
            }
        }

        const scalar cycleBefore = runCyc_;
        dtPhysStep_ = std::min(dtPhys_, std::max(stopTime - runTime_, scalar(0.0)));
        if (dtPhysStep_ <= tiny) break;

        preparePhysicalHistoryForStep();
        runTime_ += dtPhysStep_;
        if (stopTime - runTime_ < stopTolerance) runTime_ = stopTime;
        runCyc_ = runTime_ / cfg_.period;
        updateTimeDependentBoundaries();

        std::array<scalar, 4> rr{{0.0, 0.0, 0.0, 0.0}};
        std::array<scalar, 5> g0{{1.0, 1.0, 1.0, 1.0, 1.0}};
        bool haveG0 = false;
        const bool writeInnerThisStep = shouldWriteInnerResidual();
        const bool trackInnerResidual = writeInnerThisStep || cfg_.innerRelResTol > 0.0;
        std::vector<InnerResidualSample> innerResidualRows;
        label lastInnerStored = -1;
        label usedInner = 0;

        for (label inner = 0; inner < cfg_.nInnerDT; ++inner) {
            solvePoisson();
            computeFluxesAD();
            computeDt();
            updateChemistryJacobian();
            computeJouleJacobianAD();
            assembleRHS_DT();
            rr = residualLog10();
            bool converged = innerResidualConverged(rr);
            if (trackInnerResidual) {
                const auto g = residualNorms();
                if (!haveG0) {
                    g0 = g;
                    haveG0 = true;
                }
                if (writeInnerThisStep) {
                    innerResidualRows.push_back(InnerResidualSample{inner, dtMin_, g, g0});
                    lastInnerStored = inner;
                }
                converged = converged || innerRelativeResidualConverged(g, g0);
            }
            if (converged) {
                usedInner = inner;
                break;
            }
            storeOld();
            forwardSweepImplicit();
            backwardSweepImplicit();
            updateSolution();
            usedInner = inner + 1;
        }
        innerIter_ = usedInner;

        solvePoisson();
        computeFluxesAD();
        computeDt();
        updateChemistryJacobian();
        computeJouleJacobianAD();
        assembleRHS_DT();
        rr = residualLog10();

        if (writeInnerThisStep && usedInner != lastInnerStored) {
            const auto g = residualNorms();
            if (!haveG0) g0 = g;
            innerResidualRows.push_back(InnerResidualSample{usedInner, dtMin_, g, g0});
        }
        if (writeInnerThisStep && !innerResidualRows.empty()) {
            writeInnerResidualZone(innerResidualRows);
        }

        const auto t1 = std::chrono::high_resolution_clock::now();
        const scalar wall = std::chrono::duration<scalar>(t1 - t0).count();
        writeResidual(rr, wall);

        const bool finalStep = runCyc_ >= cfg_.stopPeriod;
        if (step_ % cfg_.printStep == 0 || finalStep) infoResiduals(rr, wall);
        if (!allFinite()) {
            throw std::runtime_error("non-finite solution at physical step " + std::to_string(step_));
        }

        writeIterSolution();
        writeFinalSolution();
        writeRestart();
        writeDueDTAverages(cycleBefore, wall);
        writeDueDTFields(cycleBefore);
        ++physicalStep_;
    }

    writeIterSolution(true);
    writeFinalSolution(true);
    writeRestart(true);
    std::cout << "Done. Total physical steps = " << physicalStep_ << "\n";
}

}
