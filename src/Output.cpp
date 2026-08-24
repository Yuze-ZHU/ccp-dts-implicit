#include "ccp-2d.hpp"

namespace ccp2d {

std::array<scalar, 4> Solver::residualLog10() const
{
    const auto g = residualNorms();
    auto lg = [](scalar v) { return std::log10(std::max(v, 1.0e-30)); };
    return {lg(g[1]), lg(g[2]), lg(g[3]), lg(g[4])};
}

std::array<scalar, 5> Solver::residualNorms() const
{
    scalar sNe = 0.0;
    scalar sNi = 0.0;
    scalar sEe = 0.0;
    scalar sPhi = 0.0;
    const scalar dof = static_cast<scalar>(std::max(nCells_ * cfg_.numT, 1));
    for (const Cell& c : cells_) {
        for (label t = 0; t < cfg_.numT; ++t) {
            sNe += c.resNe(t) * c.resNe(t);
            sNi += c.resNi(t) * c.resNi(t);
            sEe += c.resEe(t) * c.resEe(t);
            sPhi += c.resPhi(t) * c.resPhi(t);
        }
    }
    const scalar gNe = std::sqrt(sNe / dof);
    const scalar gNi = std::sqrt(sNi / dof);
    const scalar gEe = std::sqrt(sEe / dof);
    const scalar gPhi = std::sqrt(sPhi / dof);
    const scalar gFluid = std::sqrt((sNe + sNi + sEe) / (3.0 * dof));
    return {gFluid, gNe, gNi, gEe, gPhi};
}

void Solver::writeResidual(const std::array<scalar, 4>& rr, scalar wall)
{
    if (cfg_.resWriteStep <= 0 || step_ % cfg_.resWriteStep != 0) return;
    static bool firstCall = true;
    const std::string fname = cfg_.outputDir + "/ResidualHistory.dat";
    std::ofstream out(fname, firstCall ? std::ios::trunc : std::ios::app);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << '\n';
        return;
    }
    if (firstCall) {
        out << "TITLE = \"Dual-time residual history\"\n"
            << "VARIABLES = \"step\", \"inner\", \"wallTime\", \"runTime\", \"runCyc\", "
            << "\"dtPhys\", \"dtau\", \"logGNe\", \"logGNi\", \"logGEe\", \"logGPhi\"\n"
            << "ZONE T=\"Residual\", F=POINT\n";
        firstCall = false;
    }
    out << step_ << '\t'
        << innerIter_ << '\t'
        << std::scientific << std::setprecision(8)
        << wall << '\t'
        << runTime_ * cfg_.tRef << '\t'
        << runCyc_ << '\t'
        << dtPhysStep_ * cfg_.tRef << '\t'
        << dtMin_ * cfg_.tRef << '\t'
        << std::fixed << std::setprecision(8)
        << rr[0] << '\t' << rr[1] << '\t' << rr[2] << '\t' << rr[3] << '\n';
}
void Solver::infoResiduals(const std::array<scalar, 4>& rr, scalar wall) const
{
    std::cout << "step=" << std::setw(7) << step_
              << "  inner=" << innerIter_
              << "  wall=" << std::fixed << std::setprecision(1) << wall << "s"
              << "  runTime=" << std::scientific << std::setprecision(4)
              << runTime_ * cfg_.tRef << "s"
              << "  runCyc=" << std::fixed << std::setprecision(6) << runCyc_
              << "  dtPhys=" << std::scientific << std::setprecision(2)
              << dtPhysStep_ * cfg_.tRef << "s"
              << "  dtau=" << dtMin_ * cfg_.tRef << "s"
              << "  logG=(" << std::fixed << std::setprecision(4)
              << rr[0] << ", " << rr[1] << ", " << rr[2] << ", " << rr[3] << ")"
              << std::endl;
}

bool Solver::shouldWriteInnerResidual() const
{
    return cfg_.writeInnerResidual && cfg_.innerResWriteStep > 0
        && step_ % cfg_.innerResWriteStep == 0;
}

void Solver::writeInnerResidualZone(const std::vector<InnerResidualSample>& rows) const
{
    if (rows.empty()) return;
    static bool firstCall = true;
    const std::string fname = cfg_.outputDir + "/InnerResidualHistory.dat";
    std::ofstream out(fname, firstCall ? std::ios::trunc : std::ios::app);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << '\n';
        return;
    }
    if (firstCall) {
        out << "TITLE = \"DT pseudo-time residual history\"\n"
            << "VARIABLES = \"step\", \"physicalStep\", \"inner\", \"runTime\", \"runCyc\", "
            << "\"pseudoCFL\", \"dtPhys\", \"dtau\", \"Gfluid\", \"GNe\", \"GNi\", "
            << "\"GEe\", \"GPhi\", \"logRelFluid\", \"logRelNe\", \"logRelNi\", "
            << "\"logRelEe\", \"logRelPhi\"\n";
        firstCall = false;
    }

    auto logRel = [](scalar gm, scalar gInitial) {
        const scalar denom = std::max(gInitial, 1.0e-300);
        return std::log10(std::max(gm / denom, 1.0e-300));
    };

    out << "ZONE T=\"step_" << step_ << "_phys_" << physicalStep_
        << "\", I=" << rows.size()
        << ", SOLUTIONTIME=" << std::scientific << std::setprecision(8) << runCyc_
        << ", F=POINT\n";

    for (const InnerResidualSample& row : rows) {
        const auto& g = row.g;
        const auto& g0 = row.g0;
        out << step_ << '\t'
            << physicalStep_ << '\t'
            << row.inner << '\t'
            << std::scientific << std::setprecision(8)
            << runTime_ * cfg_.tRef << '\t'
            << runCyc_ << '\t'
            << cfg_.pseudoCFL << '\t'
            << dtPhysStep_ * cfg_.tRef << '\t'
            << row.dtau * cfg_.tRef << '\t'
            << g[0] << '\t' << g[1] << '\t' << g[2] << '\t' << g[3] << '\t' << g[4] << '\t'
            << logRel(g[0], g0[0]) << '\t'
            << logRel(g[1], g0[1]) << '\t'
            << logRel(g[2], g0[2]) << '\t'
            << logRel(g[3], g0[3]) << '\t'
            << logRel(g[4], g0[4]) << '\n';
    }
}

scalar Solver::cellVal(const Cell& c, label t, label v) const
{
    if (v == 0) return c.Ne(t);
    if (v == 1) return c.Ni(t);
    return c.Ee(t);
}

scalar Solver::minVar(label v) const
{
    scalar out = std::numeric_limits<scalar>::max();
    for (const Cell& c : cells_) {
        for (label t = 0; t < cfg_.numT; ++t) out = std::min(out, cellVal(c, t, v));
    }
    return out;
}

scalar Solver::maxVar(label v) const
{
    scalar out = -std::numeric_limits<scalar>::max();
    for (const Cell& c : cells_) {
        for (label t = 0; t < cfg_.numT; ++t) out = std::max(out, cellVal(c, t, v));
    }
    return out;
}

bool Solver::allFinite() const
{
    for (const Cell& c : cells_) {
        for (label t = 0; t < cfg_.numT; ++t) {
            if (!std::isfinite(c.Ne(t)) || !std::isfinite(c.Ni(t)) ||
                !std::isfinite(c.Ee(t)) || !std::isfinite(c.Phi(t)) ||
                c.Ne(t) <= 0.0 || c.Ni(t) <= 0.0 || c.Ee(t) <= 0.0) {
                return false;
            }
        }
    }
    return true;
}

void Solver::writeConnectivity(std::ostream& out) const
{
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            out << nid(i, j) + 1 << " "
                << nid(i + 1, j) + 1 << " "
                << nid(i + 1, j + 1) + 1 << " "
                << nid(i, j + 1) + 1 << "\n";
        }
    }
}

void Solver::writeMeshGeometry() const
{
    const std::string fname = cfg_.outputDir + "/mesh_geom.dat";
    std::ofstream out(fname, std::ios::trunc);
    if (!out) {
        std::cerr << "[mesh] cannot write " << fname << "\n";
        return;
    }
    out << std::scientific << std::setprecision(8);
    out << "TITLE = \"Mesh geometry (per cell)\"\n"
        << "VARIABLES = \"X[m]\" \"Y[m]\" \"blockID\" \"cellID\" \"volume\"\n"
        << "ZONE T=\"block_0\" NODES=" << nNodes_
        << " ELEMENTS=" << nCells_
        << " DATAPACKING=BLOCK ZONETYPE=FEQUADRILATERAL"
        << " VARLOCATION=([3-5]=CELLCENTERED)\n";
    writeNodeBlock(out, nNodes_, [&](label k) { return nodes_[k].x * cfg_.LRef; });
    writeNodeBlock(out, nNodes_, [&](label k) { return nodes_[k].y * cfg_.LRef; });
    writeCellBlock(out, [&](const Cell&, label) { return 0.0; });
    writeCellBlock(out, [&](const Cell&, label cID) { return scalar(cID); });
    writeCellBlock(out, [&](const Cell& c, label) { return c.vol * cfg_.LRef * cfg_.LRef; });
    writeConnectivity(out);
    std::cout << "[mesh] wrote " << fname << "\n";
}

void Solver::writeFieldTecplotMergedAllIT(const std::string& fname,
                                          std::vector<label> iTList,
                                          bool append,
                                          const std::string& zoneTag) const
{
    if (iTList.empty()) {
        iTList.resize(cfg_.numT);
        for (label k = 0; k < cfg_.numT; ++k) iTList[k] = k;
    }
    std::ofstream out(fname, append ? std::ios::app : std::ios::trunc);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << "\n";
        return;
    }
    out << std::scientific << std::setprecision(8);
    if (!append) {
        out << "TITLE = \"plasma2d merged field (SI)\"\n"
            << "VARIABLES = \"X[m]\" \"Y[m]\" \"blockID\" "
               "\"Ne[1/m^3]\" \"Ni[1/m^3]\" \"Ee[J/m^3]\" \"Te[K]\" "
               "\"Phi[V]\" \"Ec_x[V/m]\" \"Ec_y[V/m]\"\n";
    }
    const label nNode = nNodes_;
    const scalar eRefSI = cfg_.EeRef;
    const scalar ErefSI = cfg_.phiRef / cfg_.LRef;
    for (label z = 0; z < static_cast<label>(iTList.size()); ++z) {
        const label t = iTList[z];
        out << "ZONE T=\"";
        if (!zoneTag.empty()) out << zoneTag << "_";
        out << "iT=" << t << "\""
            << " NODES=" << nNode << " ELEMENTS=" << nCells_
            << " DATAPACKING=BLOCK ZONETYPE=FEQUADRILATERAL"
            << " VARLOCATION=([4-10]=CELLCENTERED)\n";
        writeNodeBlock(out, nNode, [&](label k) { return nodes_[k].x * cfg_.LRef; });
        writeNodeBlock(out, nNode, [&](label k) { return nodes_[k].y * cfg_.LRef; });
        writeNodeBlock(out, nNode, [&](label) { return 0.0; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Ne(t) * cfg_.nRef; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Ni(t) * cfg_.nRef; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Ee(t) * eRefSI; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Te(t) * cfg_.TeRef; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Phi(t) * cfg_.phiRef; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Ex(t) * ErefSI; });
        writeCellBlock(out, [&](const Cell& c, label) { return c.Ey(t) * ErefSI; });
        writeConnectivity(out);
    }
}

void Solver::writeIterSolution(bool force) const
{
    if (!force && (cfg_.writeStep <= 0 || step_ % cfg_.writeStep != 0)) return;
    static bool firstCall = true;
    const bool append = !firstCall;
    firstCall = false;
    writeFieldTecplotMergedAllIT(cfg_.outputDir + "/iterSolution.dat", {},
                                 append, "Step_" + std::to_string(step_));
}

void Solver::writeFinalSolution(bool force) const
{
    if (!force && (cfg_.writeStep <= 0 || step_ % cfg_.writeStep != 0)) return;
    writeFieldTecplotMergedAllIT(cfg_.outputDir + "/finalSolution.dat", {}, false, "");
}

void Solver::writeRestart(bool force) const
{
    if (!force && (cfg_.writeStep <= 0 || step_ % cfg_.writeStep != 0)) return;
    const std::string fname = cfg_.outputDir + "/restart.dat";
    std::ofstream out(fname, std::ios::trunc);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << '\n';
        return;
    }
    out << "TITLE = \"Restart Solution\"\n"
        << "VARIABLES = \"cellID\", \"x\", \"y\", \"Ne\", \"Ni\", \"Ee\", \"Te\", \"Phi\"\n";
    for (label t = 0; t < cfg_.numT; ++t) {
        out << "ZONE T=\"Init_nT=" << t << "\", I=" << cfg_.nCellX
            << ", J=" << cfg_.nCellY << ", F=POINT\n";
        for (label j = 0; j < cfg_.nCellY; ++j) {
            for (label i = 0; i < cfg_.nCellX; ++i) {
                const label cID = cid(i, j);
                const Cell& c = cells_[cID];
                out << cID << '\t' << std::setprecision(10) << c.x << '\t' << c.y << '\t'
                    << c.Ne(t) << '\t' << c.Ni(t) << '\t' << c.Ee(t) << '\t'
                    << c.Te(t) << '\t' << c.Phi(t) << '\n';
            }
        }
    }
}


void Solver::writeAverageNeNiEeDT(scalar wall, label targetCycle, scalar theta) const
{
    static bool firstCall = true;
    const std::string fname = cfg_.outputDir + "/NeNiEeAveDT.dat";
    std::ofstream out(fname, firstCall ? std::ios::trunc : std::ios::app);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << '\n';
        return;
    }
    if (firstCall) {
        out << "TITLE = \"DT cycle-start spatial average\"\n"
            << "VARIABLES = \"step\", \"wallTime\", \"targetCycle\", "
               "\"runTime\", \"runCyc\", \"NeAve\", \"NiAve\", \"EeAve\"\n"
            << "ZONE T=\"Average_phase0\", F=POINT\n";
        firstCall = false;
    }
    scalar Ne = 0.0;
    scalar Ni = 0.0;
    scalar Ee = 0.0;
    scalar vol = 0.0;
    const scalar th = std::max(scalar(0.0), std::min(scalar(1.0), theta));
    auto lerp = [&](scalar a, scalar b) { return a + th * (b - a); };
    for (const Cell& c : cells_) {
        Ne += lerp(c.NePhysOld(0), c.Ne(0)) * c.vol;
        Ni += lerp(c.NiPhysOld(0), c.Ni(0)) * c.vol;
        Ee += lerp(c.EePhysOld(0), c.Ee(0)) * c.vol;
        vol += c.vol;
    }
    const scalar invVol = 1.0 / std::max(vol, tiny);
    const scalar targetTime = scalar(targetCycle) * cfg_.period;
    out << step_ << '\t'
        << std::fixed << std::setprecision(8) << wall << '\t'
        << targetCycle << '\t'
        << std::scientific << std::setprecision(8)
        << targetTime * cfg_.tRef << '\t'
        << std::fixed << std::setprecision(8) << scalar(targetCycle) << '\t'
        << std::scientific << std::setprecision(8)
        << Ne * invVol * cfg_.nRef << '\t'
        << Ni * invVol * cfg_.nRef << '\t'
        << Ee * invVol * cfg_.EeRef << '\n';
}


void Solver::writeUnsteadyFlowFieldDT(label sampleID, scalar phase, scalar targetCycle,
                                      scalar theta, bool append) const
{
    const std::string fname = cfg_.outputDir + "/UnsteadyFlowFieldDT.dat";
    std::ofstream out(fname, append ? std::ios::app : std::ios::trunc);
    if (!out) {
        std::cerr << "[Solver] cannot write " << fname << '\n';
        return;
    }
    out << std::scientific << std::setprecision(8);
    if (!append) {
        out << "TITLE = \"DT Interpolated Flow Field\"\n"
            << "VARIABLES = \"X_m\", \"Y_m\", \"blockID\", \"cellID\", "
               "\"Ne\", \"Ni\", \"Ee\", \"Te_K\", \"Te_eV\", \"Phi\"\n";
    }
    const label nNode = nNodes_;
    const scalar kToEv = cfg_.kB / cfg_.eV2J;
    const scalar th = std::max(scalar(0.0), std::min(scalar(1.0), theta));
    auto lerp = [&](scalar a, scalar b) { return a + th * (b - a); };
    out << "ZONE T=\"Step_" << step_
        << "_sample=" << sampleID
        << "_phase=" << phase
        << "_cycle=" << targetCycle
        << "\", NODES=" << nNode << " ELEMENTS=" << nCells_
        << " DATAPACKING=BLOCK ZONETYPE=FEQUADRILATERAL"
        << " VARLOCATION=([4-10]=CELLCENTERED)"
        << " STRANDID=2 SOLUTIONTIME=" << phase << "\n";
    writeNodeBlock(out, nNode, [&](label k) { return nodes_[k].x * cfg_.LRef; });
    writeNodeBlock(out, nNode, [&](label k) { return nodes_[k].y * cfg_.LRef; });
    writeNodeBlock(out, nNode, [&](label) { return 0.0; });
    writeCellBlock(out, [&](const Cell&, label cID) { return scalar(cID); });
    writeCellBlock(out, [&](const Cell& c, label) {
        return lerp(c.NePhysOld(0), c.Ne(0)) * cfg_.nRef;
    });
    writeCellBlock(out, [&](const Cell& c, label) {
        return lerp(c.NiPhysOld(0), c.Ni(0)) * cfg_.nRef;
    });
    writeCellBlock(out, [&](const Cell& c, label) {
        return lerp(c.EePhysOld(0), c.Ee(0)) * cfg_.EeRef;
    });
    writeCellBlock(out, [&](const Cell& c, label) {
        const scalar Ne = lerp(c.NePhysOld(0), c.Ne(0));
        const scalar Ee = lerp(c.EePhysOld(0), c.Ee(0));
        return Ee / std::max(Ne, cfg_.qFloor) * cfg_.TeRef;
    });
    writeCellBlock(out, [&](const Cell& c, label) {
        const scalar Ne = lerp(c.NePhysOld(0), c.Ne(0));
        const scalar Ee = lerp(c.EePhysOld(0), c.Ee(0));
        return Ee / std::max(Ne, cfg_.qFloor) * cfg_.TeRef * kToEv;
    });
    writeCellBlock(out, [&](const Cell& c, label) {
        return lerp(c.PhiPhysOld(0), c.Phi(0)) * cfg_.phiRef;
    });
    writeConnectivity(out);
}





void Solver::writeDueDTAverages(scalar cycleBefore, scalar wall)
{
    constexpr scalar tol = 1.0e-12;
    const scalar denom = runCyc_ - cycleBefore;
    while (scalar(nextAverageCycle_) <= runCyc_ + tol) {
        if (scalar(nextAverageCycle_) > cycleBefore + tol) {
            const scalar theta = std::abs(denom) > tiny
                               ? (scalar(nextAverageCycle_) - cycleBefore) / denom
                               : 1.0;
            writeAverageNeNiEeDT(wall, nextAverageCycle_, theta);
        }
        ++nextAverageCycle_;
    }
}

void Solver::writeDueDTFields(scalar cycleBefore)
{
    if (dtSamplePhases_.empty()) return;
    constexpr scalar tol = 1.0e-12;
    const scalar lastCycleStart = std::max(scalar(0.0), cfg_.stopPeriod - scalar(1.0));
    const scalar denom = runCyc_ - cycleBefore;
    while (nextDTSample_ < static_cast<label>(dtSamplePhases_.size())) {
        const scalar phase = dtSamplePhases_[static_cast<size_t>(nextDTSample_)];
        const scalar targetCycle = lastCycleStart + phase;
        if (targetCycle < cycleBefore - tol) {
            ++nextDTSample_;
            continue;
        }
        if (targetCycle <= runCyc_ + tol) {
            const scalar theta = std::abs(denom) > tiny
                               ? (targetCycle - cycleBefore) / denom
                               : 1.0;
            writeUnsteadyFlowFieldDT(nextDTSample_, phase, targetCycle, theta, dtFieldWritten_);
            dtFieldWritten_ = true;
            ++nextDTSample_;
            continue;
        }
        break;
    }
}

}
