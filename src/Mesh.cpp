#include "ccp-2d.hpp"

namespace ccp2d {

scalar dist2D(const Node& a, const Node& b)
{
    const scalar dx = a.x - b.x;
    const scalar dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

scalar quadArea(const Node& a, const Node& b, const Node& c, const Node& d)
{
    const scalar s =
        a.x * b.y - b.x * a.y +
        b.x * c.y - c.x * b.y +
        c.x * d.y - d.x * c.y +
        d.x * a.y - a.x * d.y;
    return 0.5 * std::abs(s);
}

bool Solver::readPlot3DMesh(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) return false;
    std::vector<std::string> tk;
    std::string line, word;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream iss(line);
        while (iss >> word) tk.push_back(word);
    }
    if (tk.empty()) throw std::runtime_error("[Plot3D] empty mesh file: " + filename);

    auto toLabel = [](const std::string& s) { return static_cast<label>(std::stoll(s)); };
    auto toScalar = [](const std::string& s) { return std::stod(s); };
    const label nblocks = toLabel(tk[0]);
    if (nblocks != 1) {
        throw std::runtime_error("[Plot3D] src_2d_ad currently supports one structured block");
    }

    auto tryLayout = [&](label dimsPer, label coordsPer, label& ni, label& nj) {
        if (tk.size() < static_cast<size_t>(1 + dimsPer)) return false;
        try {
            ni = toLabel(tk[1]);
            nj = toLabel(tk[2]);
            label nk = 1;
            if (dimsPer == 3) nk = toLabel(tk[3]);
            if (ni < 2 || nj < 2 || nk != 1) return false;
            const size_t p0 = 1 + static_cast<size_t>(dimsPer);
            const size_t need = p0 + static_cast<size_t>(coordsPer) * ni * nj;
            if (need > tk.size()) return false;
            for (size_t p = p0; p < need; ++p) (void)toScalar(tk[p]);
            return true;
        } catch (...) {
            return false;
        }
    };

    label ni = 0;
    label nj = 0;
    label dimsPer = 0;
    label coordsPer = 0;
    if (tryLayout(3, 3, ni, nj)) {
        dimsPer = 3;
        coordsPer = 3;
    } else if (tryLayout(3, 2, ni, nj)) {
        dimsPer = 3;
        coordsPer = 2;
    } else if (tryLayout(2, 2, ni, nj)) {
        dimsPer = 2;
        coordsPer = 2;
    } else {
        throw std::runtime_error("[Plot3D] could not auto-detect mesh layout: " + filename);
    }

    cfg_.nNodeX = ni;
    cfg_.nNodeY = nj;
    cfg_.nCellX = cfg_.nNodeX - 1;
    cfg_.nCellY = cfg_.nNodeY - 1;
    nCells_ = cfg_.nCellX * cfg_.nCellY;
    nNodes_ = cfg_.nNodeX * cfg_.nNodeY;
    nodes_.assign(static_cast<size_t>(nNodes_), Node{});
    size_t p = 1 + static_cast<size_t>(dimsPer);
    for (label j = 0; j < nj; ++j) {
        for (label i = 0; i < ni; ++i) nodes_[nid(i, j)].x = toScalar(tk[p++]) / cfg_.LRef;
    }
    for (label j = 0; j < nj; ++j) {
        for (label i = 0; i < ni; ++i) nodes_[nid(i, j)].y = toScalar(tk[p++]) / cfg_.LRef;
    }
    if (coordsPer == 3) p += static_cast<size_t>(ni * nj);

    std::cout << "[Plot3D] loaded 1 block(s) from " << filename
              << "  (form " << (dimsPer == 3 ? "2D-as-3D" : "pure-2D")
              << (coordsPer == 3 ? " with z" : "") << ")\n"
              << "         block 0:  " << ni << " x " << nj
              << " nodes  (" << cfg_.nCellX << " x " << cfg_.nCellY << " cells)\n"
              << "         total: " << ni * nj << " nodes\n";
    return true;
}

Node Solver::cellCenter(label i, label j) const
{
    const Node& n0 = nodes_[nid(i, j)];
    const Node& n1 = nodes_[nid(i + 1, j)];
    const Node& n2 = nodes_[nid(i + 1, j + 1)];
    const Node& n3 = nodes_[nid(i, j + 1)];
    return {(n0.x + n1.x + n2.x + n3.x) * 0.25,
            (n0.y + n1.y + n2.y + n3.y) * 0.25};
}

Node Solver::faceCenterX(label i, label j) const
{
    const Node& a = nodes_[nid(i, j)];
    const Node& b = nodes_[nid(i, j + 1)];
    return {(a.x + b.x) * 0.5, (a.y + b.y) * 0.5};
}

Node Solver::faceCenterY(label i, label j) const
{
    const Node& a = nodes_[nid(i, j)];
    const Node& b = nodes_[nid(i + 1, j)];
    return {(a.x + b.x) * 0.5, (a.y + b.y) * 0.5};
}

void Solver::applyBoundaryFile()
{
    if (cfg_.bccFile.empty()) {
        throw std::runtime_error("bccFile is required; 2-D boundary conditions must come from the BCC file");
    }
    std::ifstream in(cfg_.bccFile);
    if (!in) {
        throw std::runtime_error("cannot open BCC file: " + cfg_.bccFile);
    }

    auto wallTeFromEv = [&](scalar teEv) {
        return teEv * cfg_.eV2J / cfg_.kB / cfg_.TeRef;
    };
    auto wallTeToEv = [&](scalar teNondim) {
        return teNondim * cfg_.TeRef * cfg_.kB / cfg_.eV2J;
    };
    auto wallModelName = [](WallModel model) {
        if (model == WallModel::Lymberopoulos) return std::string("LYWALL");
        if (model == WallModel::Sakiyama) return std::string("SAWALL");
        return std::string("NONE");
    };
    auto parseWallModel = [](const std::string& raw) {
        const std::string value = upper(raw);
        if (value == "LYWALL" || value == "LY" ||
            value == "LYMBEROPOULOS" || value == "LYMBEROPOULOSBC") {
            return WallModel::Lymberopoulos;
        }
        if (value == "SAWALL" || value == "SA" ||
            value == "SAKIYAMA" || value == "SAKIYAMABC") {
            return WallModel::Sakiyama;
        }
        throw std::runtime_error("[BCC] invalid wallBC: " + raw);
    };

    auto setPatch = [&](Face& f, FaceKind dirichletKind,
                        const std::string& rawType,
                        scalar phiAmpRaw, bool hasAmp,
                        scalar phaseRaw, bool hasPhase,
                        WallModel wallModel,
                        scalar wallTeRaw, bool hasWallTe) {
        const std::string type = upper(rawType);
        f.bcSet = true;
        f.phiBC.setZero();
        f.timeDependentPhi = false;
        f.phiAmp = 0.0;
        f.phiPhase = 0.0;
        f.wallModel = WallModel::None;
        f.wallTe = 0.0;
        if (type == "SYMMETRY" || type == "AXIS") {
            f.kind = FaceKind::Symmetry;
            return;
        }
        f.wallModel = wallModel;
        f.wallTe = hasWallTe ? wallTeFromEv(wallTeRaw) : 0.0;
        if (type == "PHI_GROUND" || type == "GROUND") {
            f.kind = dirichletKind;
            return;
        }
        if (type == "PHI_RF") {
            f.kind = dirichletKind;
            const scalar amp = hasAmp ? phiAmpRaw / cfg_.phiRef : cfg_.PhiRF;
            const scalar phase = hasPhase ? phaseRaw * pi / 180.0 : cfg_.phase;
            f.timeDependentPhi = true;
            f.phiAmp = amp;
            f.phiPhase = phase;
            for (label t = 0; t < cfg_.numT; ++t) {
                f.phiBC(t) = amp * std::sin(phase);
            }
            return;
        }
        if (type == "PHI_AR") {
            f.kind = dirichletKind;
            f.phiBC.setConstant(cfg_.PhiAR);
            return;
        }

        throw std::runtime_error("[BCC] unsupported patch type: " + rawType);
    };

    bool inPatches = false;
    label patchRecords = 0;
    label assignedFaces = 0;
    std::vector<std::string> patchSummaries;
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            inPatches = upper(line) == "[PATCHES]";
            continue;
        }
        if (!inPatches) continue;

        std::istringstream iss(line);
        std::vector<std::string> tk;
        std::string s;
        while (iss >> s) tk.push_back(s);
        if (tk.size() < 6) continue;
        ++patchRecords;

        const std::string side = upper(tk[2]);
        const label i0 = std::max<label>(0, static_cast<label>(std::stoll(tk[3])));
        const label i1 = std::max<label>(i0, static_cast<label>(std::stoll(tk[4])));
        const std::string type = tk[5];
        scalar amp = 0.0;
        scalar phase = 0.0;
        scalar wallTe = 0.0;
        WallModel wallModel = WallModel::None;
        bool hasAmp = false;
        bool hasPhase = false;
        bool hasWallTe = false;
        bool hasWallModel = false;
        for (size_t k = 6; k < tk.size(); ++k) {
            const auto eq = tk[k].find('=');
            if (eq == std::string::npos) continue;
            const std::string key = upper(tk[k].substr(0, eq));
            const std::string rawValue = tk[k].substr(eq + 1);
            if (key == "PHIAMP" || key == "AMP") {
                const scalar value = std::stod(rawValue);
                amp = value;
                hasAmp = true;
            } else if (key == "PHASE") {
                const scalar value = std::stod(rawValue);
                phase = value;
                hasPhase = true;
            } else if (key == "WALLTE" || key == "TEWALL" ||
                       key == "ELECTRONTE" || key == "TEE" || key == "TE") {
                const scalar value = std::stod(rawValue);
                wallTe = value;
                hasWallTe = true;
            } else if (key == "WALLBC" || key == "WALLMODEL" || key == "FLUXBC") {
                wallModel = parseWallModel(rawValue);
                hasWallModel = true;
            }
        }

        const std::string patchName = tk[0];
        const std::string typeUpper = upper(type);
        const bool isSymmetry = typeUpper == "SYMMETRY" || typeUpper == "AXIS";
        if (isSymmetry && (hasWallModel || hasWallTe)) {
            throw std::runtime_error("[BCC] patch " + patchName +
                                     " is symmetry/axis and must not specify wallBC or wallTe");
        }
        if (!isSymmetry && !hasWallModel) {
            throw std::runtime_error("[BCC] patch " + patchName +
                                     " must specify wallBC=LYWALL or wallBC=SAWALL");
        }
        if (wallModel == WallModel::Lymberopoulos) {
            if (!hasWallTe) {
                throw std::runtime_error("[BCC] patch " + patchName +
                                         " uses wallBC=LYWALL and must specify wallTe=<eV>");
            }
            if (cfg_.Ks <= 0.0) {
                throw std::runtime_error("[BCC] patch " + patchName +
                                         " uses wallBC=LYWALL and requires positive Ks in the case file");
            }
        } else if (wallModel == WallModel::Sakiyama && hasWallTe) {
            throw std::runtime_error("[BCC] patch " + patchName +
                                     " uses wallBC=SAWALL and must not specify wallTe");
        }

        const label assignedBefore = assignedFaces;
        if (side == "IMIN") {
            for (label j = i0; j < std::min(i1, cfg_.nCellY); ++j) {
                setPatch(xFaces_[xf(0, j)], FaceKind::MinusBC, type,
                         amp, hasAmp, phase, hasPhase, wallModel, wallTe, hasWallTe);
                ++assignedFaces;
            }
        } else if (side == "IMAX") {
            for (label j = i0; j < std::min(i1, cfg_.nCellY); ++j) {
                setPatch(xFaces_[xf(cfg_.nCellX, j)], FaceKind::PlusBC, type,
                         amp, hasAmp, phase, hasPhase, wallModel, wallTe, hasWallTe);
                ++assignedFaces;
            }
        } else if (side == "JMIN") {
            for (label i = i0; i < std::min(i1, cfg_.nCellX); ++i) {
                setPatch(yFaces_[yf(i, 0)], FaceKind::MinusBC, type,
                         amp, hasAmp, phase, hasPhase, wallModel, wallTe, hasWallTe);
                ++assignedFaces;
            }
        } else if (side == "JMAX") {
            for (label i = i0; i < std::min(i1, cfg_.nCellX); ++i) {
                setPatch(yFaces_[yf(i, cfg_.nCellY)], FaceKind::PlusBC, type,
                         amp, hasAmp, phase, hasPhase, wallModel, wallTe, hasWallTe);
                ++assignedFaces;
            }
        } else {
            throw std::runtime_error("[BCC] unsupported patch side: " + side);
        }

        const label nAssignedPatch = assignedFaces - assignedBefore;
        std::ostringstream os;
        os << std::scientific << std::setprecision(6)
           << "      " << tk[0] << ": side=" << side
           << " type=" << typeUpper
           << " faces=" << nAssignedPatch;
        if (typeUpper == "SYMMETRY" || typeUpper == "AXIS") {
            os << "  dPhi/dn=0, normal flux=0";
        } else if (typeUpper == "PHI_GROUND" || typeUpper == "GROUND") {
            os << "  phi=0.000000e+00 V, phi*=0.000000e+00";
        } else if (typeUpper == "PHI_AR") {
            os << "  phi=" << cfg_.PhiAR * cfg_.phiRef << " V"
               << ", phi*=" << cfg_.PhiAR;
        } else if (typeUpper == "PHI_RF") {
            const scalar ampDim = hasAmp ? amp : cfg_.PhiRF * cfg_.phiRef;
            const scalar ampND = ampDim / cfg_.phiRef;
            const scalar phaseDeg = hasPhase ? phase : cfg_.phase * 180.0 / pi;
            const scalar phaseRad = phaseDeg * pi / 180.0;
            os << "  amp=" << ampDim << " V"
               << ", amp*=" << ampND
               << ", phase=" << phaseDeg << " deg"
               << " (" << phaseRad << " rad)";
        }
        if (!isSymmetry) {
            os << ", wallBC=" << wallModelName(wallModel);
            if (wallModel == WallModel::Lymberopoulos) {
                const scalar teND = wallTeFromEv(wallTe);
                os << ", wallTe=" << wallTeToEv(teND) << " eV"
                   << ", wallTe*=" << teND;
            }
        }
        patchSummaries.push_back(os.str());
    }

    std::cout << "[BCC] loaded boundary conditions from " << cfg_.bccFile << "\n"
              << "      patch records: " << patchRecords
              << ", boundary faces assigned: " << assignedFaces << "\n";
    for (const std::string& summary : patchSummaries) {
        std::cout << summary << "\n";
    }

    for (label j = 0; j < cfg_.nCellY; ++j) {
        if (!xFaces_[xf(0, j)].bcSet) {
            throw std::runtime_error("[BCC] IMIN boundary is not fully specified");
        }
        if (!xFaces_[xf(cfg_.nCellX, j)].bcSet) {
            throw std::runtime_error("[BCC] IMAX boundary is not fully specified");
        }
    }
    for (label i = 0; i < cfg_.nCellX; ++i) {
        if (!yFaces_[yf(i, 0)].bcSet) {
            throw std::runtime_error("[BCC] JMIN boundary is not fully specified");
        }
        if (!yFaces_[yf(i, cfg_.nCellY)].bcSet) {
            throw std::runtime_error("[BCC] JMAX boundary is not fully specified");
        }
    }
}

void Solver::initGridAndFaces()
{
    makeDir(cfg_.outputDir);
    if (cfg_.meshFile.empty()) {
        throw std::runtime_error("meshFile is required");
    }
    if (!readPlot3DMesh(cfg_.meshFile)) {
        throw std::runtime_error("cannot open mesh file: " + cfg_.meshFile);
    }

    cells_.assign(static_cast<size_t>(nCells_), Cell(cfg_.numT));
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            Cell& c = cells_[cid(i, j)];
            const Node cc = cellCenter(i, j);
            c.x = cc.x;
            c.y = cc.y;
            c.vol = quadArea(nodes_[nid(i, j)], nodes_[nid(i + 1, j)],
                             nodes_[nid(i + 1, j + 1)], nodes_[nid(i, j + 1)]);
            c.Ne.setConstant(cfg_.Ne0);
            c.Ni.setConstant(cfg_.Ni0);
            c.Te.setConstant(cfg_.Te0);
            c.Ee.setConstant(std::max(cfg_.Ne0 * cfg_.Te0, cfg_.qFloor));
            c.Phi.setConstant(cfg_.Phi0);
            c.NeOld = c.Ne;
            c.NiOld = c.Ni;
            c.EeOld = c.Ee;
            c.PhiOld = c.Phi;
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

    xFaces_.assign(static_cast<size_t>((cfg_.nCellX + 1) * cfg_.nCellY), Face(cfg_.numT));
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i <= cfg_.nCellX; ++i) {
            Face& f = xFaces_[xf(i, j)];
            f.area = dist2D(nodes_[nid(i, j)], nodes_[nid(i, j + 1)]);
            const Node fc = faceCenterX(i, j);
            if (i == 0) {
                f.cR = cid(0, j);
                f.dist = dist2D(fc, {cells_[f.cR].x, cells_[f.cR].y});
            } else if (i == cfg_.nCellX) {
                f.cL = cid(cfg_.nCellX - 1, j);
                f.dist = dist2D(fc, {cells_[f.cL].x, cells_[f.cL].y});
            } else {
                f.kind = FaceKind::Interior;
                f.cL = cid(i - 1, j);
                f.cR = cid(i, j);
                f.dist = dist2D({cells_[f.cL].x, cells_[f.cL].y},
                                {cells_[f.cR].x, cells_[f.cR].y});
            }
            f.dist = std::max(f.dist, 1.0e-30);
        }
    }

    yFaces_.assign(static_cast<size_t>(cfg_.nCellX * (cfg_.nCellY + 1)), Face(cfg_.numT));
    for (label j = 0; j <= cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            Face& f = yFaces_[yf(i, j)];
            f.area = dist2D(nodes_[nid(i, j)], nodes_[nid(i + 1, j)]);
            const Node fc = faceCenterY(i, j);
            if (j == 0 || j == cfg_.nCellY) {
                if (j > 0) f.cL = cid(i, j - 1);
                if (j < cfg_.nCellY) f.cR = cid(i, j);
                const label ref = j == 0 ? cid(i, 0) : cid(i, cfg_.nCellY - 1);
                f.dist = dist2D(fc, {cells_[ref].x, cells_[ref].y});
            } else {
                f.kind = FaceKind::Interior;
                f.cL = cid(i, j - 1);
                f.cR = cid(i, j);
                f.dist = dist2D({cells_[f.cL].x, cells_[f.cL].y},
                                {cells_[f.cR].x, cells_[f.cR].y});
            }
            f.dist = std::max(f.dist, 1.0e-30);
        }
    }

    applyBoundaryFile();
    writeMeshGeometry();
}

}
