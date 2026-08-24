#include "ccp-2d.hpp"

namespace ccp2d {

scalar Solver::semiImplicitPoissonFaceCoeff(label c, label nb, label t) const
{
    return cfg_.e * cfg_.muE * 0.5 * (cells_[c].Ne(t) + cells_[nb].Ne(t)) * poissonTimeStep(cells_[c]);
}

scalar Solver::semiImplicitPoissonBoundaryCoeff(label c, label t) const
{
    return cfg_.e * cfg_.muE * cells_[c].Ne(t) * poissonTimeStep(cells_[c]);
}

scalar Solver::poissonTimeStep(const Cell& c) const
{
    return c.dt;
}

bool Solver::hasPoissonDirichletBoundary() const
{
    auto hasDirichlet = [](const Face& f) {
        return f.kind == FaceKind::MinusBC || f.kind == FaceKind::PlusBC;
    };
    for (const Face& f : xFaces_) {
        if (hasDirichlet(f)) return true;
    }
    for (const Face& f : yFaces_) {
        if (hasDirichlet(f)) return true;
    }
    return false;
}

void Solver::assemblePhysicalPoisson(label t, std::vector<Eigen::Triplet<scalar>>& triplets, Vec& b) const
{



    for (label cID = 0; cID < nCells_; ++cID) {
        const label i = ci(cID);
        const label j = cj(cID);
        const Face& fL = xFaces_[xf(i, j)];
        const Face& fR = xFaces_[xf(i + 1, j)];
        const Face& fB = yFaces_[yf(i, j)];
        const Face& fT = yFaces_[yf(i, j + 1)];
        scalar diag = 0.0;

        auto addInterior = [&](label nb, const Face& f) {
            if (f.kind == FaceKind::Symmetry) return;
            const scalar a = cfg_.eps0 * f.area / f.dist;
            diag += a;
            triplets.emplace_back(cID, nb, -a);
        };

        auto addBoundary = [&](const Face& f) {
            if (f.kind == FaceKind::Symmetry) return;
            const scalar a = cfg_.eps0 * f.area / f.dist;
            diag += a;
            b(cID) += a * f.phiBC(t);
        };




        if (i > 0) addInterior(cid(i - 1, j), fL);
        else addBoundary(fL);
        if (i + 1 < cfg_.nCellX) addInterior(cid(i + 1, j), fR);
        else addBoundary(fR);
        if (j > 0) addInterior(cid(i, j - 1), fB);
        else addBoundary(fB);
        if (j + 1 < cfg_.nCellY) addInterior(cid(i, j + 1), fT);
        else addBoundary(fT);

        const Cell& cc = cells_[cID];
        b(cID) += cfg_.e * (cc.Ni(t) - cc.Ne(t)) * cc.vol;
        triplets.emplace_back(cID, cID, diag);
    }
}

void Solver::addSemiImplicitPoissonCorrection(label t, std::vector<Eigen::Triplet<scalar>>& triplets, Vec& b) const
{




    for (label cID = 0; cID < nCells_; ++cID) {
        const label i = ci(cID);
        const label j = cj(cID);
        const Face& fL = xFaces_[xf(i, j)];
        const Face& fR = xFaces_[xf(i + 1, j)];
        const Face& fB = yFaces_[yf(i, j)];
        const Face& fT = yFaces_[yf(i, j + 1)];
        const Cell& cc = cells_[cID];
        scalar diag = 0.0;

        auto addInterior = [&](label nb, const Face& f) {
            if (f.kind == FaceKind::Symmetry) return;
            const scalar a = semiImplicitPoissonFaceCoeff(cID, nb, t) * f.area / f.dist;
            diag += a;
            triplets.emplace_back(cID, nb, -a);
        };

        auto addBoundary = [&](const Face& f) {
            if (f.kind == FaceKind::Symmetry) return;
            const scalar a = semiImplicitPoissonBoundaryCoeff(cID, t) * f.area / f.dist;
            diag += a;
            b(cID) += a * f.phiBC(t);
        };

        if (i > 0) addInterior(cid(i - 1, j), fL);
        else addBoundary(fL);
        if (i + 1 < cfg_.nCellX) addInterior(cid(i + 1, j), fR);
        else addBoundary(fR);
        if (j > 0) addInterior(cid(i, j - 1), fB);
        else addBoundary(fB);
        if (j + 1 < cfg_.nCellY) addInterior(cid(i, j + 1), fT);
        else addBoundary(fT);

        const scalar gradXL = i > 0 ? (cc.Ne(t) - cells_[cid(i - 1, j)].Ne(t)) / fL.dist : 0.0;
        const scalar gradXR = i + 1 < cfg_.nCellX ? (cells_[cid(i + 1, j)].Ne(t) - cc.Ne(t)) / fR.dist : 0.0;
        const scalar gradYB = j > 0 ? (cc.Ne(t) - cells_[cid(i, j - 1)].Ne(t)) / fB.dist : 0.0;
        const scalar gradYT = j + 1 < cfg_.nCellY ? (cells_[cid(i, j + 1)].Ne(t) - cc.Ne(t)) / fT.dist : 0.0;
        const scalar diffDiv = cfg_.De * (gradXR * fR.area - gradXL * fL.area
                                        + gradYT * fT.area - gradYB * fB.area);
        b(cID) += -cfg_.e * poissonTimeStep(cc) * diffDiv;
        triplets.emplace_back(cID, cID, diag);
    }
}

void Solver::solvePoisson()
{
    using Sparse = Eigen::SparseMatrix<scalar>;
    using Triplet = Eigen::Triplet<scalar>;

    if (!hasPoissonDirichletBoundary()) {
        throw std::runtime_error("Poisson needs at least one potential boundary; pure symmetry is singular");
    }

    for (label t = 0; t < cfg_.numT; ++t) {
        std::vector<Triplet> triplets;
        Vec b = Vec::Zero(nCells_);
        triplets.reserve(static_cast<size_t>(nCells_) * 10);

        assemblePhysicalPoisson(t, triplets, b);
        {
            addSemiImplicitPoissonCorrection(t, triplets, b);
        }

        Sparse A(nCells_, nCells_);
        A.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseLU<Sparse> lu;
        lu.compute(A);
        if (lu.info() != Eigen::Success) throw std::runtime_error("Poisson factorization failed");
        const Vec x = lu.solve(b);
        if (lu.info() != Eigen::Success) throw std::runtime_error("Poisson solve failed");
        for (label cID = 0; cID < nCells_; ++cID) cells_[cID].Phi(t) = x(cID);
    }
    updateElectricField();
}

void Solver::updateElectricField()
{
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            Cell& c = cells_[cid(i, j)];
            const Face& fL = xFaces_[xf(i, j)];
            const Face& fR = xFaces_[xf(i + 1, j)];
            const Face& fB = yFaces_[yf(i, j)];
            const Face& fT = yFaces_[yf(i, j + 1)];
            for (label t = 0; t < cfg_.numT; ++t) {
                auto minusBoundaryE = [&](const Face& f) {
                    return f.kind == FaceKind::Symmetry ? 0.0 : -(c.Phi(t) - f.phiBC(t)) / f.dist;
                };
                auto plusBoundaryE = [&](const Face& f) {
                    return f.kind == FaceKind::Symmetry ? 0.0 : -(f.phiBC(t) - c.Phi(t)) / f.dist;
                };

                if (cfg_.nCellX == 1) {
                    c.Ex(t) = 0.5 * (minusBoundaryE(fL) + plusBoundaryE(fR));
                } else if (i == 0) {
                    c.Ex(t) = minusBoundaryE(fL);
                } else if (i + 1 == cfg_.nCellX) {
                    c.Ex(t) = plusBoundaryE(fR);
                }
                else {
                    const Cell& l = cells_[cid(i - 1, j)];
                    const Cell& r = cells_[cid(i + 1, j)];
                    const scalar d = dist2D({l.x, l.y}, {r.x, r.y});
                    c.Ex(t) = -(r.Phi(t) - l.Phi(t)) / std::max(d, 1.0e-30);
                }

                if (cfg_.nCellY == 1) {
                    c.Ey(t) = 0.5 * (minusBoundaryE(fB) + plusBoundaryE(fT));
                } else if (j == 0) {
                    c.Ey(t) = minusBoundaryE(fB);
                } else if (j + 1 == cfg_.nCellY) {
                    c.Ey(t) = plusBoundaryE(fT);
                }
                else {
                    const Cell& b = cells_[cid(i, j - 1)];
                    const Cell& tt = cells_[cid(i, j + 1)];
                    const scalar d = dist2D({b.x, b.y}, {tt.x, tt.y});
                    c.Ey(t) = -(tt.Phi(t) - b.Phi(t)) / std::max(d, 1.0e-30);
                }
            }
        }
    }
}

scalar Solver::poissonResidual(label cID, label t) const
{
    const label i = ci(cID);
    const label j = cj(cID);
    const Cell& c = cells_[cID];
    const Face& fL = xFaces_[xf(i, j)];
    const Face& fR = xFaces_[xf(i + 1, j)];
    const Face& fB = yFaces_[yf(i, j)];
    const Face& fT = yFaces_[yf(i, j + 1)];
    scalar r = cfg_.e * (c.Ni(t) - c.Ne(t)) * c.vol;

    auto add = [&](scalar phiNb, scalar eps, const Face& f) {
        if (f.kind == FaceKind::Symmetry) return;
        r += eps * (phiNb - c.Phi(t)) * f.area / f.dist;
    };

    if (i == 0) add(fL.phiBC(t), cfg_.eps0, fL);
    else add(cells_[cid(i - 1, j)].Phi(t), cfg_.eps0, fL);
    if (i + 1 == cfg_.nCellX) add(fR.phiBC(t), cfg_.eps0, fR);
    else add(cells_[cid(i + 1, j)].Phi(t), cfg_.eps0, fR);
    if (j > 0) add(cells_[cid(i, j - 1)].Phi(t), cfg_.eps0, fB);
    else add(fB.phiBC(t), cfg_.eps0, fB);
    if (j + 1 < cfg_.nCellY) add(cells_[cid(i, j + 1)].Phi(t), cfg_.eps0, fT);
    else add(fT.phiBC(t), cfg_.eps0, fT);
    return r;
}

}
