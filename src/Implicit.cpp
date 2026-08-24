#include "ccp-2d.hpp"

namespace ccp2d {


void Solver::assembleRHS_DT()
{
    const scalar invDtPhys = 1.0 / std::max(dtPhysStep_, tiny);
    const bool useBDF2 = useBDF2DT();
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            Cell& c = cells_[cid(i, j)];
            const Face& fL = xFaces_[xf(i, j)];
            const Face& fR = xFaces_[xf(i + 1, j)];
            const Face& fB = yFaces_[yf(i, j)];
            const Face& fT = yFaces_[yf(i, j + 1)];
            for (label t = 0; t < cfg_.numT; ++t) {
                const scalar netNe = fR.area * fR.flux[t].ne - fL.area * fL.flux[t].ne
                                   + fT.area * fT.flux[t].ne - fB.area * fB.flux[t].ne;
                const scalar netNi = fR.area * fR.flux[t].ni - fL.area * fL.flux[t].ni
                                   + fT.area * fT.flux[t].ni - fB.area * fB.flux[t].ni;
                const scalar netEe = fR.area * fR.flux[t].ee - fL.area * fL.flux[t].ee
                                   + fT.area * fT.flux[t].ee - fB.area * fB.flux[t].ee;
                const scalar netJ = fR.area * fR.flux[t].j - fL.area * fL.flux[t].j
                                  + fT.area * fT.flux[t].j - fB.area * fB.flux[t].j;
                const scalar Sion = c.kIon(t) * cfg_.N * c.Ne(t) * c.vol;
                const scalar Sexc = cfg_.useExcitationLoss
                                  ? c.kExc(t) * cfg_.N * c.Ne(t) * c.vol
                                  : 0.0;
                const scalar chemNe = Sion;
                const scalar chemNi = Sion;
                const scalar chemEe = -cfg_.Hion * Sion - cfg_.Hexc * Sexc;
                const scalar joule = c.Phi(t) * netJ;
                scalar physNe = -c.vol * (c.Ne(t) - c.NePhysOld(t)) * invDtPhys;
                scalar physNi = -c.vol * (c.Ni(t) - c.NiPhysOld(t)) * invDtPhys;
                scalar physEe = -c.vol * (c.Ee(t) - c.EePhysOld(t)) * invDtPhys;
                if (useBDF2) {
                    physNe = -c.vol * (1.5 * c.Ne(t) - 2.0 * c.NePhysOld(t)
                                      + 0.5 * c.NePhysOlder(t)) * invDtPhys;
                    physNi = -c.vol * (1.5 * c.Ni(t) - 2.0 * c.NiPhysOld(t)
                                      + 0.5 * c.NiPhysOlder(t)) * invDtPhys;
                    physEe = -c.vol * (1.5 * c.Ee(t) - 2.0 * c.EePhysOld(t)
                                      + 0.5 * c.EePhysOlder(t)) * invDtPhys;
                }

                const scalar rNe = -netNe + chemNe + physNe;
                const scalar rNi = -netNi + chemNi + physNi;
                const scalar rEe = -netEe - joule + chemEe + physEe;
                const scalar coeff = c.dt / std::max(c.vol, tiny);
                c.rhsNe(t) = coeff * rNe;
                c.rhsNi(t) = coeff * rNi;
                c.rhsEe(t) = coeff * rEe;
                c.resNe(t) = rNe;
                c.resNi(t) = rNi;
                c.resEe(t) = rEe;
                c.resPhi(t) = poissonResidual(cid(i, j), t);
            }
        }
    }
}

Mat3 Solver::buildBlockB(label cID, label t) const
{
    const Cell& c = cells_[cID];
    const label i = ci(cID);
    const label j = cj(cID);
    const Face& fL = xFaces_[xf(i, j)];
    const Face& fR = xFaces_[xf(i + 1, j)];
    const Face& fB = yFaces_[yf(i, j)];
    const Face& fT = yFaces_[yf(i, j + 1)];

    const scalar cD = cfg_.cD;
    const std::array<const Face*, 4> faces{{&fL, &fR, &fB, &fT}};
    scalar lamE = 0.0;
    scalar lamI = 0.0;
    scalar lamEe = 0.0;
    for (const Face* fp : faces) {
        const Face& f = *fp;
        if (f.kind == FaceKind::Symmetry) continue;
        const scalar areaOverVol = f.area / std::max(c.vol, tiny);
        const scalar h = std::max(f.dist, 1.0e-30);
        const scalar absEf = std::abs(f.flux[t].ef);
        lamE += areaOverVol * (0.5 * cfg_.muE * absEf + cD * cfg_.De / h);
        lamI += areaOverVol * (0.5 * cfg_.muI * absEf + cD * cfg_.Di / h);
        lamEe += areaOverVol * ((5.0 / 6.0) * cfg_.muE * absEf
                              + (5.0 / 3.0) * cD * cfg_.De / h);
    }
    const scalar coeff = cfg_.epsImplicit * c.dt / std::max(c.vol, tiny);

    Mat3 B = Mat3::Identity();
    B(0, 0) += cfg_.epsImplicit * c.dt * lamE;
    B(1, 1) += cfg_.epsImplicit * c.dt * lamI;
    B(2, 2) += cfg_.epsImplicit * c.dt * lamEe;
    B.diagonal().array() += cfg_.epsImplicit * c.dt * physicalTimeDerivativeCoeff();
    const Mat3 chemJ = c.chemJ[t].topLeftCorner<3, 3>();
    B += coeff * chemJ;
    for (label row = 0; row < 3; ++row) {
        B(row, row) += coeff * cfg_.chemRowSumDamping * chemJ.row(row).cwiseAbs().sum();
    }
    return B;
}

void Solver::forwardSweepImplicit()
{

    for (Cell& c : cells_) {
        c.dQhNe.setZero();
        c.dQhNi.setZero();
        c.dQhEe.setZero();
    }

    for (label cID = 0; cID < nCells_; ++cID) {
        const label i = ci(cID);
        const label j = cj(cID);
        Cell& c = cells_[cID];
        for (label t = 0; t < cfg_.numT; ++t) {
            Vec3 rhs(c.rhsNe(t), c.rhsNi(t), c.rhsEe(t));
            const scalar coeff = cfg_.epsImplicit * c.dt / std::max(c.vol, tiny);

            if (i > 0) {
                const label nb = cid(i - 1, j);
                const Face& f = xFaces_[xf(i, j)];
                Vec3 dq(cells_[nb].dQhNe(t), cells_[nb].dQhNi(t), cells_[nb].dQhEe(t));
                const Mat3 A = f.area * f.JL[t].topLeftCorner<3, 3>()
                             - c.jouleL[t].topLeftCorner<3, 3>();
                rhs += coeff * A * dq;
            }

            if (j > 0) {
                const label nb = cid(i, j - 1);
                const Face& f = yFaces_[yf(i, j)];
                Vec3 dq(cells_[nb].dQhNe(t), cells_[nb].dQhNi(t), cells_[nb].dQhEe(t));
                const Mat3 A = f.area * f.JL[t].topLeftCorner<3, 3>()
                             - c.jouleB[t].topLeftCorner<3, 3>();
                rhs += coeff * A * dq;
            }

            const Vec3 dq = buildBlockB(cID, t).partialPivLu().solve(rhs);
            c.dQhNe(t) = dq(0);
            c.dQhNi(t) = dq(1);
            c.dQhEe(t) = dq(2);
        }
    }
}

void Solver::backwardSweepImplicit()
{

    for (Cell& c : cells_) {
        c.dQNe.setZero();
        c.dQNi.setZero();
        c.dQEe.setZero();
    }

    for (label cID = nCells_ - 1; cID >= 0; --cID) {
        const label i = ci(cID);
        const label j = cj(cID);
        Cell& c = cells_[cID];
        for (label t = 0; t < cfg_.numT; ++t) {
            const scalar coeff = cfg_.epsImplicit * c.dt / std::max(c.vol, tiny);
            Vec3 down = Vec3::Zero();

            if (i + 1 < cfg_.nCellX) {
                const label nb = cid(i + 1, j);
                const Face& f = xFaces_[xf(i + 1, j)];
                Vec3 dq(cells_[nb].dQNe(t), cells_[nb].dQNi(t), cells_[nb].dQEe(t));
                const Mat3 A = f.area * f.JR[t].topLeftCorner<3, 3>()
                             + c.jouleR[t].topLeftCorner<3, 3>();
                down += coeff * A * dq;
            }

            if (j + 1 < cfg_.nCellY) {
                const label nb = cid(i, j + 1);
                const Face& f = yFaces_[yf(i, j + 1)];
                Vec3 dq(cells_[nb].dQNe(t), cells_[nb].dQNi(t), cells_[nb].dQEe(t));
                const Mat3 A = f.area * f.JR[t].topLeftCorner<3, 3>()
                             + c.jouleT[t].topLeftCorner<3, 3>();
                down += coeff * A * dq;
            }

            const Vec3 half(c.dQhNe(t), c.dQhNi(t), c.dQhEe(t));
            const Vec3 dq = half - buildBlockB(cID, t).partialPivLu().solve(down);
            c.dQNe(t) = dq(0);
            c.dQNi(t) = dq(1);
            c.dQEe(t) = dq(2);
        }
        if (cID == 0) break;
    }
}

void Solver::updateSolution()
{
    for (Cell& c : cells_) {
        for (label t = 0; t < cfg_.numT; ++t) {
            const scalar fNe = std::max(cfg_.qFloor, cfg_.relFloor * std::abs(c.NeOld(t)));
            const scalar fNi = std::max(cfg_.qFloor, cfg_.relFloor * std::abs(c.NiOld(t)));
            const scalar fEe = std::max(cfg_.qFloor, cfg_.relFloor * std::abs(c.EeOld(t)));
            const scalar NeOld = std::max(c.NeOld(t), fNe);
            const scalar NiOld = std::max(c.NiOld(t), fNi);
            const scalar EeOld = std::max(c.EeOld(t), fEe);
            const scalar dNe = c.dQNe(t);
            const scalar dNi = c.dQNi(t);
            const scalar dEe = c.dQEe(t);

            scalar alpha = cfg_.omegaImplicit;
            auto positivity = [&](scalar qOld, scalar dQ, scalar qMin) {
                if (dQ < 0.0) {
                    alpha = std::min(alpha, 0.9 * std::max(qOld - qMin, scalar(0.0)) / (-dQ));
                }
            };
            auto relLimit = [&](scalar qOld, scalar dQ, scalar qFloor) {
                if (std::abs(dQ) > tiny) {
                    const scalar qRef = std::max(std::abs(qOld), qFloor);
                    alpha = std::min(alpha, cfg_.dQmaxRel * qRef / std::abs(dQ));
                }
            };

            positivity(NeOld, dNe, fNe);
            positivity(NiOld, dNi, fNi);
            positivity(EeOld, dEe, fEe);
            relLimit(NeOld, dNe, fNe);
            relLimit(NiOld, dNi, fNi);
            relLimit(EeOld, dEe, fEe);
            alpha = std::max(scalar(0.0), std::min(alpha, cfg_.omegaImplicit));

            scalar Ne = std::max(fNe, NeOld + alpha * dNe);
            scalar Ni = std::max(fNi, NiOld + alpha * dNi);
            scalar Ee = std::max(fEe, EeOld + alpha * dEe);
            scalar Te = Ee / std::max(Ne, cfg_.qFloor);
            if (cfg_.TeMin > 0.0 && Te < cfg_.TeMin) {
                Ee = std::max(fEe, Ne * cfg_.TeMin);
                Te = Ee / std::max(Ne, cfg_.qFloor);
            }
            if (std::isfinite(cfg_.TeMax) && Te > cfg_.TeMax) {
                Ee = std::max(fEe, Ne * cfg_.TeMax);
                Te = Ee / std::max(Ne, cfg_.qFloor);
            }

            c.Ne(t) = Ne;
            c.Ni(t) = Ni;
            c.Ee(t) = Ee;
            c.Te(t) = Te;
        }
    }
}

void Solver::computeDt()
{
    dtMin_ = std::numeric_limits<scalar>::max();
    const scalar cD = cfg_.cD;
    const scalar cfl = cfg_.pseudoCFL;
    for (label cID = 0; cID < nCells_; ++cID) {
        Cell& c = cells_[cID];
        const label i = ci(cID);
        const label j = cj(cID);
        const Face& fL = xFaces_[xf(i, j)];
        const Face& fR = xFaces_[xf(i + 1, j)];
        const Face& fB = yFaces_[yf(i, j)];
        const Face& fT = yFaces_[yf(i, j + 1)];
        const std::array<const Face*, 4> faces{{&fL, &fR, &fB, &fT}};
        scalar lamE = 0.0;
        scalar lamI = 0.0;
        scalar lamEe = 0.0;
        scalar NeMax = 0.0;
        scalar NiMax = 0.0;
        for (label t = 0; t < cfg_.numT; ++t) {
            NeMax = std::max(NeMax, c.Ne(t));
            NiMax = std::max(NiMax, c.Ni(t));
            scalar lamEt = 0.0;
            scalar lamIt = 0.0;
            scalar lamEet = 0.0;
            for (const Face* fp : faces) {
                const Face& f = *fp;
                if (f.kind == FaceKind::Symmetry) continue;
                const scalar areaOverVol = f.area / std::max(c.vol, tiny);
                const scalar h = std::max(f.dist, 1.0e-30);
                const scalar absEf = std::abs(f.flux[t].ef);
                lamEt += areaOverVol * (0.5 * cfg_.muE * absEf + cD * cfg_.De / h);
                lamIt += areaOverVol * (0.5 * cfg_.muI * absEf + cD * cfg_.Di / h);
                lamEet += areaOverVol * ((5.0 / 6.0) * cfg_.muE * absEf
                                        + (5.0 / 3.0) * cD * cfg_.De / h);
            }
            lamE = std::max(lamE, lamEt);
            lamI = std::max(lamI, lamIt);
            lamEe = std::max(lamEe, lamEet);
        }
        const scalar invTauM = cfg_.e * (cfg_.muE * NeMax + cfg_.muI * NiMax)
                             / std::max(cfg_.eps0, tiny);
        const scalar specRad = std::max(lamEe, std::max(lamE, lamI));
        c.dt = cfl / (specRad + invTauM + tiny);
        dtMin_ = std::min(dtMin_, c.dt);
    }
    for (Cell& c : cells_) c.dt = dtMin_;
}

bool Solver::innerResidualConverged(const std::array<scalar, 4>& rr) const
{
    if (cfg_.innerResTol <= 0.0) return false;
    const scalar target = std::log10(std::max(cfg_.innerResTol, 1.0e-300));
    return std::max(std::max(rr[0], rr[1]), std::max(rr[2], rr[3])) <= target;
}

bool Solver::innerRelativeResidualConverged(const std::array<scalar, 5>& g,
                                            const std::array<scalar, 5>& g0) const
{
    if (cfg_.innerRelResTol <= 0.0) return false;
    const scalar denom = std::max(g0[0], 1.0e-300);
    return g[0] / denom <= cfg_.innerRelResTol;
}

}
