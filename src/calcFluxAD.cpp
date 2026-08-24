#include "ccp-2d.hpp"

namespace ccp2d {

void Solver::computeFaceAD(Face& f, label t)
{
    f.JL[t].setZero();
    f.JR[t].setZero();
    f.flux[t] = Flux{};
    if (f.kind == FaceKind::Symmetry) return;

    auto copyResult = [&](const auto& r) {
        f.flux[t].ne = r.ne.val();
        f.flux[t].ni = r.ni.val();
        f.flux[t].ee = r.ee.val();
        f.flux[t].j = r.j.val();
        f.flux[t].phi = r.phi.val();
        f.flux[t].ef = r.ef.val();
    };

    if (f.kind == FaceKind::MinusBC) {
        const Cell& R = cells_[f.cR];
        constexpr label nDofs = 4;
        ADScalar NeR(nDofs, 0, R.Ne(t));
        ADScalar NiR(nDofs, 1, R.Ni(t));
        ADScalar EeR(nDofs, 2, R.Ee(t));
        ADScalar PhiR(nDofs, 3, R.Phi(t));
        const auto r = minusBCFluxKernel(NeR, NiR, EeR, PhiR,
                                         f.phiBC(t), f.wallModel, f.wallTe, f.dist, cfg_);
        copyResult(r);
        for (label k = 0; k < nDofs; ++k) {
            f.JR[t](0, k) = r.ne.dx(k);
            f.JR[t](1, k) = r.ni.dx(k);
            f.JR[t](2, k) = r.ee.dx(k);
            f.JR[t](3, k) = r.phi.dx(k);
        }
        return;
    }

    if (f.kind == FaceKind::PlusBC) {
        const Cell& L = cells_[f.cL];
        constexpr label nDofs = 4;
        ADScalar NeL(nDofs, 0, L.Ne(t));
        ADScalar NiL(nDofs, 1, L.Ni(t));
        ADScalar EeL(nDofs, 2, L.Ee(t));
        ADScalar PhiL(nDofs, 3, L.Phi(t));
        const auto r = plusBCFluxKernel(NeL, NiL, EeL, PhiL,
                                        f.phiBC(t), f.wallModel, f.wallTe, f.dist, cfg_);
        copyResult(r);
        for (label k = 0; k < nDofs; ++k) {
            f.JL[t](0, k) = r.ne.dx(k);
            f.JL[t](1, k) = r.ni.dx(k);
            f.JL[t](2, k) = r.ee.dx(k);
            f.JL[t](3, k) = r.phi.dx(k);
        }
        return;
    }

    const Cell& L = cells_[f.cL];
    const Cell& R = cells_[f.cR];
    constexpr label nDofs = 8;
    ADScalar NeL(nDofs, 0, L.Ne(t));
    ADScalar NiL(nDofs, 1, L.Ni(t));
    ADScalar EeL(nDofs, 2, L.Ee(t));
    ADScalar PhiL(nDofs, 3, L.Phi(t));
    ADScalar NeR(nDofs, 4, R.Ne(t));
    ADScalar NiR(nDofs, 5, R.Ni(t));
    ADScalar EeR(nDofs, 6, R.Ee(t));
    ADScalar PhiR(nDofs, 7, R.Phi(t));
    const auto r = internalFluxKernel(NeL, NiL, EeL, PhiL, NeR, NiR, EeR, PhiR, f.dist, cfg_);
    copyResult(r);
    for (label k = 0; k < 4; ++k) {
        f.JL[t](0, k) = r.ne.dx(k);
        f.JL[t](1, k) = r.ni.dx(k);
        f.JL[t](2, k) = r.ee.dx(k);
        f.JL[t](3, k) = r.phi.dx(k);
        f.JR[t](0, k) = r.ne.dx(k + 4);
        f.JR[t](1, k) = r.ni.dx(k + 4);
        f.JR[t](2, k) = r.ee.dx(k + 4);
        f.JR[t](3, k) = r.phi.dx(k + 4);
    }
}

void Solver::computeFluxesAD()
{
    for (Face& f : xFaces_) {
        for (label t = 0; t < cfg_.numT; ++t) computeFaceAD(f, t);
    }
    for (Face& f : yFaces_) {
        for (label t = 0; t < cfg_.numT; ++t) computeFaceAD(f, t);
    }
}

void Solver::addJRow(Mat4& M, scalar scale, const Mat4& J, label t)
{
    (void)t;
    for (label k = 0; k < 4; ++k) M(2, k) += scale * cfg_.e * J(0, k);
}

void Solver::computeJouleJacobianAD()
{
    for (label j = 0; j < cfg_.nCellY; ++j) {
        for (label i = 0; i < cfg_.nCellX; ++i) {
            Cell& c = cells_[cid(i, j)];
            const Face& fL = xFaces_[xf(i, j)];
            const Face& fR = xFaces_[xf(i + 1, j)];
            const Face& fB = yFaces_[yf(i, j)];
            const Face& fT = yFaces_[yf(i, j + 1)];
            for (label t = 0; t < cfg_.numT; ++t) {
                c.jouleC[t].setZero();
                c.jouleL[t].setZero();
                c.jouleR[t].setZero();
                c.jouleB[t].setZero();
                c.jouleT[t].setZero();

                const scalar phi = c.Phi(t);
                const scalar netJ = fR.area * fR.flux[t].j - fL.area * fL.flux[t].j
                                  + fT.area * fT.flux[t].j - fB.area * fB.flux[t].j;
                c.jouleC[t](2, 3) += netJ;

                if (i > 0) {
                    addJRow(c.jouleL[t], -phi * fL.area, fL.JL[t], t);
                    addJRow(c.jouleC[t], -phi * fL.area, fL.JR[t], t);
                } else {
                    addJRow(c.jouleC[t], -phi * fL.area, fL.JR[t], t);
                }

                if (i + 1 < cfg_.nCellX) {
                    addJRow(c.jouleC[t], phi * fR.area, fR.JL[t], t);
                    addJRow(c.jouleR[t], phi * fR.area, fR.JR[t], t);
                } else {
                    addJRow(c.jouleC[t], phi * fR.area, fR.JL[t], t);
                }

                if (j > 0) {
                    addJRow(c.jouleB[t], -phi * fB.area, fB.JL[t], t);
                    addJRow(c.jouleC[t], -phi * fB.area, fB.JR[t], t);
                } else {
                    addJRow(c.jouleC[t], -phi * fB.area, fB.JR[t], t);
                }

                if (j + 1 < cfg_.nCellY) {
                    addJRow(c.jouleC[t], phi * fT.area, fT.JL[t], t);
                    addJRow(c.jouleT[t], phi * fT.area, fT.JR[t], t);
                } else {
                    addJRow(c.jouleC[t], phi * fT.area, fT.JL[t], t);
                }
            }
        }
    }
}

}
