#include "ccp-2d.hpp"

namespace ccp2d {

ChemTable::ChemTable(const std::string& name)
{
    std::ifstream in(name);
    if (!in) throw std::runtime_error("cannot open chemistry csv: " + name);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        std::stringstream ss(line);
        std::string a, b;
        if (std::getline(ss, a, ',') && std::getline(ss, b, ',')) {
            rows_.push_back({std::stod(trim(a)), std::stod(trim(b))});
        }
    }
    if (rows_.empty()) throw std::runtime_error("empty chemistry csv: " + name);
}













scalar ChemTable::interpolate(scalar eV) const
{
    if (rows_.empty()) return 0.0;

    if (eV <= rows_.front().energy) return rows_.front().rate;
    if (eV >= rows_.back().energy) return rows_.back().rate;

    auto it = std::lower_bound(
        rows_.begin(), rows_.end(), eV,
        [](const ChemPoint& p, scalar x) { return p.energy < x; });

    const ChemPoint& hi = *it;
    const ChemPoint& lo = *(it - 1);

    const scalar rateFloor = 1.0e-60;
    const scalar k0 = std::max(lo.rate, rateFloor);
    const scalar k1 = std::max(hi.rate, rateFloor);

    const scalar w =
        (eV - lo.energy) / std::max(hi.energy - lo.energy, tiny);

    const scalar logk =
        std::log(k0) + w * (std::log(k1) - std::log(k0));

    return std::exp(logk);
}

void Solver::updateChemistryJacobian()
{
    auto meanEnergyEv = [&](scalar Ne, scalar Ee) {
        const scalar NeSafe = std::max(Ne, 1.0e-12);
        const scalar EePhys = Ee * cfg_.EeRef;
        const scalar NePhys = NeSafe * cfg_.nRef;


        return std::max(EePhys / (NePhys * cfg_.eV2J), scalar(0.0));
    };

    auto tableRate = [&](const ChemTable& table, scalar energyEv) {
        return table.interpolate(energyEv) / cfg_.klRef;
    };

    for (Cell& c : cells_)
    {
        for (label t = 0; t < cfg_.numT; ++t)
        {
            c.chemJ[t].setZero();
            const scalar Ne = c.Ne(t);
            const scalar Ee = c.Ee(t);
            const scalar energyEv = meanEnergyEv(Ne, Ee);
            const scalar kIon = tableRate(chemIon_, energyEv);
            const scalar kExc = cfg_.useExcitationLoss ? tableRate(chemExc_, energyEv) : 0.0;
            c.kIon(t) = kIon;
            c.kExc(t) = kExc;
            c.kl(t) = kIon;

            const scalar dSion_dNe = kIon * cfg_.N * c.vol;
            const scalar dSexc_dNe = kExc * cfg_.N * c.vol;
            c.chemJ[t](0, 0) = -dSion_dNe;
            c.chemJ[t](1, 0) = -dSion_dNe;
            c.chemJ[t](2, 0) = cfg_.Hion * dSion_dNe + cfg_.Hexc * dSexc_dNe;
        }
    }
}

}
