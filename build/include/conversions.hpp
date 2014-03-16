#ifndef __CONVERSIONS__
#define __CONVERSIONS__

/* Conversion factors to and from atomic units (a.u. or AU)
 *
 * Sources:
 * Weinhold, F.; Landis, C.R. Discovering Chemistry with Natural Bond Orbitals, Appendix E
 * http://en.wikipedia.org/wiki/Atomic_units, retrieved 2013-09-23
 *
 * USAGE
 * These are multiplicative factors.
 * To convert from a.u. (Hartree) to eV, for instance:
 * energy_ev = energy_au * AU2EV;
 */

// Energy (Hartree in eV)
const double AU2EV = 27.2114;
const double EV2AU = 1.0/AU2EV;

// Mass (m_e in kg)
const double AU2KG = 9.10939e-31;
const double KG2AU = 1.0/AU2KG;

// Charge (q in Coulombs)
const double AU2COUL = 1.602188e-19;
const double COUL2AU = 1.0/AU2COUL;

// Angular momentum (hbar in J*s)
const double AU2ANGMOM = 1.05457e-34;
const double ANGMOM2AU = 1.0/AU2ANGMOM;

// Length (Bohr radius in m)
const double AU2M = 5.29177e-11;
const double M2AU = 1.0/AU2M;

// Time (s)
const double AU2S = 2.41888e-17;
const double S2AU = 1.0/AU2S;

// Electric dipole moment (C*m)
const double AU2CM = 8.47836-30;
const double CM2AU = 1.0/AU2CM;

// Magnetic dipole moment (Bohr magneton J/T)
const double AU2BOHRMAG = 9.27402e-24;
const double BOHRMAG2AU = 1.0/AU2BOHRMAG;

// Temperature (K)
const double AU2K = 3.157746455e5;
const double K2AU = 1.0/AU2K;

// Velocity (m/s)
const double AU2VEL = 2.187691263373e6;
const double VEL2AU = 1.0/AU2VEL;

// Pressure (Pa)
const double AU2PA = 2.942191219e13;
const double PA2AU = 1.0/AU2PA;

// Electric Field (V/m)
const double AU2EFIELD = 5.1422065211e11;
const double EFIELD2AU = 1.0/AU2EFIELD;

#endif
