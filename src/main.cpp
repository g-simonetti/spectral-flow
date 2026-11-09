// Adapted from Grid/tests/lanczos/Test_wilson_specflow.cc.
// Same Grid IRL calculation; adjoint links and parameters supplied by XML.
#include <Grid/Grid.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace Grid;
using FermionOp = WilsonAdjFermionD;
using FermionField = FermionOp::FermionField;

struct SpectralFlowPar : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(SpectralFlowPar,
    std::string, cfg_prefix,
    int, cfg_start,
    int, cfg_end,
    int, cfg_step,
    double, M5_start,
    double, M5_end,
    double, M5_step,
    int, Nstop,
    int, Nk,
    int, Np,
    int, MaxIt,
    double, resid,
    double, ChebyLow,
    double, ChebyHigh,
    int, ChebyOrder,
    std::vector<int>, boundary_conditions,
    std::vector<int>, rng_seeds);
};

void run(const std::string &parameterFile) {
  SpectralFlowPar p;
  XmlReader reader(parameterFile);
  read(reader, "spectralflow", p);
  if (p.cfg_prefix.empty() || p.cfg_step <= 0 || p.cfg_end < p.cfg_start ||
      !std::isfinite(p.M5_start) || !std::isfinite(p.M5_end) ||
      !std::isfinite(p.M5_step) || p.M5_step <= 0 || p.M5_end < p.M5_start ||
      p.Nstop < 1 || p.Nk < p.Nstop || p.Np < 1 || p.MaxIt < 1 ||
      !std::isfinite(p.resid) || p.resid <= 0 ||
      !std::isfinite(p.ChebyLow) || !std::isfinite(p.ChebyHigh) ||
      p.ChebyHigh <= p.ChebyLow || p.ChebyOrder < 2 ||
      p.boundary_conditions.size() != 4 || p.rng_seeds.empty())
    throw std::runtime_error("Invalid configuration range, M5 scan, Lanczos parameters, boundaries, or seeds.");
  for (int b : p.boundary_conditions)
    if (b != 1 && b != -1) throw std::runtime_error("Boundary phases must be +1 or -1.");
  std::cout << GridLogMessage << p << std::endl;
  std::cout << GridLogMessage << "Kernel: SU adjoint; Nc=" << Nc << std::endl;

  std::unique_ptr<GridCartesian> UGrid(SpaceTimeGrid::makeFourDimGrid(
    GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi()));
  std::unique_ptr<GridRedBlackCartesian> UrbGrid(SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid.get()));
  GridCartesian *FGrid = UGrid.get();
  GridRedBlackCartesian *FrbGrid = UrbGrid.get();
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(p.rng_seeds);
  LatticeGaugeField Umu(UGrid.get());
  AdjointRepresentation adj(UGrid.get());
  std::vector<Complex> boundary;
  for (int b : p.boundary_conditions) boundary.emplace_back(double(b), 0.);
  FermionOp::ImplParams impl(boundary);

  const int Nstop = p.Nstop;
  const int Nk = p.Nk;
  const int Nm = Nk + p.Np;
  const int MaxIt = p.MaxIt;
  const RealD resid = p.resid;

  for (long long cfg = p.cfg_start; cfg <= p.cfg_end; cfg += p.cfg_step) {
    FieldMetaData header;
    std::string file = p.cfg_prefix + "." + std::to_string(cfg);
    NerscIO::readConfiguration(Umu, header, file);
    adj.update_representation(Umu);

    // The original test decreases mass by 0.1. Here M5 increases by the
    // XML step, with exactly the same Wilson mass convention: mass = -M5.
    for (long long step = 0; ; ++step) {
      RealD M5 = p.M5_start + step * p.M5_step;
      if (M5 > p.M5_end + 1.e-12 * std::max(1., std::abs(p.M5_end))) break;
      RealD mass = -M5;
      FermionOp WilsonOperator(adj.U, *FGrid, *FrbGrid, mass, impl);
      MdagMLinearOperator<FermionOp,FermionField> HermOp(WilsonOperator);
      Gamma5HermitianLinearOperator<FermionOp,FermionField> HermOp2(WilsonOperator);

      Chebyshev<FermionField> Cheby(p.ChebyLow, p.ChebyHigh, p.ChebyOrder);
      FunctionHermOp<FermionField> OpCheby(Cheby,HermOp);
      PlainHermOp<FermionField> Op2(HermOp2);
      ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op2, Nstop, Nk, Nm, resid, MaxIt);

      std::vector<RealD> eval(Nm);
      FermionField src(FGrid);
      gaussian(RNG5, src);
      std::vector<FermionField> evec(Nm, FGrid);
      int Nconv;
      IRL.calc(eval, evec, src, Nconv);

      std::cout << GridLogMessage << "cfg=" << cfg << " M5=" << M5
                << " mass=" << mass << " Nconv=" << Nconv << std::endl;
      std::cout << mass << " : " << eval << std::endl;
      Gamma g5(Gamma::Algebra::Gamma5);
      FermionField tmp(FGrid);
      for (int i = 0; i < Nstop; ++i) {
        tmp = g5*evec[i];
        ComplexD dot = innerProduct(tmp,evec[i]);
        std::cout << mass << " : " << eval[i] << " "
                  << real(dot) << " " << imag(dot) << std::endl;
      }
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2 || std::string(argv[1]) == "--help") {
    std::cout << "Usage: " << argv[0] << " parameters.xml [Grid options]\n";
    return argc < 2 ? EXIT_FAILURE : EXIT_SUCCESS;
  }
  std::string parameterFile = argv[1];
  for (int i = 1; i < argc-1; ++i) argv[i] = argv[i+1];
  argv[--argc] = nullptr;
  Grid_init(&argc, &argv);
  int status = EXIT_SUCCESS;
  try { run(parameterFile); }
  catch (const std::exception &e) {
    std::cerr << "spectral_flow: " << e.what() << std::endl;
    status = EXIT_FAILURE;
  }
  Grid_finalize();
  return status;
}
