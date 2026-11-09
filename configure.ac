AC_PREREQ([2.63])
AC_INIT([spectral-flow], [0.1.0])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_ARG_WITH([grid],
  [AS_HELP_STRING([--with-grid=PREFIX], [installed Grid prefix containing bin/grid-config])],
  [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"])
AC_PATH_PROG([GRID_CONFIG], [grid-config])
AS_IF([test -z "$GRID_CONFIG"], [AC_MSG_ERROR([grid-config not found; use --with-grid=PREFIX])])
AS_IF([test -z "$CXX"], [CXX=`"$GRID_CONFIG" --cxx`])
AC_PROG_CXX
AC_LANG([C++])
AC_ARG_VAR([CXXLD], [C++ linker command (defaults to grid-config --cxxld)])
AS_IF([test -z "$CXXLD"], [CXXLD=`"$GRID_CONFIG" --cxxld`])
AC_SUBST([CXXLD])
CXXFLAGS="$CXXFLAGS `"$GRID_CONFIG" --cxxflags`"
LDFLAGS="$LDFLAGS `"$GRID_CONFIG" --ldflags`"
LIBS="`"$GRID_CONFIG" --libs` $LIBS"
AC_MSG_CHECKING([whether Grid with double-precision adjoint Wilson support links])
AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <Grid/Grid.h>
int main(int argc, char **argv) {
  Grid::Grid_init(&argc, &argv);
  using Action = Grid::WilsonAdjFermionD;
  auto *g = Grid::SpaceTimeGrid::makeFourDimGrid(Grid::GridDefaultLatt(),
    Grid::GridDefaultSimd(Grid::Nd, Grid::vComplexD::Nsimd()), Grid::GridDefaultMpi());
  auto *rb = Grid::SpaceTimeGrid::makeFourDimRedBlackGrid(g);
  Grid::LatticeGaugeField u(g);
  Grid::SU<Grid::Nc>::ColdConfiguration(u);
  Grid::AdjointRepresentation adj(g);
  adj.update_representation(u);
  Action d(adj.U, *g, *rb, -1.0);
  Action::FermionField in(g), out(g);
  in = Grid::Zero();
  d.M(in, out);
  Grid::Grid_finalize();
  return 0;
}
]])], [AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])
 AC_MSG_ERROR([Cannot link adjoint Wilson Grid program; inspect config.log and ensure Grid was built with fermion representations.])])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
