// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Copyright (c) 2018 by Norwegian University of Science and Technology
// PythonBindings.cpp
// This file is part of the Garamon for c4ga.
// Authors: Stephane Breuils, Vincent Nozick, and Lars Tingelstad
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file PythonBindings.cpp
/// \author Lars Tingelstad
/// \brief Python bindings using pybind11.

#include "c4ga/Mvec.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

namespace py = pybind11;


/*!
 * @namespace c4ga
 */
namespace c4ga {

PYBIND11_MODULE(c4ga_py, m) {

  m.attr("E0") = 1;
m.attr("E1") = 2;
m.attr("E2") = 4;
m.attr("E3") = 8;
m.attr("E4") = 16;
m.attr("Ei") = 32;
m.attr("E01") = 3;
m.attr("E02") = 5;
m.attr("E03") = 9;
m.attr("E04") = 17;
m.attr("E0i") = 33;
m.attr("E12") = 6;
m.attr("E13") = 10;
m.attr("E14") = 18;
m.attr("E1i") = 34;
m.attr("E23") = 12;
m.attr("E24") = 20;
m.attr("E2i") = 36;
m.attr("E34") = 24;
m.attr("E3i") = 40;
m.attr("E4i") = 48;
m.attr("E012") = 7;
m.attr("E013") = 11;
m.attr("E014") = 19;
m.attr("E01i") = 35;
m.attr("E023") = 13;
m.attr("E024") = 21;
m.attr("E02i") = 37;
m.attr("E034") = 25;
m.attr("E03i") = 41;
m.attr("E04i") = 49;
m.attr("E123") = 14;
m.attr("E124") = 22;
m.attr("E12i") = 38;
m.attr("E134") = 26;
m.attr("E13i") = 42;
m.attr("E14i") = 50;
m.attr("E234") = 28;
m.attr("E23i") = 44;
m.attr("E24i") = 52;
m.attr("E34i") = 56;
m.attr("E0123") = 15;
m.attr("E0124") = 23;
m.attr("E012i") = 39;
m.attr("E0134") = 27;
m.attr("E013i") = 43;
m.attr("E014i") = 51;
m.attr("E0234") = 29;
m.attr("E023i") = 45;
m.attr("E024i") = 53;
m.attr("E034i") = 57;
m.attr("E1234") = 30;
m.attr("E123i") = 46;
m.attr("E124i") = 54;
m.attr("E134i") = 58;
m.attr("E234i") = 60;
m.attr("E01234") = 31;
m.attr("E0123i") = 47;
m.attr("E0124i") = 55;
m.attr("E0134i") = 59;
m.attr("E0234i") = 61;
m.attr("E1234i") = 62;
m.attr("E01234i") = 63;

  
  m.def("metric", [](){return metric;});

  m.attr("scalar") = 0;
  m.def("e0", &e0<double>);
m.def("e1", &e1<double>);
m.def("e2", &e2<double>);
m.def("e3", &e3<double>);
m.def("e4", &e4<double>);
m.def("ei", &ei<double>);
m.def("e01", &e01<double>);
m.def("e02", &e02<double>);
m.def("e03", &e03<double>);
m.def("e04", &e04<double>);
m.def("e0i", &e0i<double>);
m.def("e12", &e12<double>);
m.def("e13", &e13<double>);
m.def("e14", &e14<double>);
m.def("e1i", &e1i<double>);
m.def("e23", &e23<double>);
m.def("e24", &e24<double>);
m.def("e2i", &e2i<double>);
m.def("e34", &e34<double>);
m.def("e3i", &e3i<double>);
m.def("e4i", &e4i<double>);
m.def("e012", &e012<double>);
m.def("e013", &e013<double>);
m.def("e014", &e014<double>);
m.def("e01i", &e01i<double>);
m.def("e023", &e023<double>);
m.def("e024", &e024<double>);
m.def("e02i", &e02i<double>);
m.def("e034", &e034<double>);
m.def("e03i", &e03i<double>);
m.def("e04i", &e04i<double>);
m.def("e123", &e123<double>);
m.def("e124", &e124<double>);
m.def("e12i", &e12i<double>);
m.def("e134", &e134<double>);
m.def("e13i", &e13i<double>);
m.def("e14i", &e14i<double>);
m.def("e234", &e234<double>);
m.def("e23i", &e23i<double>);
m.def("e24i", &e24i<double>);
m.def("e34i", &e34i<double>);
m.def("e0123", &e0123<double>);
m.def("e0124", &e0124<double>);
m.def("e012i", &e012i<double>);
m.def("e0134", &e0134<double>);
m.def("e013i", &e013i<double>);
m.def("e014i", &e014i<double>);
m.def("e0234", &e0234<double>);
m.def("e023i", &e023i<double>);
m.def("e024i", &e024i<double>);
m.def("e034i", &e034i<double>);
m.def("e1234", &e1234<double>);
m.def("e123i", &e123i<double>);
m.def("e124i", &e124i<double>);
m.def("e134i", &e134i<double>);
m.def("e234i", &e234i<double>);
m.def("e01234", &e01234<double>);
m.def("e0123i", &e0123i<double>);
m.def("e0124i", &e0124i<double>);
m.def("e0134i", &e0134i<double>);
m.def("e0234i", &e0234i<double>);
m.def("e1234i", &e1234i<double>);
m.def("e01234i", &e01234i<double>);

  m.def("I", &I<double>);


  // Class definition
  auto mvec = py::class_<Mvec<double>>(m, "Mvec");
  // Constructors
  mvec.def(py::init<>());
  // Get/Set
  mvec.def("__setitem__",
           [](Mvec<double>& mv, int idx, double value) { mv[idx] = value; });
  mvec.def("__getitem__",
           [](Mvec<double>& mv, int idx) { return mv[idx]; });
  // Operators
  mvec.def(py::self + py::self)
      .def(py::self + float())
      .def(float() + py::self)
      .def(py::self += py::self)
      .def(py::self - py::self)
      .def(py::self - float())
      .def(float() - py::self)
      .def(py::self -= py::self)
      .def(py::self * py::self)
      .def(py::self * float())
      .def(float() * py::self)
      .def(py::self *= py::self)
      .def(py::self / py::self)
      .def(py::self / float())
      .def(float() / py::self)
      .def(py::self /= py::self)
      .def(py::self < py::self)
      .def(py::self < float())
      .def(float() < py::self)
      .def(py::self > py::self)
      .def(py::self > float())
      .def(float() > py::self)
      .def("__invert__",
           [](const Mvec<double>& a) { return ~a; })
      .def("__eq__",
           [](Mvec<double>& a, const Mvec<double>& b) { return a == b; })
      .def("__neq__",
           [](Mvec<double>& a, const Mvec<double>& b) { return a != b; })
      .def("__or__",
           [](const Mvec<double>& a, const Mvec<double>& b) { return a | b; })
      .def("__or__",
           [](const Mvec<double>& a, double b) { return a | b; })
      .def("__ror__",
           [](const Mvec<double>& a, double b) { return a | b; })
      .def("__ior__",
           [](Mvec<double>& a, const Mvec<double>& b) { a |= b; return a; })
      .def("__xor__",
           [](const Mvec<double>& a, const Mvec<double>& b) { return a ^ b; })
      .def("__xor__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__rxor__",
           [](const Mvec<double>& a, double b) { return a ^ b; })
      .def("__ixor__",
           [](Mvec<double>& a, const Mvec<double>& b) { a ^= b; return a; });

  // Print
  mvec.def("__repr__",
           [](const Mvec<double>& mv) {
             std::stringstream ss;
             ss << mv;
             return ss.str();
           })
      .def("norm", &Mvec<double>::norm)
      .def("quadratic_norm", &Mvec<double>::quadraticNorm)
      .def("reverse", &Mvec<double>::reverse)
      .def("display", &Mvec<double>::display,
        py::call_guard<py::scoped_ostream_redirect,
                       py::scoped_estream_redirect>());


  mvec.def("outer_primal_dual", &Mvec<double>::outerPrimalDual);
  mvec.def("outer_dual_primal", &Mvec<double>::outerDualPrimal);
  mvec.def("outer_dual_dual", &Mvec<double>::outerDualDual);
  mvec.def("dual", &Mvec<double>::dual);


  mvec.def("scalar_product", &Mvec<double>::scalarProduct);
  mvec.def("dot_product", &Mvec<double>::dotProduct);
  mvec.def("inv", &Mvec<double>::inv);

  mvec.def("grades", &Mvec<double>::grades);
  mvec.def("grade", [](const Mvec<double>& a){return a.grade();});
  mvec.def("grade", [](const Mvec<double>& a, const int i){return a.grade(i);});
  mvec.def("clear", &Mvec<double>::clear);

}

}  // namespace c4ga
