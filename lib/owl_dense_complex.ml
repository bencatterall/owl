(*
 * OWL - an OCaml math library for scientific computing
 * Copyright (c) 2016 Liang Wang <liang.wang@cl.cam.ac.uk>
 *)

(** [ Dense complex matrix ]  *)

open Bigarray
open Owl_types

type mat = Gsl.Matrix_complex.matrix

type elt = Complex.t

let shape x = (Array2.dim1 x, Array2.dim2 x)

let create m n v = Gsl.Matrix_complex.create ~init:v m n

let empty m n = Gsl.Matrix.create m n

let zeros m n = create m n Complex.zero

let pp_dsmat_complex x = Format.printf "%a\n" Owl_pretty.Toplevel.pp_cmat x
