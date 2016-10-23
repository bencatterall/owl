(*
 * OWL - an OCaml math library for scientific computing
 * Copyright (c) 2016 Liang Wang <liang.wang@cl.cam.ac.uk>
 *)

(** [ Dense complex matrix ]  *)

open Bigarray
open Owl_types


module CommonImpl = struct

  type elt = Complex.t

  type prc = complex64_elt

  type mat = Gsl.Matrix_complex.matrix

  let const_0 = Complex.zero

  let const_1 = Complex.one

  let empty m n = Gsl.Matrix_complex.create m n

  let create m n v = Gsl.Matrix_complex.create ~init:v m n

end

include CommonImpl
include Owl_matrix.Common (CommonImpl)

let sequential m n =
  let x = empty m n and c = ref const_0 in
  for i = 0 to m - 1 do
    for j = 0 to n - 1 do
      c := Complex.(add !c const_1);
      x.{i,j} <- !c
    done
  done; x

let linspace a b n =
  let x = empty 1 n in
  let c = ((b -. a) /. (float_of_int (n - 1))) in
  for i = 0 to n - 1 do
    let d = a +. c *. (float_of_int i) in
     x.{0,i} <- Complex.({re = d; im = 0.})
  done; x


let pp_dsmat_complex x = Format.printf "%a\n" Owl_pretty.Toplevel.pp_cmat x
