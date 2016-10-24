(*
 * OWL - an OCaml math library for scientific computing
 * Copyright (c) 2016 Liang Wang <liang.wang@cl.cam.ac.uk>
 *)

(** [ Complex dense matrix ]  *)

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

  let swap_rows = Gsl.Matrix_complex.swap_rows

  let swap_cols = Gsl.Matrix_complex.swap_columns

  let swap_rowcol = Gsl.Matrix_complex.swap_rowcol

  let transpose x =
    let y = empty (Array2.dim2 x) (Array2.dim1 x) in
    Gsl.Matrix_complex.transpose y x; y

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


(* matrix mathematical operations *)

let add x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix_complex.add x3 x2; x3

let sub x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix_complex.sub x3 x2; x3

let mul x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix_complex.mul_elements x3 x2; x3

let div x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix_complex.div_elements x3 x2; x3

let dot x1 x2 =
  let open Gsl.Blas.Complex in
  let x3 = empty (row_num x1) (col_num x2) in
  let _ = gemm ~ta:NoTrans ~tb:NoTrans ~alpha:const_1 ~beta:const_0 ~a:x1 ~b:x2 ~c:x3
  in x3

let power x c = map (fun y -> Complex.pow y c) x

let abs x = map (fun y -> Complex.({re = norm y; im = 0.})) x

let neg x =
  let y = clone x in
  Gsl.Matrix_complex.scale y Complex.({re = -1.; im = -1.}); y

let sum x =
  let y = ones 1 (row_num x) in
  let z = ones (col_num x) 1 in
  (dot (dot y x) z).{0,0}

let sum_cols x =
  let y = ones (col_num x) 1 in
  dot x y

let sum_rows x =
  let y = ones 1 (row_num x) in
  dot y x

let average x =
  let c = float_of_int (numel x) in
  Complex.(div (sum x) {re = c; im = 0.})

let average_cols x =
  let m, n = shape x in
  let c = 1. /. (float_of_int n) in
  let y = create n 1 Complex.({re = c; im = 0.}) in
  dot x y

let average_rows x =
  let m, n = shape x in
  let c = 1. /. (float_of_int m) in
  let y = create 1 m Complex.({re = c; im = 0.}) in
  dot y x

let ( +@ ) = add

let ( -@ ) = sub

let ( *@ ) = mul

let ( /@ ) = div

let ( $@ ) = dot

let ( **@ ) = power

let ( +$ ) x a =
  let y = clone x in
  Gsl.Matrix_complex.add_constant y a; y

let ( $+ ) a x = ( +$ ) x a

let ( -$ ) x a = ( +$ ) x (Complex.neg a)

let ( $- ) a x = ( -$ ) x a

let ( *$ ) x a =
  let y = clone x in
  Gsl.Matrix_complex.scale y a; y

let ( $* ) a x = ( *$ ) x a

let ( /$ ) x a = ( *$ ) x (Complex.inv a)

let ( $/ ) a x = ( /$ ) x a

let add_scalar = ( +$ )

let sub_scalar = ( -$ )

let mul_scalar = ( *$ )

let div_scalar = ( /$ )



let pp_dsmat_complex x = Format.printf "%a\n" Owl_pretty.Toplevel.pp_cmat x
