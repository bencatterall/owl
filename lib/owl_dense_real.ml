(*
 * OWL - an OCaml math library for scientific computing
 * Copyright (c) 2016 Liang Wang <liang.wang@cl.cam.ac.uk>
 *)

(** [ Dense matrix ]  *)

open Bigarray
open Owl_types


module CommonImpl = struct

  type elt = float

  type prc = float64_elt

  type mat = Gsl.Matrix.matrix

  let const_0 = 0.

  let const_1 = 1.

  let empty m n = Gsl.Matrix.create m n

  let create m n v = Gsl.Matrix.create ~init:v m n

  let swap_rows = Gsl.Matrix.swap_rows

  let swap_cols = Gsl.Matrix.swap_columns

  let swap_rowcol = Gsl.Matrix.swap_rowcol

  let transpose x =
    let y = empty (Array2.dim2 x) (Array2.dim1 x) in
    Gsl.Matrix.transpose y x; y

end

include CommonImpl
include Owl_matrix.Common (CommonImpl)


(* matrix creation operations *)

let sequential m n =
  let x = empty m n and c = ref 0 in
  for i = 0 to m - 1 do
    for j = 0 to n - 1 do
      c := !c + 1; x.{i,j} <- (float_of_int !c)
    done
  done; x

let linspace a b n =
  let x = empty 1 n in
  let c = ((b -. a) /. (float_of_int (n - 1))) in
  for i = 0 to n - 1 do
    x.{0,i} <- a +. c *. (float_of_int i)
  done; x


(* matrix mathematical operations *)

let add x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix.add x3 x2; x3

let ( +@ ) = add

let sub x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix.sub x3 x2; x3

let ( -@ ) = sub

(* TODO: far too slow! need to find a solution! *)
let dot x1 x2 =
  let open Gsl.Blas in
  let x3 = empty (row_num x1) (col_num x2) in
  gemm ~ta:NoTrans ~tb:NoTrans ~alpha:1. ~beta:0. ~a:x1 ~b:x2 ~c:x3; x3

let ( $@ ) = dot

let mul x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix.mul_elements x3 x2; x3

let ( *@ ) = mul

let div x1 x2 =
  let x3 = clone x1 in
  Gsl.Matrix.div_elements x3 x2; x3

let ( /@ ) = div

let power x c = map (fun y -> y ** c) x

let ( **@ ) = power

let abs x = map abs_float x

let neg x =
  let y = clone x in
  Gsl.Matrix.scale y (-1.); y

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

let average x = (sum x) /. (float_of_int (numel x))

let average_cols x =
  let m, n = shape x in
  let y = create n 1 (1. /. (float_of_int n)) in
  dot x y

let average_rows x =
  let m, n = shape x in
  let y = create 1 m (1. /. (float_of_int m)) in
  dot y x

let stderr = None

let stderr_cols = None

let stderr_rows = None

let is_equal x1 x2 =
  let open Owl_matrix_foreign in
  let x1 = mat_to_matptr x1 in
  let x2 = mat_to_matptr x2 in
  (gsl_matrix_equal x1 x2) = 1

let ( =@ ) = is_equal

let is_unequal x1 x2 = not (is_equal x1 x2)

let ( <>@ ) = is_unequal

let is_greater x1 x2 =
  let open Owl_matrix_foreign in
  let x3 = sub x1 x2 in
  let x3 = mat_to_matptr x3 in
  (gsl_matrix_ispos x3) = 1

let ( >@ ) = is_greater

let is_smaller x1 x2 = is_greater x2 x1

let ( <@ ) = is_smaller

let equal_or_greater x1 x2 =
  let open Owl_matrix_foreign in
  let x3 = sub x1 x2 in
  let x3 = mat_to_matptr x3 in
  (gsl_matrix_isnonneg x3) = 1

let ( >=@ ) = equal_or_greater

let equal_or_smaller x1 x2 = equal_or_greater x2 x1

let ( <=@ ) = equal_or_smaller

let is_zero x =
  let open Owl_matrix_foreign in
  let x = mat_to_matptr x in
  (gsl_matrix_isnull x) = 1

let is_positive x =
  let open Owl_matrix_foreign in
  let x = mat_to_matptr x in
  (gsl_matrix_ispos x) = 1

let is_negative x =
  let open Owl_matrix_foreign in
  let x = mat_to_matptr x in
  (gsl_matrix_isneg x) = 1

let is_nonnegative x =
  let open Owl_matrix_foreign in
  let x = mat_to_matptr x in
  (gsl_matrix_isnonneg x) = 1

let min x =
  let open Owl_matrix_foreign in
  let open Ctypes in
  let x = mat_to_matptr x in
  let i = allocate int 0 in
  let j = allocate int 0 in
  let r = gsl_matrix_min x in
  let _ = gsl_matrix_min_index x i j in
  r, !@i, !@j

let min_cols x =
  mapi_cols (fun j v ->
    let r, i, _ = min v in r, i, j
  ) x

let min_rows x =
  mapi_rows (fun i v ->
    let r, _, j = min v in r, i, j
  ) x

let max x =
  let open Owl_matrix_foreign in
  let open Ctypes in
  let x = mat_to_matptr x in
  let i = allocate int 0 in
  let j = allocate int 0 in
  let r = gsl_matrix_max x in
  let _ = gsl_matrix_max_index x i j in
  r, !@i, !@j

let max_cols x =
  mapi_cols (fun j v ->
    let r, i, _ = max v in r, i, j
  ) x

let max_rows x =
  mapi_rows (fun i v ->
    let r, _, j = max v in r, i, j
  ) x

let minmax x =
  let xmin, irow, icol = min x in
  let xmax, arow, acol = max x in
  xmin, xmax, irow, icol, arow, acol

let ( +$ ) x a =
  let y = clone x in
  Gsl.Matrix.add_constant y a; y

let ( $+ ) a x = ( +$ ) x a

let ( -$ ) x a = ( +$ ) x (-1. *. a)

let ( $- ) a x = ( -$ ) x a

let ( *$ ) x a =
  let y = clone x in
  Gsl.Matrix.scale y a; y

let ( $* ) a x = ( *$ ) x a

let ( /$ ) x a = ( *$ ) x (1. /. a)

let ( $/ ) a x = ( /$ ) x a

let add_scalar = ( +$ )

let sub_scalar = ( -$ )

let mul_scalar = ( *$ )

let div_scalar = ( /$ )

(* advanced matrix methematical operations *)

let log x = map log x

let log10 x = map log10 x

let exp x = map exp x

let sigmoid x = map (fun y -> 1. /. (1. +. (Pervasives.exp (-1. *. y)))) x

let diag x =
  let m = Pervasives.min (row_num x) (col_num x) in
  let y = empty 1 m in
  for i = 0 to m - 1 do y.{0,i} <- x.{i,i} done; y

let trace x = sum (diag x)

let add_diag x a =
  let m = Pervasives.min (row_num x) (col_num x) in
  let y = clone x in
  for i = 0 to m - 1 do y.{i,i} <- x.{i,i} +. a done; y


(* formatted input / output operations *)

let to_array x = Gsl.Matrix.to_array x

let to_arrays x = Gsl.Matrix.to_arrays x

let of_array x m n = Gsl.Matrix.of_array x m n

let of_arrays x = Gsl.Matrix.of_arrays x

let save_txt x f =
  let h = open_out f in
  iter_rows (fun y ->  (* TODO: 64-bit -> 16 digits *)
    iter (fun z -> Printf.fprintf h "%.8f\t" z) y;
    Printf.fprintf h "\n"
  ) x;
  close_out h

let load_txt f =
  let h = open_in f in
  let s = input_line h in
  let n = List.length(Str.split (Str.regexp "\t") s) in
  let m = ref 1 in (* counting lines in the input file *)
  let _ = try while true do ignore(input_line h); m := !m + 1
    done with End_of_file -> () in
  let x = zeros !m n in seek_in h 0;
  for i = 0 to !m - 1 do
    let s = Str.split (Str.regexp "\t") (input_line h) in
    List.iteri (fun j y -> x.{i,j} <- float_of_string y) s
  done;
  close_in h; x

let save x f =
  let s = Marshal.to_string x [] in
  let h = open_out f in
  output_string h s;
  close_out h

let load f =
  let h = open_in f in
  let s = really_input_string h (in_channel_length h) in
  Marshal.from_string s 0

let print x = let open Owl_pretty in
  Format.printf "%a\n" Owl_pretty.pp_fmat x;;

let pp_dsmat x = let open Owl_pretty in
  Format.printf "%a\n" Toplevel.pp_fmat x;;

(* some other uncategorised functions *)

let uniform_int ?(a=0) ?(b=99) m n =
  let x = empty m n in
  iteri (fun i j _ -> x.{i,j} <-
    float_of_int (Owl_stats.Rnd.uniform_int ~a ~b ())
  ) x; x

let uniform ?(scale=1.) m n =
  let x = empty m n in
  iteri (fun i j _ ->
    x.{i,j} <- Owl_stats.Rnd.uniform () *. scale
  ) x; x

let gaussian ?(sigma=1.) m n =
  let x = empty m n in
  iteri (fun i j _ -> x.{i,j} <- Owl_stats.Rnd.gaussian ~sigma ()) x; x

let vector_uniform n = uniform 1 n

let semidef n =
  let x = uniform n n in
  dot (transpose x) x

let draw_rows ?(replacement=true) x c =
  let a = Array.init (row_num x - 1) (fun i -> i) in
  let l = match replacement with
    | true  -> Owl_stats.sample a c
    | false -> Owl_stats.choose a c
  in rows x l, l

let draw_cols ?(replacement=true) x c =
  let a = Array.init (col_num x - 1) (fun i -> i) in
  let l = match replacement with
    | true  -> Owl_stats.sample a c
    | false -> Owl_stats.choose a c
  in cols x l, l

let shuffle_rows x =
  let y = clone x in
  let m, n = shape x in
  for i = 0 to m - 1 do
    swap_rows y i (Owl_stats.Rnd.uniform_int ~a:0 ~b:(m-1) ())
  done; y

let shuffle_cols x =
  let y = clone x in
  let m, n = shape x in
  for i = 0 to n - 1 do
    swap_cols y i (Owl_stats.Rnd.uniform_int ~a:0 ~b:(n-1) ())
  done; y

let shuffle x = x |> shuffle_rows |> shuffle_cols

let reshape m n x = of_array (to_array x) m n

let meshgrid xa xb ya yb xn yn =
  let u = linspace xa xb xn in
  let v = linspace ya yb yn in
  let x = map_by_row (fun _ -> u) (empty yn xn) in
  let y = map_by_row (fun _ -> v) (empty xn yn) in
  x, transpose y

let meshup x y =
  let xn = numel x in
  let yn = numel y in
  let x = map_by_row (fun _ -> x) (empty yn xn) in
  let y = map_by_row (fun _ -> y) (empty xn yn) in
  x, transpose y

let ( @@ ) f x = map f x  (* TODO: experimental *)

(* TODO: use this to replace col function, faster *)
let gsl_col x i =
  let open Owl_matrix_foreign in
  let y = allocate_col_vecptr (row_num x) in
  let _ = gsl_matrix_get_col y.vptr (mat_to_matptr x) i in
  y.vdata
