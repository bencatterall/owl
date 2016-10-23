

open Bigarray


module type MatrixSig = sig

  type elt
  type prc
  type mat = (elt, prc, c_layout) Array2.t

  val const_0 : elt

  val const_1 : elt

  val empty : int -> int -> mat

  val create : int -> int -> elt -> mat

end


module Common (MatrixImpl : MatrixSig) = struct

  open MatrixImpl

  type area = { a : int; b : int; c : int; d : int }

  let shape x = (Array2.dim1 x, Array2.dim2 x)

  let row_num x = fst (shape x)

  let col_num x = snd (shape x)

  let numel x = (row_num x) * (col_num x)

  let zeros m n = create m n const_0

  let ones m n = create m n const_1

  let eye n =
    let x = zeros n n in
    for i = 0 to n - 1 do
      x.{i,i} <- const_1
    done; x

  let vector n = empty 1 n

  let vector_ones n = ones 1 n

  let vector_zeros n = zeros 1 n

  let same_shape x1 x2 = shape x1 = shape x2

  let area a b c d = { a = a; b = b; c = c; d = d }

  let area_of x =
    let m, n = shape x in
    { a = 0; b = 0; c = m - 1; d = n - 1 }

  let area_of_row x i = area i 0 i (col_num x - 1)

  let area_of_col x i = area 0 i (row_num x - 1) i

  let equal_area r1 r2 =
    ((r1.c-r1.a = r2.c-r2.a) && (r1.d-r1.b = r2.d-r2.b))

  let same_area r1 r2 = r1 = r2

  let set = Array2.unsafe_set

  let get = Array2.unsafe_get

  let row x i =
    let y = Array2.slice_left x i in
    reshape_2 (genarray_of_array1 y) 1 (col_num x)

  let col x j =
    let m, n = shape x in
    let y = empty m 1 in
    for i = 0 to m - 1 do
      set y i 0 (get x i j)
    done; y

  let copy_area_to x1 r1 x2 r2 =
    if not (equal_area r1 r2) then
      failwith "Error: area mismatch"
    else
      for i = 0 to r1.c - r1.a do
        for j = 0 to r1.d - r1.b do
          set x2 (r2.a + i) (r2.b + j)
          (get x1 (r1.a + i) (r1.b + j))
        done
      done

  let copy_to x1 x2 = Array2.blit x1 x2

  let ( >> ) = copy_to

  let ( << ) x1 x2 = copy_to x2 x1

  let clone_area x r =
    let y = empty (r.c - r.a + 1) (r.d - r.b + 1) in
    copy_area_to x r y (area_of y)

end
