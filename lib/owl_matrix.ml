

open Bigarray


module type MatrixSig = sig
  type mat
  type elt

  val create : int -> int -> elt -> mat

end


module Matrix (MatrixImpl : MatrixSig) = struct
  type mat = MatrixImpl.mat
  type elt = MatrixImpl.elt

  let shape x = (Array2.dim1 x, Array2.dim2 x)

end
