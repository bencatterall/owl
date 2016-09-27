(** [ Sparse Matrix ]  *)

type spmat

(** sparse matrix creation and basic functions *)

val empty : int -> int -> spmat

val shape : spmat -> int * int

val row_num : spmat -> int

val col_num : spmat -> int

val nnz : spmat -> int

val density : spmat -> float

val eye : int -> spmat

val uniform : ?scale:float -> int -> int -> spmat

val uniform_int : ?a:int -> ?b:int -> int -> int -> spmat

(** matrix manipulations *)

val set : spmat -> int -> int -> float -> unit

val get : spmat -> int -> int -> float

val reset : spmat -> unit

val clone : spmat -> spmat

val transpose : spmat -> spmat

(** matrix interation functions *)

val row : spmat -> int -> spmat

val col : spmat -> int -> spmat

val rows : spmat -> int array -> spmat

val cols : spmat -> int array -> spmat

val iteri : (int -> int -> float -> 'a) -> spmat -> unit

val iter : (float -> 'a) -> spmat -> unit

val mapi : (int -> int -> float -> float) -> spmat -> spmat

val map : (float -> float) -> spmat -> spmat

val fold : ('a -> float -> 'a) -> 'a -> spmat -> 'a

val filteri : (int -> int -> float -> bool) -> spmat -> (int * int) array

val filter : (float -> bool) -> spmat -> (int * int) array

val iteri_rows : (int -> spmat -> unit) -> spmat -> unit

val iter_rows : (spmat -> unit) -> spmat -> unit

val iteri_cols : (int -> spmat -> unit) -> spmat -> unit

val iter_cols : (spmat -> unit) -> spmat -> unit

val mapi_rows : (int -> spmat -> spmat) -> spmat -> spmat

val mapi_cols : (int -> spmat -> spmat) -> spmat -> spmat

val fold_rows : ('a -> spmat -> 'a) -> 'a -> spmat -> 'a

val fold_cols : ('a -> spmat -> 'a) -> 'a -> spmat -> 'a

val iteri_nz : (int -> int -> float -> 'a) -> spmat -> unit

val iter_nz : (float -> 'a) -> spmat -> unit

val mapi_nz : (int -> int -> float -> float) -> spmat -> spmat

val map_nz : (float -> float) -> spmat -> spmat

val fold_nz : ('a -> float -> 'a) -> 'a -> spmat -> 'a

val filteri_nz : (int -> int -> float -> bool) -> spmat -> (int * int) array

val filter_nz : (float -> bool) -> spmat -> (int * int) array

val iter_rows_nz : (spmat -> unit) -> spmat -> unit

val iteri_rows_nz : (int -> spmat -> unit) -> spmat -> unit

val iter_cols_nz : (spmat -> unit) -> spmat -> unit

val iteri_cols_nz : (int -> spmat -> unit) -> spmat -> unit

(** TODO: val mapi_rows_nz : *)

(** TODO: val mapi_cols_nz : *)

val fold_rows_nz : ('a -> spmat -> 'a) -> 'a -> spmat -> 'a

val fold_cols_nz : ('a -> spmat -> 'a) -> 'a -> spmat -> 'a

val exists : (float -> bool) -> spmat -> bool

val not_exists : (float -> bool) -> spmat -> bool

val for_all : (float -> bool) -> spmat -> bool

val exists_nz : (float -> bool) -> spmat -> bool

val not_exists_nz : (float -> bool) -> spmat -> bool

val for_all_nz :  (float -> bool) -> spmat -> bool

(** matrix mathematical operations *)

val mul_scalar : spmat -> float -> spmat

val div_scalar : spmat -> float -> spmat

val add : spmat -> spmat -> spmat

val sub : spmat -> spmat -> spmat

val dot : spmat -> spmat -> spmat

val mul : spmat -> spmat -> spmat

val div : spmat -> spmat -> spmat

val abs : spmat -> spmat

val neg : spmat -> spmat

val sum : spmat -> float

val average : spmat -> float

val is_zero : spmat -> bool

val is_positive : spmat -> bool

val is_negative : spmat -> bool

val is_nonnegative : spmat -> bool

val min : spmat -> float

val max : spmat -> float

val minmax : spmat -> float * float

(** compare two matrices *)

val is_equal : spmat -> spmat -> bool

val is_unequal : spmat -> spmat -> bool

val is_greater : spmat -> spmat -> bool

val is_smaller : spmat -> spmat -> bool

val equal_or_greater : spmat -> spmat -> bool

val equal_or_smaller : spmat -> spmat -> bool

(** advanced matrix methematical operations *)

val diag : spmat -> spmat

val trace : spmat -> float

(** TODO: val svd : *)

(** TODO: val qr : *)

(** TODO: val lu : *)

(** transform to and from different types *)

val to_dense : spmat -> Dense.dsmat

val of_dense : Dense.dsmat -> spmat

(** some other uncategorised functions *)

val print : spmat -> unit

val pp_spmat : spmat -> unit

val draw_rows : ?replacement:bool -> spmat -> int -> spmat

val draw_cols : ?replacement:bool -> spmat -> int -> spmat

val permutation_matrix : int -> spmat

val shuffle_rows : spmat -> spmat

val shuffle_cols : spmat -> spmat

val shuffle : spmat -> spmat

(** short-hand infix operators *)

val ( +@ ) : spmat -> spmat -> spmat

val ( -@ ) : spmat -> spmat -> spmat

val ( *@ ) : spmat -> spmat -> spmat

val ( /@ ) : spmat -> spmat -> spmat

val ( $@ ) : spmat -> spmat -> spmat

val ( **@ ) : spmat -> float -> spmat

val ( *$ ) : spmat -> float -> spmat

val ( /$ ) : spmat -> float -> spmat

val ( $* ) : float -> spmat -> spmat

val ( $/ ) : float -> spmat -> spmat

val ( =@ ) : spmat -> spmat -> bool

val ( >@ ) : spmat -> spmat -> bool

val ( <@ ) : spmat -> spmat -> bool

val ( <>@ ) : spmat -> spmat -> bool

val ( >=@ ) : spmat -> spmat -> bool

val ( <=@ ) : spmat -> spmat -> bool

(* val ( @@ ) : (float -> float) -> dsmat -> dsmat *)