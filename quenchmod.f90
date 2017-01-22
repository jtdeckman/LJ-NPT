MODULE QUENCHMOD
  USE ISO_C_BINDING
  INTERFACE
    FUNCTION QUENCHF(Q, N, MI, NIT) BIND(C,name="c_quench")
      USE ISO_C_BINDING
      INTEGER(c_int), value :: N
      INTEGER(c_int), value :: MI
      INTEGER(c_int) :: NIT(*)
      real(c_double) :: Q(*)
      real(c_double) :: QUENCHF
    END FUNCTION QUENCHF
 END INTERFACE
END MODULE QUENCHMOD
