MODULE R_interfaces
  IMPLICIT NONE
  ! Declare R's built-in print functions as external
  EXTERNAL :: intpr, dblepr
END MODULE R_interfaces