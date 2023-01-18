Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,gradP, gradPError, V, divV, divVError, laplacianP, laplacianPError, T, TExact, TError)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P, divV, divVError, laplacianP, laplacianPError, T, TError, TExact
  Real,Dimension(0:NI,0:NJ, 2)::gradP, gradPError, V

  Write(IO,*) 'VARIABLES = "X", "Y", "P", "gradP_x", "gradP_y", "V_x", "V_y", "divV", "T", "TExact", "TError"' 
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') gradP(1:NI-1,1:NJ-1, 1)
  Write(IO,'(100F14.7)') gradP(1:NI-1,1:NJ-1, 2)
  !Write(IO,'(100F14.7)') gradPError(1:NI-1,1:NJ-1, 1)
  !Write(IO,'(100F14.7)') gradPError(1:NI-1,1:NJ-1, 2)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1, 1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1, 2)
  Write(IO,'(100F14.7)') divV(1:NI-1,1:NJ-1)
  !Write(IO,'(100F14.7)') divVError(1:NI-1,1:NJ-1)
  !Write(IO,'(100F14.7)') laplacianP(1:NI-1,1:NJ-1)
  !Write(IO,'(100F14.7)') laplacianPError(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') TExact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') TError(1:NI-1,1:NJ-1)
End Subroutine 
