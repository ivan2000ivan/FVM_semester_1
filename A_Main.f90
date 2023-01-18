Program Main

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' ! names of input and output files
  character MeshFile*30 , ctmp       ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume, divV, divVError, divVExact, laplacianP, laplacianPExact, laplacianPError, T, TExact, TError! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector, gradP, gradPLast, gradPExact, gradPError, &
									  V! vector arrays
  integer iterationNum, k, scheme
  real eps, maxRes, tmp

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  !Read (IO, *) scheme  ! div(pV) scheme
  READ(IO,*) scheme, MeshFile  ! read name of file with computational mesh
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(T(0:NI,0:NJ))   ! Pressure
  allocate(TExact(0:NI,0:NJ)) 
  allocate(gradP(0:NI,0:NJ, 2))   ! Pressure grad
  allocate(gradPLast(0:NI,0:NJ, 2))   ! Pressure grad
  allocate(gradPExact(0:NI,0:NJ, 2))   ! Pressure grad exact
  allocate(gradPError(0:NI,0:NJ, 2))   ! Pressure grad error
  allocate(V(0:NI,0:NJ, 2))   ! div V
  allocate(divV(0:NI,0:NJ))   ! div V error
  allocate(divVExact(0:NI,0:NJ))   ! Pressure grad exact
  allocate(divVError(0:NI,0:NJ))   ! Pressure grad error
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces
  allocate(laplacianP(0:NI,0:NJ))  !laplacian P
  allocate(laplacianPError(0:NI,0:NJ))  !laplacian P error
  allocate(laplacianPExact(0:NI,0:NJ))  !laplacian P exact

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),tmp,I=1,NI),J=1,NJ)
  CLOSE(IO)
  
!==== READ SOLUTION ===
  write(*, *) 'Read solution'
  open(io, file = 'caverna.plt')
  read(io, *) ctmp
  read(io, *) ctmp
  read(io, *) ((tmp, tmp, V(i, j, 1), V(i, j, 2), tmp, P(i, j), TExact(i,j), tmp, tmp, i = 0, NI), j = 0, NJ)
  close(io)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      !P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
	  !call Velocity(V(I, J, :), CellCenter(I,J,1),CellCenter(i,j,2))
    ENDDO
  ENDDO

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives' 
  gradP = 0.0 
  gradPLast = 0.0  
  iterationNum = 100
  eps = 1.0E-10
  maxRes = 1E6
  k = 1;
  do while ((maxRes > eps) .AND. (k < iterationNum))
	Call B_CalcGradient(NI, NJ, P, gradP, gradPLast, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
	maxRes = MAXVAL(abs(gradP - gradPLast))
	WRITE(*,*) 'Iteration num = ', k, ', maxRes = ', maxRes
	k = k + 1
  enddo
  Call CalcGradPExact(NI, NJ, gradPExact, CellCenter)
  gradPError = abs(gradP - gradPExact) / gradPExact
  
!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate div, scheme:', scheme 
  divV = 0.0
  call B_CalcDiv(scheme, NI, NJ, V, P, gradP, divV, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
  call CalcDivVExact(NI, NJ, divVExact, CellCenter)
  divVError = abs(divV - divVExact) / divVExact
  
!=== CALCULATE LAPLACIAN ===
  WRITE(*,*) 'Calculate laplacian'
  laplacianP = 0
  call B_CalcLaplacian(NI, NJ, P, gradP, laplacianP, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
  call CalcLaplacianPExact(NI, NJ, laplacianPExact, CellCenter)
  laplacianPError = abs (laplacianP - laplacianPExact) / laplacianPExact
 
!=== CALCULATE T === 
  WRITE(*,*) 'Calculate T'
  call solveT(NI, NJ, V, T, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
  TError = abs(T - TExact) / TExact

!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,gradP, gradPError, V, divV, divVError, laplacianP, laplacianPError, T, TExact,TError)
  Close(IO)
  
  deallocate(X, Y, P, gradP, CellVolume, CellCenter, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector, gradPExact, gradPError, gradPLast, &
			 V, divV, divVError, divVExact, laplacianP, laplacianPError, laplacianPExact)
END PROGRAM Main  
