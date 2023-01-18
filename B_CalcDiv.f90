Subroutine B_CalcDiv(scheme, NI, NJ, V, P, gradP, divV, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
	integer NI, NJ, scheme	
	Real, Dimension(0:NI,0:NJ,2)::V
	Real, Dimension(0:NI,0:NJ)::P
	Real, Dimension(0:NI,0:NJ)::divV
	Real, Dimension(NI-1,NJ-1):: CellVolume ! Cell Volumes
	Real, Dimension(0:NI,0:NJ,2):: CellCenter, gradP
	Real, Dimension(NI,NJ-1,2):: IFaceCenter
	Real, Dimension(NI,NJ-1,2):: IFaceVector
	Real, Dimension(NI-1,NJ,2):: JFaceCenter
	Real, Dimension(NI-1,NJ,2):: JFaceVector
    
	integer i, j, IFace, I_N, J_N
	real linearInterpolation, pf, v_sigma
	Real, Dimension(2):: faceVector, faceCenter, vM
	
	do j = 1, NJ - 1
		do i = 1, NI - 1
			
			do IFace = 1, 4
				if (IFace .EQ. 1) then 
					I_N = i
					J_N = j - 1
					faceVector = -JFaceVector(i, j, :)
					faceCenter = JFaceCenter(i, j, :)
				else if (IFace .EQ. 2) then
					I_N = i
					J_N = j + 1
					faceVector = JFaceVector(i, j + 1, :)
					faceCenter = JFaceCenter(i, j + 1, :)
				else if (IFace .EQ. 3) then
					I_N = i - 1
					J_N = j
					faceVector = -IFaceVector(i, j, :)
					faceCenter = IFaceCenter(i, j, :)
				else if (IFace .EQ. 4) then
					I_N = i + 1
					J_N = j
					faceVector = IFaceVector(i + 1, j, :)
					faceCenter = IFaceCenter(i + 1, j, :)
				endif
				
				vM(1) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), V(i, j, 1), CellCenter(I_N, J_N, :), V(I_N, J_N, 1))
				vM(2) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), V(i, j, 2), CellCenter(I_N, J_N, :), V(I_N, J_N, 2))
				
				v_sigma = DOT_PRODUCT(vM, faceVector)
				
				if (scheme .EQ. 1) then
					pf = linearInterpolation(faceCenter(:), CellCenter(i, j, :), P(i, j), CellCenter(I_N, J_N, :), P(I_N, J_N))
				else if (scheme .EQ. 2) then
					if (v_sigma >= 0) then
						pf = p(i, j)
					else 
						pf = p(I_N, J_N) 
					endif
				else if (scheme .EQ. 3) then 
					if (v_sigma >= 0) then
						pf = p(i, j) + DOT_PRODUCT(gradP(i, j, :), faceCenter - CellCenter(i, j, :))
					else 
						pf = p(I_N, J_N) + DOT_PRODUCT(gradP(I_N, J_N, :), faceCenter - CellCenter(I_N, J_N, :))
					endif
				endif
				
				divV(i, j) = divV(i, j) + pf * v_sigma
			enddo
			divV(i, j) = divV(i, j) / CellVolume(i, j)
		enddo
	enddo
End Subroutine 
