Subroutine B_CalcGradient(NI, NJ, P, gradP, gradPLast, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
	integer NI, NJ	
	Real, Dimension(0:NI,0:NJ)::P
	Real,Dimension(0:NI,0:NJ,2)::gradP, gradPLast
	Real, Dimension(NI-1,NJ-1):: CellVolume ! Cell Volumes
	Real, Dimension(0:NI,0:NJ,2):: CellCenter
	Real, Dimension(NI,NJ-1,2):: IFaceCenter
	Real, Dimension(NI,NJ-1,2):: IFaceVector
	Real, Dimension(NI-1,NJ,2):: JFaceCenter
	Real, Dimension(NI-1,NJ,2):: JFaceVector
    
	integer i, j, IFace, I_N, J_N
	real pE, pM, linearInterpolation
	Real, Dimension(2):: faceVector, faceCenter, E, EM, gradPE
	!gradP(:, :, :) = 0.0
	gradPLast = gradP
	gradP(:, :, :) = 0.0
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
				
				pE = linearInterpolation(faceCenter(:), CellCenter(i, j, :), P(i, j), CellCenter(I_N, J_N, :), P(I_N, J_N))
				E(1) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), CellCenter(i, j, 1), CellCenter(I_N, J_N, :), CellCenter(I_N, J_N, 1))
				E(2) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), CellCenter(i, j, 2), CellCenter(I_N, J_N, :), CellCenter(I_N, J_N, 2))
				EM = faceCenter - E
				gradPE(1) = linearInterpolation(E(:), CellCenter(i, j, :), gradPLast(i, j, 1), CellCenter(I_N, J_N, :), gradPLast(I_N, J_N, 1))
				gradPE(2) = linearInterpolation(E(:), CellCenter(i, j, :), gradPLast(i, j, 2), CellCenter(I_N, J_N, :), gradPLast(I_N, J_N, 2))
				pM = pE + DOT_PRODUCT(gradPE, EM) !iteration
				
				!d1 = Norm2(CellCenter(i, j, :) - faceCenter(:))
				!d2 = Norm2(CellCenter(I_N, J_N, :) - faceCenter(:))
				!pM = (P(i, j) * d2 + P(I_N, J_N) * d1) / (d1 + d2)
				!pM = linearInterpolation(faceCenter(:), CellCenter(i, j, :), P(i, j), CellCenter(I_N, J_N, :), P(I_N, J_N)) !no iteration
				gradP(i, j, :) = gradP(i, j, :) + pM * faceVector(:)
			enddo
			gradP(i, j, :) = gradP(i, j, :) / CellVolume(i, j)
		enddo
	enddo
End Subroutine 
