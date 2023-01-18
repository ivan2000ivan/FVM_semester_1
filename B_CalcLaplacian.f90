Subroutine B_CalcLaplacian(NI, NJ, P, gradP, laplacianP, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
	integer NI, NJ
	Real, Dimension(0:NI,0:NJ)::P
	Real, Dimension(0:NI,0:NJ)::laplacianP
	Real, Dimension(NI-1,NJ-1):: CellVolume ! Cell Volumes
	Real, Dimension(0:NI,0:NJ,2):: CellCenter, gradP
	Real, Dimension(NI,NJ-1,2):: IFaceCenter
	Real, Dimension(NI,NJ-1,2):: IFaceVector
	Real, Dimension(NI-1,NJ,2):: JFaceCenter
	Real, Dimension(NI-1,NJ,2):: JFaceVector
    
	integer i, j, IFace, I_N, J_N
	real linearInterpolation, dp_dn
	Real, Dimension(2):: faceVector, faceCenter, i_ksi, gradPf
	
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
				
				if ((i .EQ. 1) .AND. (IFace .EQ. 3)) then 
					dp_dn = (5.0/3.0) * (P(I_N, J_N) - P(i, j)) / Norm2(CellCenter(i, j, :) - CellCenter(I_N, J_N, :)) - (2.0/3.0) * DOT_PRODUCT(gradP(i, j, :), faceVector(:) / Norm2(faceVector(:))) 
				else if ((i .EQ. NI-1) .AND. (IFace .EQ. 4)) then 
					dp_dn = (5.0/3.0) * (P(I_N, J_N) - P(i, j)) / Norm2(CellCenter(i, j, :) - CellCenter(I_N, J_N, :)) - (2.0/3.0) * DOT_PRODUCT(gradP(i, j, :), faceVector(:) / Norm2(faceVector(:)))
				else if ((j .EQ. 1) .AND. (IFace .EQ. 1)) then 
					dp_dn = (5.0/3.0) * (P(I_N, J_N) - P(i, j)) / Norm2(CellCenter(i, j, :) - CellCenter(I_N, J_N, :)) - (2.0/3.0) * DOT_PRODUCT(gradP(i, j, :), faceVector(:) / Norm2(faceVector(:)))
				else if ((j .EQ. NJ-1) .AND. (IFace .EQ. 2)) then 
					dp_dn = (5.0/3.0) * (P(I_N, J_N) - P(i, j)) / Norm2(CellCenter(i, j, :) - CellCenter(I_N, J_N, :)) - (2.0/3.0) * DOT_PRODUCT(gradP(i, j, :), faceVector(:) / Norm2(faceVector(:)))
				else 
					dp_dn = (P(I_N, J_N) - P(i, j)) / Norm2(CellCenter(i, j, :) - CellCenter(I_N, J_N, :))
					i_ksi = CellCenter(I_N, J_N, :) - CellCenter(i, j, :)
				    i_ksi = i_ksi / Norm2(i_ksi)
				    gradPf(1) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), gradP(i, j, 1), CellCenter(I_N, J_N, :), gradP(I_N, J_N, 1))
				    gradPf(2) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), gradP(i, j, 2), CellCenter(I_N, J_N, :), gradP(I_N, J_N, 2))
				    dp_dn = dp_dn + DOT_PRODUCT((faceVector / Norm2(faceVector(:))  - i_ksi), gradPf)
				endif
				
				
				
				!i_ksi = CellCenter(I_N, J_N, :) - CellCenter(i, j, :)
				!i_ksi = i_ksi / Norm2(i_ksi)
				!gradPf(1) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), gradP(i, j, 1), CellCenter(I_N, J_N, :), gradP(I_N, J_N, 1))
				!gradPf(2) = linearInterpolation(faceCenter(:), CellCenter(i, j, :), gradP(i, j, 2), CellCenter(I_N, J_N, :), gradP(I_N, J_N, 2))
				!dp_dn = dp_dn + DOT_PRODUCT((faceVector / Norm2(faceVector(:))  - i_ksi), gradPf)
				
				laplacianP(i, j) = laplacianP(i, j) + dp_dn * Norm2(faceVector(:))
			enddo
			laplacianP(i, j) = laplacianP(i, j) / CellVolume(i, j)
		enddo
	enddo
End Subroutine 
