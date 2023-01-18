Subroutine solveT(NI, NJ, V, T, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
	integer NI, NJ, scheme	
	Real, Dimension(0:NI,0:NJ,2)::V, gradT, gradTLast
	Real, Dimension(0:NI,0:NJ)::T, TLast, divTV, lapT, res
	Real, Dimension(NI-1,NJ-1):: CellVolume ! Cell Volumes
	Real, Dimension(0:NI,0:NJ,2):: CellCenter
	Real, Dimension(NI,NJ-1,2):: IFaceCenter
	Real, Dimension(NI,NJ-1,2):: IFaceVector
	Real, Dimension(NI-1,NJ,2):: JFaceCenter
	Real, Dimension(NI-1,NJ,2):: JFaceVector
	
	Real, Dimension(0:NI,0:NJ)::tau
	Real, Dimension(2)::gradTF, VF, faceVector, faceCenter
	integer m, mMax, i, j, IFace, k
	real sum_f, maxRes, Re, Pr, maxResT, TF, v_sigma, eps, dtConvX, dtConvY, dtDiffX, dtDiffY, CFL, VNM, dyi, dxi, Uxi, Uyi
	sum_f = 0
	tau = 0.001
	mMax = 500
	Re = 100
	Pr = 0.7
	T = 0
	divTV = 0
	lapT = 0
	res = 0
	eps = 0.1
	maxResT = 1E4
	CFL = 0.3
	VNM = 0.3
	m = 0
	
	do i = 1, NI-1
		do j = 1, NJ - 1
			dxi = Norm2(IFaceCenter(i + 1, j, :) - IFaceCenter(i, j, :))
			dyi = Norm2(IFaceCenter(i, j + 1, :) - IFaceCenter(i, j, :))
			Ui = Norm2(V(i, j, :))
			dtConvX = CFL * dxi / Ui
			dtDiffX = VNM * dxi * dxi / 2.0 * Re * Pr
			dtConvY = CFL * dyi / Ui
			dtDiffY = VNM * dyi * dyi / 2.0 * Re * Pr
			tau(i, j) = min(1.0 / (1.0 / dtConvX + 1.0 / dtDiffX), 1.0 / (1.0 / dtConvY + 1.0 / dtDiffY));
		enddo
	enddo
	
	print*, 'tau min/max:', MINVAL(tau), MAXVAL(tau)
	
	do while (maxResT > 0.000001)
	
		T(0, :) = 0.0
		T(NI, :) = 1.0
		T(:, 0) = T(:, 1)
		T(:, NJ) = T(:, NJ - 1)
	
		gradT = 0.0 
		gradTLast = 0.0  
		iterationNum = 100
		eps = 1.0E-10
		maxRes = 1E6
		k = 1;
		!'--------Calc Grad--------'
		do while ((maxRes > eps) .AND. (k < iterationNum))
			Call B_CalcGradient(NI, NJ, T, gradT, gradTLast, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
			maxRes = MAXVAL(abs(gradT - gradTLast))
			!WRITE(*,*) 'Iteration num = ', k, ', maxRes = ', maxRes
			k = k + 1
		enddo
		!'--------End Calc Grad--------'
		
		lapT = 0
		divTV = 0
		call B_CalcDiv(3, NI, NJ, V, T, gradT, divTV, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
		call B_CalcLaplacian(NI, NJ, T, gradT, lapT, CellVolume, CellCenter, IFaceVector, JFaceVector, IFaceCenter, JFaceCenter)
		res = (divTV - lapT / (Re * Pr))
		T = T - tau * res
		
		maxResT = MAXVAL(abs(res))
		m = m + 1
		WRITE(*,*) 'Iteration num = ', m, ', maxRes = ', maxResT
		
	enddo
	
end Subroutine