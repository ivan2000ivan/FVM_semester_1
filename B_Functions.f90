Function Pressure(X,Y)
  Pressure = x**2 + y**2
End Function

Subroutine Velocity(V, X, Y)
	Real, Dimension(2) :: V
	real x, y
	V(1) = 1.0 * x
	V(2) = 1.0 * y
End Subroutine

Subroutine CalcGradPExact(NI, NJ, gradPExact, CellCenter)
	Real, Dimension(0:NI,0:NJ,2) :: gradPExact
	Real, Dimension(0:NI,0:NJ,2) :: CellCenter
	gradPExact(:, :, 1) = 2.0 * CellCenter(:, :, 1) 
	gradPExact(:, :, 2) = 1.0 * CellCenter(:, :, 2) 
end subroutine

Subroutine CalcDivVExact(NI, NJ, divVExact, CellCenter)
	Real, Dimension(0:NI,0:NJ) :: divVExact
	Real, Dimension(0:NI,0:NJ,2) :: CellCenter
	divVExact = 3.0 * CellCenter(:, :, 1) + 3.0 * CellCenter(:, :, 2)
end subroutine

subroutine CalcLaplacianPExact(NI, NJ, laplacianPExact, CellCenter)
	Real, Dimension(0:NI,0:NJ) :: laplacianPExact
	Real, Dimension(0:NI,0:NJ,2) :: CellCenter
	laplacianPExact = 4
end subroutine

function linearInterpolation(pointTarget, point1, value1, point2, value2)
	real, dimension(2) :: pointTarget, point1, point2
	real value1, value2, d1, d2, linearInterpolation
	
	d1 = Norm2(point1(:) - pointTarget(:))
	d2 = Norm2(point2(:) - pointTarget(:))
	
	linearInterpolation = (value1 * d2 + value2 * d1) / (d1 + d2)
end function
