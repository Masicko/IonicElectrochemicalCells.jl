function half_electrolyte_1D(h,dmin) 
	X = ExtendableGrids.geomspace(0,h , dmin, 2e-1*h)
	grid = ExtendableGrids.simplexgrid(X)
    btol = 1e-2*dmin
	ExtendableGrids.cellmask!(grid, [0.0], [h], Ω_YSZ, tol = btol)
	ExtendableGrids.bfacemask!(grid, [0.0], [0.0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [h], [h], 2, tol=btol)
    grid
end

function full_electrolyte_1D(h,dmin) 
	Xleft = ExtendableGrids.geomspace(0,h , dmin, 2e-1*h)
	Xright= ExtendableGrids.geomspace(h,2*h, 2e-1*h, dmin)
	X = ExtendableGrids.glue(Xleft,Xright)
	grid = ExtendableGrids.simplexgrid(X)
    btol = 1e-2*dmin
	ExtendableGrids.cellmask!(grid, [0.0], [2*h], Ω_YSZ, tol = btol)
	ExtendableGrids.bfacemask!(grid, [0.0], [0.0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [2*h], [2*h], 2, tol=btol)
    grid
end

# merge cell1D and simplefullcell1D

function cell1D(hEL, hElec, hER; dmin) 
	x0 = 0.0
	x1 = hEL
	x2 = hEL + 0.5*hElec
	x3 = hEL + hElec
	x4 = hEL + hElec + hER
	eldaLeft = ExtendableGrids.geomspace(x0, x1, 1e-1*hEL, dmin)
	elecLeft = ExtendableGrids.geomspace(x1, x2, dmin, 1e-1*hElec)
	elecRight= ExtendableGrids.geomspace(x2, x3, 1e-1*hElec, dmin)
	eldaRight= ExtendableGrids.geomspace(x3, x4, dmin, 1e-1*hER)
	X1 = ExtendableGrids.glue(eldaLeft,elecLeft)
	X2 = ExtendableGrids.glue(X1,elecRight)
	X3 = ExtendableGrids.glue(X2,eldaRight)
	grid = ExtendableGrids.simplexgrid(X3)
    btol = 1e-2*dmin
	ExtendableGrids.cellmask!(grid, [x0], [x1], Ω_Aul, tol=btol)
	ExtendableGrids.cellmask!(grid, [x1], [x3], Ω_YSZ, tol=btol)
	ExtendableGrids.cellmask!(grid, [x3], [x4], Ω_Aur, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x0], [x0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x1], [x1], 3, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x3], [x3], 4, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x4], [x4], 2, tol=btol)
    grid
end

function simplefullcell1D(hElec, hElda; N=3) 
    btol = min(1e-2*N^(-1),1e-10)
	x0 = 0.0
	x1 = hElda
	x3 = hElda + hElec
	x4 = 2*hElda + hElec
	eldaLeft = ExtendableGrids.LinRange(x0, x1, N)
	elec= ExtendableGrids.LinRange(x1, x3, N)
	eldaRight= ExtendableGrids.LinRange(x3, x4, N)
	X1 = ExtendableGrids.glue(eldaLeft,elec)
	X3 = ExtendableGrids.glue(X1,eldaRight)
	grid = ExtendableGrids.simplexgrid(X3)
	ExtendableGrids.cellmask!(grid, [x0], [x1], 1, tol=btol)
	ExtendableGrids.cellmask!(grid, [x1], [x3], 2, tol=btol)
	ExtendableGrids.cellmask!(grid, [x3], [x4], 3, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x0], [x0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x1], [x1], 2, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x3], [x3], 3, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x4], [x4], 4, tol=btol)
    grid
end

function half_cell1D(hEL, hElec, hER; dmin) 
	x0 = 0.0
	x1 = hEL
	x2 = hEL + 0.5*hElec
	eldaLeft = ExtendableGrids.geomspace(x0, x1, 1e-1*hEL, dmin)
	elecLeft = ExtendableGrids.geomspace(x1, x2, dmin, 1e-1*hElec)
	X1 = ExtendableGrids.glue(eldaLeft,elecLeft)
	grid = ExtendableGrids.simplexgrid(X1)
    btol = 1e-2*dmin
	ExtendableGrids.cellmask!(grid, [x0], [x1], Ω_Aul, tol=btol)
	ExtendableGrids.cellmask!(grid, [x1], [x2], Ω_YSZ, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x0], [x0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x1], [x1], 2, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x1], [x1], 3, tol=btol)
    grid
end

function doublehalf_cell1D(hEL, hElec, hER; dmin) 
	x0 = 0.0
	x1 = hEL
	x2 = hEL + 0.5*hElec
	x3 = hEL + hElec
	x4 = hEL + hElec + hER
	eldaLeft = ExtendableGrids.geomspace(x0, x1, 1e-1*hEL, dmin)
	elecLeft = ExtendableGrids.geomspace(x1, x2, dmin, 1e-1*hElec)
	elecRight= ExtendableGrids.geomspace(x2, x3, 1e-1*hElec, dmin)
	eldaRight= ExtendableGrids.geomspace(x3, x4, dmin, 1e-1*hER)
	X1 = ExtendableGrids.glue(eldaLeft,elecLeft)
	X2 = ExtendableGrids.glue(X1,elecRight)
	X3 = ExtendableGrids.glue(X2,eldaRight)
	grid = ExtendableGrids.simplexgrid(X3)
    btol = 1e-2*dmin
	ExtendableGrids.cellmask!(grid, [x0], [x1], Ω_Aul, tol=btol)
	ExtendableGrids.cellmask!(grid, [x1], [x2], e_Ω_YSZl, tol=btol)
	ExtendableGrids.cellmask!(grid, [x3], [x4], Ω_Aur, tol=btol)
	ExtendableGrids.cellmask!(grid, [x2], [x3], e_Ω_YSZr, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x0], [x0], 1, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x1], [x1], 3, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x3], [x3], 4, tol=btol)
	ExtendableGrids.bfacemask!(grid, [x4], [x4], 2, tol=btol)
    grid
end