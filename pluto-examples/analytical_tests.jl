### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 862cbd68-d0d5-11ed-2cad-89ac2bffffc9
begin
    import Pkg
    #Pkg.activate(Base.current_project())
	Pkg.activate("..")
	#
    # activate the shared project environment
    # instantiate, i.e. make sure that all packages are downloaded
    #Pkg.instantiate()
	#Pkg.resolve()
	using Markdown, Revise
	using PlutoLinks: @revise
	using PlutoVista
	using VoronoiFVM, DataFrames, GridVisualize, Parquet, ExtendableGrids
	using Plots
end

# ╔═╡ cd111a18-02a8-4e76-ab06-7d012c1692da
begin 
	@revise using IonicElectrochemicalCells
	iec = IonicElectrochemicalCells
end

# ╔═╡ 16a755fc-ec5d-4f2e-9f3a-172f279373a9
md"""
# Analytical tests and view of solution

### Parameters of the model
#### Bulk YSZ
- AYSZ => Influence of **interaction energy** between charged species in bulk YSZ
- alpha => ratio of **mobile of oxides**
  - if $0.0$ -> no mobile oxides... 
  - if $1.0$ -> all vacancies (oxides) can move

#### Surface
- AYSZs  =>  as AYSZ but for surface YSZ lattice
- alphas => as alpha but for surface YSZ lattice
- GA     => **Gibbs energy of adsorption** of vacancies from YSZ to surface
- Ge 	 => **Gibbs energy of adsorption** of electrones from Au to surface
- SL     => **area ratio** for **left** ISR between bulk structures and surface structures  -> it influences density of surface YSZ lattice sites $\underline{n^\#_\text{YSZ}} = S_\text{L}  \ (n^\#_\text{YSZ})^\frac{2}{3}$ and similarly for Au
- SR 	 => **area ratio** for **right** ISR
"""

# ╔═╡ c5ffb51f-4591-4447-9b48-518a8c87b74c
begin
	STmy_bias = 0.5
	HALF_bool = false
end

# ╔═╡ a1219ce0-ba29-4760-83db-338b850624ff
begin	
	U_storage = []; tend_storage = []; tend = 1e-2
	
	
	if HALF_bool
		STcell = iec.get_my_top_AYALG1iBoltzmann_HALF();
	else
		#STcell = iec.get_my_top_AYALG1iBoltzmann(); 
	end 

	STcell = iec.AYA_Sl()
	
	std_dict = Dict(
			:T => 800.0,

			:AYSZ => 0.0,
			:alpha => 8.0e-2,
			#:alpha => 8.0e-5,
			
			:AYSZs => 0.0,
			:alphas => 0.025,
			#:alphas => 0.0025,
			:GA => -0.1*e0,
			:Ge => -0.0*e0,
			:SL => 20.0,
			:SR => 1.0,
			:DYSZ => 1.0e-9,
			#:DYSZ => 1.0e-7,
		
			
			#:boundary_charge_fac => 0.0, 
			#:kA => 1.0e0
			)
	STcell.U = inival(STcell)
	stationary_update!(STcell,
		std_dict
		#act_prms_dict	
	)
	if STmy_bias == 0.0
		bias_list = [0.0]
	else
		bias_list = collect(0.0 : sign(STmy_bias)*0.05 : STmy_bias)
	end
	for bias ∈ bias_list
		#@show bias
		stationary_update!(STcell, Dict(:bias => bias), tend=tend)
		push!(U_storage, copy(STcell.U))
		push!(tend_storage, bias)
	end
end

# ╔═╡ 8ee98e1a-a3d7-4973-85b0-ced03996a732
begin
iiddxxx = length(tend_storage)
#@show tend_storage[iiddxxx]
l_DL = 1.0e-9
side = "L"
if side == "L"
	view_l = iec.electrode_thickness
elseif side == "R"
	view_l = iec.electrode_thickness + iec.electrolyte_thickness
end
iec.stateplot(STcell, U_storage[iiddxxx] , 
	Plotter=Plots,
	Xscale=true, #xlim=view_l .+(-l_DL, l_DL),
	#ylim = (0.99, 1.00),
	#lim = (340, 340),
	title = !HALF_bool ? "Full cell: Au | YSZ | Au" : "Half cell : Au | YSZ"
)
end

# ╔═╡ 1c71ef2b-d4d6-44c3-a850-2fbe6c76129f
begin
plot_L = iec.stateplot(STcell, U_storage[iiddxxx], 
	Plotter=Plots,
	Xscale=true, xlim=iec.electrode_thickness .+(-l_DL, l_DL),
	title = "Left DL"
)
plot_R = iec.stateplot(STcell, U_storage[iiddxxx], 
	Plotter=Plots,
	Xscale=true, xlim=iec.electrode_thickness+ iec.electrolyte_thickness .+(-l_DL, l_DL),
	title = "Right DL"
)	
Plots.plot(plot_L, plot_R)
end

# ╔═╡ 78510040-422a-4970-a512-caa583b3d90e
#iec.get_half_DL_charge(STcell, "l")

# ╔═╡ a9335d58-7109-41a7-869b-3d71fbdfec6f
#iec.get_half_DL_charge(STcell, "r")

# ╔═╡ eb13fcd5-2af0-4707-ab99-40b0dd87ac16
# if !HALF_bool
  
#   # be careful! this capacitance_measurement changes STcell, specially its bias which is important later in evaluation of ISR_chargedensity
  
# 	p=Plots.plot(xlabel="bias (V)", ylabel="C [F/m^2]")
# 	gdf = iec.capacitance_measurement!(STcell, bend=0.3, bstep=0.01, tend=1.0e-2); iec.plot_capacitance(p, gdf)
# 	p
# end

# ╔═╡ 9ddaf4e7-0082-4ef3-b24a-245d4ec1b73b
# begin 
# 	data = STcell.system.physics.data
# 	X = STcell.system.grid.components[ExtendableGrids.Coordinates]
# 	BFNodes = STcell.system.grid.components[ExtendableGrids.BFaceNodes]
# 	testing_id = b1_id = BFNodes[3]
# 	HALF_bool ? b2_id = 0 : b2_id = BFNodes[4]
	
# 	U = U_storage[iiddxxx]
# 	U_ISR = U[:, b1_id]
	
# 	bQ_YSZ_L = iec.YSZ_surface_charge(U[3, b1_id], data.alphas, data.SL)
# 	bQ_Au_L = iec.Au_surface_charge(
# 		(1/iec.nLAus(data.SL))*iec.ISR_electrondensity(
# 			U[:, b1_id],  whatever(iec.:Γ_YSZl), data 
# 		), 
# 		data.SL
# 	)

# 	get_charge_from_sys(yV) = iec.YSZ_charge_density(
# 		iec.nVmax(data.alpha) * yV
# 	)
# 	if HALF_bool
# 		b2_id = 0
# 		phi_C = U[1, length(X)-1] # the "-1" is because of weird behavior at the very end
# 		X_list_L = X[b1_id : end-1]
# 		charge_list_L = get_charge_from_sys.(U[2, b1_id : end-1])
# 		nF_YSZ_L = iec.integrate_arrays_by_trapezoid(X_list_L, charge_list_L)

# 		nF_YSZ_R = 0
# 		bQ_YSZ_R = 0
# 		bQ_Au_R = 0
# 		(Q_Au_L, nF_YSZ_C) = (-1) .* iec.bulk_charge(STcell)[1, :]
# 		Q_Au_R = 0
# 		(bQ_total_L) = (-1) .* iec.boundary_charge(STcell)[1, 3:end][1]
	
# 	else
# 		b2_id = BFNodes[4]
# 		charge_list = get_charge_from_sys.(U[2, b1_id : b2_id])
# 		X_list = X[b1_id : b2_id]
# 		center = Int(round(length(X_list)/2))
# 		phi_C = U[1, b1_id+center]
# 		nF_YSZ_L = iec.integrate_arrays_by_trapezoid(X_list[1:center], charge_list[1:center])
# 		nF_YSZ_R = iec.integrate_arrays_by_trapezoid(X_list[center:end], charge_list[center:end])

		
# 		(Q_Au_L, B, Q_Au_R) = (-1) .* iec.bulk_charge(STcell)[1, :]
# 		(bQ_total_L, bQ_total_R) = (-1) .* iec.boundary_charge(STcell)[1, 3:end][1:2]

# 		nF_YSZ_C = 0
# 		bQ_YSZ_R = iec.YSZ_surface_charge(U[3, b2_id], data.alphas, data.SR)
# 		bQ_Au_R = iec.Au_surface_charge(
# 			(1/iec.nLAus(data.SR))*iec.ISR_electrondensity(
# 				U[:, b2_id],  whatever(iec.:Γ_YSZr), data 
# 			), 
# 			data.SR
# 		)
# 	end
# 	nothing
# end

# ╔═╡ 86b548ce-6cb7-4810-9a6d-f2895208e0a5
pch = iec.get_partial_charges(STcell)

# ╔═╡ f42113d5-6c40-415b-9463-dabd3a939df3
BUCH = iec.bulk_charge(STcell)[1,:]

# ╔═╡ 0c939603-58f1-4c10-955b-44f10d5b9d77
BACH = iec.boundary_charge(STcell)[1,:]

# ╔═╡ aff3b380-0fb9-4f20-8efd-c26c9ee07198
BUCH[1] + BUCH[2] + BACH[3]

# ╔═╡ 147e6ed3-fbb0-448b-8779-9be41e85221e
BUCH[4] + BUCH[3] + BACH[4]

# ╔═╡ 59b83092-cfbc-4929-8edb-cde267a4866b
pch[:nF_Au_L] + pch[:bQ_Au_L]

# ╔═╡ 54a9bff8-2693-4b6a-bd51-e78ffd6dfec9
pch[:nF_YSZ_L] + pch[:bQ_YSZ_L]

# ╔═╡ dce1dab7-ee31-46f0-b60d-3a180e820deb
pch[:nF_YSZ_R] + pch[:bQ_YSZ_R]

# ╔═╡ ca5ecb7b-4a2e-47ff-853a-6eaf5160edf0
pch[:nF_Au_R] + pch[:bQ_Au_R]

# ╔═╡ 42071df9-7527-4347-83a3-a6504eb48824
pch[:bQ_YSZ_L] + pch[:bQ_Au_L] - BACH[3]

# ╔═╡ 3de2f510-f33e-4845-8485-7b1eae4ea8d4
pch[:bQ_YSZ_R] + pch[:bQ_Au_R] - BACH[4]

# ╔═╡ 13f6814e-392e-4f9a-b1c2-e64c95f8a569
pch[:nF_Au_L] + pch[:bQ_Au_L] + pch[:bQ_YSZ_L] + pch[:nF_YSZ_L]

# ╔═╡ 38adbd36-94fb-4357-9dfa-fe5780f1c9f1
pch[:nF_Au_R] + pch[:bQ_Au_R] + pch[:bQ_YSZ_R] + pch[:nF_YSZ_R]

# ╔═╡ d6968767-b29e-4a4c-b19c-ee4ba76b2415
pch[:bQ_YSZ_L] + pch[:nF_YSZ_L] + pch[:nF_YSZ_R] + pch[:bQ_YSZ_R]

# ╔═╡ 653a2058-2e5c-4b1c-bd1b-ec4efd2d7989
md"""
### control quantities of numerical solution on left DL
"""

# ╔═╡ 4f348a1d-52b1-42ed-b5f9-18985c685b6a
begin
	data = STcell.system.physics.data
	X = STcell.system.grid.components[ExtendableGrids.Coordinates]
	BFNodes = STcell.system.grid.components[ExtendableGrids.BFaceNodes]
	testing_id = b1_id = BFNodes[3]
	HALF_bool ? b2_id = 0 : b2_id = BFNodes[4]
	U = U_storage[iiddxxx]
	U_ISR = U[:, b1_id]

	if HALF_bool
		phi_C = U[1, length(X)-1] # the "-1" is because of weird behavior at the very
	else
		b2_id = BFNodes[4]
		X_list = X[b1_id : b2_id]
		center = Int(round(length(X_list)/2))
		phi_C = U[1, b1_id+center]
	end

	yV = U[2, testing_id]
	phi_ISR = U[1, testing_id]
	V_ISR = phi_ISR - phi_C
	dphi = U[1, testing_id+1] - U[1, testing_id]
	dx = X[testing_id+1]-X[testing_id]
	@show V_ISR
	
	println("\n <<<< Comparable quantities to analytical solution: >>>>")
	println(" --- YSZ side:")
	@show yV dphi/dx

	phi_L = U[1, 1]
	U_Au = U[1, testing_id] - phi_L
	ye_eq = iec.BoltzmannAu_ne(data, U_Au)/iec.nLAu
	dphi_Au = U[1, testing_id] - U[1, testing_id-1]
	dx_Au = X[testing_id] - X[testing_id-1]
	println(" --- Au side:")
	@show ye_eq dphi_Au/dx_Au

	
	nVs_eq = iec.nVmaxs(data.alphas, data.SL ) * U[3, testing_id]
	nes_eq = iec.ISR_electrondensity(U_ISR, iec.dummy_bnode(iec.:Γ_YSZl), data)
	yVs_eq = U[3, testing_id]
	yes_eq = (1/iec.nLAus(data.SL))*nes_eq
	println(" --- ISR:")
	@show yVs_eq yes_eq
end;

# ╔═╡ 7926148f-39c7-405b-93e1-9ab1d0921ebd
md"""
### (semi) analytical solution of left DL:
"""

# ╔═╡ 2772b5bb-83c5-4dad-b2c0-e3d3b8cab1b7
begin
	V_tot = phi_L - phi_C

	V_ISR_anal = iec.HALF_get_analytic_V_ISR(data, V_tot)
	println("\n>>> relative V_ISR error: ", (V_ISR_anal  - V_ISR)/V_ISR)
	
	iec.HALF_test_analytical(data, V_tot, V_ISR_anal, verbose=true)	
	nothing
end

# ╔═╡ 8c35921f-ac8a-49e4-b3f5-aebb37d7dc6e
md"""
# notes

- tests passed - hoorray
"""

# ╔═╡ fc88ffb4-c87f-481a-b08d-717d098409b3
md"""
## questions:

- why Dirichlet BC on HALF model cannot be set in bulk YSZ for vacancy density?
  - isnt it connected to a weird jump on the right end of electrostatic potencial  phi in YSZ domain? The phenomennon is that phi is constant on the value cca 1.0e-6 up to the very last YSZ node, where it jumps up to zero. There is physically no reason for this behavior. You can try it by changing "HALF_bool = true" and run the Pluto cell below
  - but when using the whole cell with both electrodes on sides, it seems ok
"""

# ╔═╡ 4209eb4c-104b-46eb-afc5-251f42608203
begin 
	if HALF_bool 
		println("The strange behavior of HALF cell at the right end of YSZ: \n")
		@show U[1, end-4 : end]
	end
end

# ╔═╡ Cell order:
# ╟─862cbd68-d0d5-11ed-2cad-89ac2bffffc9
# ╟─cd111a18-02a8-4e76-ab06-7d012c1692da
# ╟─16a755fc-ec5d-4f2e-9f3a-172f279373a9
# ╠═c5ffb51f-4591-4447-9b48-518a8c87b74c
# ╠═a1219ce0-ba29-4760-83db-338b850624ff
# ╟─8ee98e1a-a3d7-4973-85b0-ced03996a732
# ╟─1c71ef2b-d4d6-44c3-a850-2fbe6c76129f
# ╟─78510040-422a-4970-a512-caa583b3d90e
# ╟─a9335d58-7109-41a7-869b-3d71fbdfec6f
# ╟─eb13fcd5-2af0-4707-ab99-40b0dd87ac16
# ╠═9ddaf4e7-0082-4ef3-b24a-245d4ec1b73b
# ╠═86b548ce-6cb7-4810-9a6d-f2895208e0a5
# ╠═f42113d5-6c40-415b-9463-dabd3a939df3
# ╠═0c939603-58f1-4c10-955b-44f10d5b9d77
# ╠═aff3b380-0fb9-4f20-8efd-c26c9ee07198
# ╠═147e6ed3-fbb0-448b-8779-9be41e85221e
# ╠═59b83092-cfbc-4929-8edb-cde267a4866b
# ╠═54a9bff8-2693-4b6a-bd51-e78ffd6dfec9
# ╠═dce1dab7-ee31-46f0-b60d-3a180e820deb
# ╠═ca5ecb7b-4a2e-47ff-853a-6eaf5160edf0
# ╠═42071df9-7527-4347-83a3-a6504eb48824
# ╠═3de2f510-f33e-4845-8485-7b1eae4ea8d4
# ╠═13f6814e-392e-4f9a-b1c2-e64c95f8a569
# ╠═38adbd36-94fb-4357-9dfa-fe5780f1c9f1
# ╠═d6968767-b29e-4a4c-b19c-ee4ba76b2415
# ╟─653a2058-2e5c-4b1c-bd1b-ec4efd2d7989
# ╟─4f348a1d-52b1-42ed-b5f9-18985c685b6a
# ╟─7926148f-39c7-405b-93e1-9ab1d0921ebd
# ╟─2772b5bb-83c5-4dad-b2c0-e3d3b8cab1b7
# ╟─8c35921f-ac8a-49e4-b3f5-aebb37d7dc6e
# ╟─fc88ffb4-c87f-481a-b08d-717d098409b3
# ╟─4209eb4c-104b-46eb-afc5-251f42608203
