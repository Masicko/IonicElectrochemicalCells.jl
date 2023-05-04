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
	STmy_bias = 1.0
	HALF_bool = false
end

# ╔═╡ a1219ce0-ba29-4760-83db-338b850624ff
begin	
	U_storage = []; tend_storage = []; tend = 1e-2
	
	
	if HALF_bool
		STcell = iec.get_my_top_AYALG1iBoltzmann_HALF();
	else
		STcell = iec.get_my_top_AYALG1iBoltzmann(); 
	end 

	STcell = iec.AYA_Sl()
	
	std_dict = Dict(
			:T => 800.0,

			:AYSZ => 0.0,
			#:alpha => 8.0e-3,
			:alpha => 8.0e-5,
			
			:AYSZs => 0.0,
			#:alphas => 0.025,
			:alphas => 0.0025,
			:GA => 0.1*e0,
			:Ge => -0.0*e0,
			:SL => 0.1,
			:SR => 0.5,
			:DYSZ => 1.0e-7
		
			
			#:boundary_charge_fac => 0.0,
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

# ╔═╡ 7131d271-145d-4949-986e-0ad1583f2060
STcell.system.physics.data.DYSZ

# ╔═╡ 8ee98e1a-a3d7-4973-85b0-ced03996a732
begin
iiddxxx = length(tend_storage)
#@show tend_storage[iiddxxx]
l_DL = 30.0e-9
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
iec.get_half_DL_charge(STcell, "l")

# ╔═╡ eb13fcd5-2af0-4707-ab99-40b0dd87ac16
if !HALF_bool
  
  # be careful! this capacitance_measurement changes STcell, specially its bias which is important later in evaluation of ISR_chargedensity
  
	p=Plots.plot(xlabel="bias (V)", ylabel="C [F/m^2]")
	gdf = iec.capacitance_measurement!(STcell, bend=0.3, bstep=0.01, tend=1.0e-2); iec.plot_capacitance(p, gdf)
	p
end

# ╔═╡ 751270c1-3517-49b0-ad3f-e8cb86b9c6e9
HALF_bool ? nothing : iec.get_stored_charge(STcell)

# ╔═╡ 32f7d00a-fd3b-4de8-bad0-4e0dfbc1c5ac
buc = iec.bulk_charge(STcell)

# ╔═╡ c7d5000c-24ab-4c10-b9ab-9bd800a3cb94
bac = iec.boundary_charge(STcell)

# ╔═╡ 88a19671-faf7-4d8a-b579-6e9a5bf8cee2
buc[1, 2] + buc[1, 4] + bac[1, 3] + bac[1, 4]

# ╔═╡ eb734e70-844c-4f41-aed0-e4c2c0548a93
begin 
mutable struct whatever
	region
end
bnode = whatever(iec.:Γ_YSZl)
end

# ╔═╡ 9ddaf4e7-0082-4ef3-b24a-245d4ec1b73b
begin 
	data = STcell.system.physics.data
	X = STcell.system.grid.components[ExtendableGrids.Coordinates]
	BFNodes = STcell.system.grid.components[ExtendableGrids.BFaceNodes]
	testing_id = b1_id = BFNodes[3]
	U = U_storage[iiddxxx]
	U_ISR = U[:, testing_id]
	

	bQ_YSZ_L = iec.YSZ_surface_charge(U[3, b1_id], data.alphas, data.SL)
	bQ_Au_L = iec.Au_surface_charge(
		(1/iec.nLAus(data.SL))*iec.ISR_electrondensity(
			U[:, b1_id],  whatever(iec.:Γ_YSZl), data 
		), 
		data.SL
	)

	get_charge_from_sys(yV) = iec.YSZ_charge_density(
		iec.nVmax(data.alpha) * yV
	)
	if HALF_bool
		b2_id = 0
		phi_C = U[1, length(X)-1] # the "-1" is because of weird behavior at the very end
		X_list_L = X[b1_id : end-1]
		charge_list_L = get_charge_from_sys.(U[2, b1_id : end-1])
		nF_YSZ_L = iec.integrate_arrays_by_trapezoid(X_list_L, charge_list_L)

		nF_YSZ_R = 0
		bQ_YSZ_R = 0
		bQ_Au_R = 0
		(Q_Au_L, nF_YSZ_C) = (-1) .* iec.bulk_charge(STcell)[1, :]
		Q_Au_R = 0
		(bQ_total_L) = (-1) .* iec.boundary_charge(STcell)[1, 3:end][1]
	
	else
		b2_id = BFNodes[4]
		charge_list = get_charge_from_sys.(U[2, b1_id : b2_id])
		X_list = X[b1_id : b2_id]
		center = Int(round(length(X_list)/2))
		phi_C = U[1, b1_id+center]
		nF_YSZ_L = iec.integrate_arrays_by_trapezoid(X_list[1:center], charge_list[1:center])
		nF_YSZ_R = iec.integrate_arrays_by_trapezoid(X_list[center:end], charge_list[center:end])

		
		(Q_Au_L, B, Q_Au_R) = (-1) .* iec.bulk_charge(STcell)[1, :]
		(bQ_total_L, bQ_total_R) = (-1) .* iec.boundary_charge(STcell)[1, 3:end][1:2]

		nF_YSZ_C = 0
		bQ_YSZ_R = iec.YSZ_surface_charge(U[3, b2_id], data.alphas, data.SR)
		bQ_Au_R = iec.Au_surface_charge(
			(1/iec.nLAus(data.SR))*iec.ISR_electrondensity(
				U[:, b2_id],  whatever(iec.:Γ_YSZr), data 
			), 
			data.SR
		)
	end
	nothing
end

# ╔═╡ a0fefc76-fa61-4963-8ce3-f33122967e8a
"charge on YSZ", sum([bQ_YSZ_L, nF_YSZ_L, nF_YSZ_R, bQ_YSZ_R])

# ╔═╡ 5ea88725-75be-4bd4-a4e8-d6f103fe3e6a
[bQ_YSZ_L, nF_YSZ_L, nF_YSZ_R, bQ_YSZ_R]

# ╔═╡ 8d2ce173-2a5c-4dde-9fe4-82d518eab8b6
bQ_YSZ_L + nF_YSZ_C

# ╔═╡ 3896ca24-fc72-4f0f-bc13-6d98280c8b8c
Q_Au_L + bQ_Au_L

# ╔═╡ 8302770a-4592-4a6d-9dd4-b07956b6a621
bQ_YSZ_L + nF_YSZ_L

# ╔═╡ 4cb142f6-ef3e-4a37-9a16-3e89aed5eb79
bQ_YSZ_R + nF_YSZ_R

# ╔═╡ 388b8711-d39b-4c69-8e79-fe68a5f0be4e
Q_Au_R + bQ_Au_R

# ╔═╡ 1809400f-7a9e-44ff-80e7-691eccc0bc98
"DL L all", [Q_Au_L, bQ_Au_L, bQ_YSZ_L, nF_YSZ_L]

# ╔═╡ 7c8f359c-e98c-4b44-958d-efca9cc63da0
"DL L l,r", sum([Q_Au_L, bQ_Au_L]), sum([bQ_YSZ_L, nF_YSZ_L])

# ╔═╡ 4798fa96-2887-4765-bff5-f1b1c37fb259
"DL L",sum([Q_Au_L, bQ_Au_L, bQ_YSZ_L, nF_YSZ_L])

# ╔═╡ b87f717e-1be7-4911-8964-de08e33af921
"DL R", sum([nF_YSZ_R, bQ_YSZ_R, bQ_Au_R, Q_Au_R])

# ╔═╡ be613a89-10c2-4ddd-9844-bf1b472585d4
"naboj na YSZ", sum([bQ_YSZ_L, nF_YSZ_L, nF_YSZ_R, bQ_YSZ_R])

# ╔═╡ 653a2058-2e5c-4b1c-bd1b-ec4efd2d7989
md"""
### control quantities of numerical solution on left DL
"""

# ╔═╡ 4f348a1d-52b1-42ed-b5f9-18985c685b6a
begin
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
	nes_eq = iec.ISR_electrondensity(U_ISR, bnode, data)
	yVs_eq = U[3, testing_id]
	yes_eq = (1/iec.nLAus(data.SR))*nes_eq
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
# ╠═7131d271-145d-4949-986e-0ad1583f2060
# ╟─8ee98e1a-a3d7-4973-85b0-ced03996a732
# ╟─1c71ef2b-d4d6-44c3-a850-2fbe6c76129f
# ╠═78510040-422a-4970-a512-caa583b3d90e
# ╠═eb13fcd5-2af0-4707-ab99-40b0dd87ac16
# ╠═9ddaf4e7-0082-4ef3-b24a-245d4ec1b73b
# ╠═a0fefc76-fa61-4963-8ce3-f33122967e8a
# ╠═5ea88725-75be-4bd4-a4e8-d6f103fe3e6a
# ╠═8d2ce173-2a5c-4dde-9fe4-82d518eab8b6
# ╠═3896ca24-fc72-4f0f-bc13-6d98280c8b8c
# ╠═8302770a-4592-4a6d-9dd4-b07956b6a621
# ╠═4cb142f6-ef3e-4a37-9a16-3e89aed5eb79
# ╠═388b8711-d39b-4c69-8e79-fe68a5f0be4e
# ╠═751270c1-3517-49b0-ad3f-e8cb86b9c6e9
# ╠═1809400f-7a9e-44ff-80e7-691eccc0bc98
# ╠═7c8f359c-e98c-4b44-958d-efca9cc63da0
# ╠═4798fa96-2887-4765-bff5-f1b1c37fb259
# ╠═b87f717e-1be7-4911-8964-de08e33af921
# ╠═be613a89-10c2-4ddd-9844-bf1b472585d4
# ╠═32f7d00a-fd3b-4de8-bad0-4e0dfbc1c5ac
# ╠═c7d5000c-24ab-4c10-b9ab-9bd800a3cb94
# ╠═88a19671-faf7-4d8a-b579-6e9a5bf8cee2
# ╟─eb734e70-844c-4f41-aed0-e4c2c0548a93
# ╟─653a2058-2e5c-4b1c-bd1b-ec4efd2d7989
# ╟─4f348a1d-52b1-42ed-b5f9-18985c685b6a
# ╟─7926148f-39c7-405b-93e1-9ab1d0921ebd
# ╟─2772b5bb-83c5-4dad-b2c0-e3d3b8cab1b7
# ╟─8c35921f-ac8a-49e4-b3f5-aebb37d7dc6e
# ╟─fc88ffb4-c87f-481a-b08d-717d098409b3
# ╟─4209eb4c-104b-46eb-afc5-251f42608203
