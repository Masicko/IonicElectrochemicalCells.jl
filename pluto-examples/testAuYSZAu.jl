### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0880ba1e-1f25-465e-b2e5-cc34acc7cea1
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
	using VoronoiFVM, DataFrames, GridVisualize
end

# ╔═╡ 3ad982f8-81a4-4f05-a2a2-edb374a5b347
@revise using IonicElectrochemicalCells

# ╔═╡ fd4f4ac3-f921-4ac9-a271-53a370950eec
begin
	alphas = [.30, 0.5]
	pwd()
end

# ╔═╡ 02041fc3-1078-458b-bfb2-8b5251948b0d
md"""
# Bulk vacancies, bulk electrons, no interface (MODVAL)
"""

# ╔═╡ 2c2c845a-4a96-4d50-92e8-baa0bc54cd35
Bcell = AYABoltzmann();

# ╔═╡ 69c19db2-23cd-42dc-be00-f591ae44af94
gridplot(Bcell.system.grid, Plotter=PlutoVista)

# ╔═╡ 3530836f-a613-49ac-a0c0-bc3223e39037
begin
	modval_slider = @bind modvalbiasindex html"<input type=range min=1 max=21>"
	md"""
	**Parameters**
	
	Modval bias: $(modval_slider)

	"""
end

# ╔═╡ 2fe99ebb-ef41-4d75-bde2-85867f7fccca
md"""
# Bulk vacancies, fast electrons (MODVAL with e- in equilbrium)
"""

# ╔═╡ b5a1c520-a787-444e-a087-9015e20b22c7
BNcell = AYALGBoltzmann();

# ╔═╡ cb51fb7c-ec4e-4bbb-9722-b9d24151dd44
md"""
- [ ] the weird bump around the middle (negative bias)
"""

# ╔═╡ 89f67db7-78bb-4427-893e-5e5ebe25ed23
begin
	LG_slider = @bind LGbiasindex html"<input type=range min=1 max=21 start=14>"
	md"""
	**Parameters**
	
	LG bias: $(LG_slider)

	"""
end

# ╔═╡ 94b4126d-da32-4204-b550-4ab9cbbf7b3f
md"""
## Comparison (MODVAL vs. electrons in equilibrium)
"""

# ╔═╡ 3c75a7eb-fc7d-4d3c-824c-7b8f6eb43068
begin
	comparison_slider = @bind 🐶 html"<input type=range min=1 max=21>"
	md"""
	**Parameters**
	
	Comparsion bias: $(comparison_slider)

	"""
end

# ╔═╡ 44f8cf91-ca35-4081-a28e-79ee679c3f1b
md"""
**Conclusion:** for the investigated bias range, we compared the electrostatic potential and the vacancy density of the MODVAL model and it's simplified version with equilibrated electrons. The results differ to at most 1e-14. However, this requires `tend > 1.0`.
"""

# ╔═╡ 13825791-9960-496c-ad4e-9cef09aeb7da
md"""
# Interface vacancies, fast electrons
"""

# ╔═╡ 2b75684e-3182-42a6-bb12-7b8cada8f4d9
BLG1icell = AYALG1iBoltzmann()

# ╔═╡ c58123ca-4637-4188-8509-d209f428fa54
BLG1icellpars=Dict(
		:T => 800.0,
		:alpha => 1e-1,	
		:alphas => 1e-1,
		:AYSZ => 0.38,
		:As => 0.49,
		:GA => -1.0*e0,
		:Ge => 0.0*e0,
		:SL => 1.0e0
)

# ╔═╡ 3ddeccf0-9382-442c-a987-c50447c96a9c
begin
	update_parameters!(BLG1icell, BLG1icellpars)
	BLG1icell.U = inival(BLG1icell)
	df = biasshow!(BLG1icell, tend=1.0e2)
	# cellinival = inival(BLG1icell);
end;

# ╔═╡ 7cb18a48-4012-4e52-88ec-a0ef10f1ecd2
gr = collect(1:length(BLG1icell.U[1,:]));

# ╔═╡ 58fafa95-7562-45cd-9a55-ec934615a5a9
function stateplot(cell::AYABoltzmann, U)
	p = scalarplot(gr, U[2,:], color=:orange ,Plotter=PlutoVista)
	plot!(p, gr, U[1,:], color=:blue,linewidth = 2)
	plot!(p, gr, U[3,:], linewidth = 2, color=:green,markertype=:circle, markersize=10)
	plot!(p, gr, U[4,:], color=:orange)
	return p
end

# ╔═╡ e17cc3cb-ad06-4f58-a698-8ed590ee7d6e
function stateplot(cell::AYALGBoltzmann, U)
	grid = gr #cell.system.grid
	p = scalarplot(grid, U[2,:], color=:orange ,Plotter=PlutoVista, title = "Coverage")
	plot!(p, grid ,U[1,:], color=:blue,linewidth = 2)
	#plot!(p, grid, U[3,:], linewidth = 10, color=:green,markertype=:circle, markersize=10)
	return p
end

# ╔═╡ 4942295d-7aaa-43ac-89bf-857655378809
function compareplot(U1,U2)
	p = scalarplot(gr, (U1[1,:]-U2[1,:]), color=:orange ,Plotter=PlutoVista)
	plot!(p, gr, (U1[2,:] - U2[3,:]), color=:blue,linewidth = 2)
	#plot!(p, gr, U1[2,:], linewidth = 10, color=:green,markertype=:circle, markersize=10)
	return p
end

# ╔═╡ 67b8e041-48f9-465b-a4e8-516ea4b377aa
md"""
**Observations:**

safe envelope
```
alpha, alphas
-1.0 eV < :GA < 1.0 eV
-0.05 eV < :Ge < 0.2 eV
-1.0 eV < :AYSZ < 1.0 eV
-1.0 eV < :As < 0.2 eV
1e-2 < :SL, :SR < 1e1
```
**Further investigation**

"""

# ╔═╡ 6aae8169-96cf-411c-ac28-a8e9164e415d
begin
	LG1i_slider = @bind LG1ibiasindex html"<input type=range min=1 max=21 start=11>"
	md"""
	**Parameters**
	
	LG1i bias: $(LG1i_slider)

	"""
end

# ╔═╡ 3d601f9a-e718-441f-a676-7652bdd048bd
LG1ibias = df.bias[LG1ibiasindex]

# ╔═╡ e929ca9d-b7d8-4868-b4bf-bbe44de041b4
function stateplot(cell::AYALG1iBoltzmann, U)
	p = scalarplot(gr, U[1,:], color=:blue ,Plotter=PlutoVista, title = "Coverage")
    #plot!(pp, gr, cellinival[2,:], linestyle=:dash, color=:gray)
	plot!(p, gr, U[2,:], color=:orange,linewidth = 2)
	plot!(p, gr, U[3,:], linewidth = 10, color=:green,markertype=:circle, markersize=10)
	#plot!(pp, gr, cellinival[3,:], linewidth = 10, color=:gray,markertype=:circle, markersize=10)
	return p
end

# ╔═╡ d389d5ec-923f-42c9-9903-f43d941b5a9b
stateplot(BLG1icell, df.U[LG1ibiasindex])

# ╔═╡ b5e3ac3a-627c-4ee5-84a6-6f5391044739
stateplot(BLG1icell, df.U[LG1ibiasindex])

# ╔═╡ d6d5bdfe-d730-445c-80bd-6c5a4d57cbb6
function biascapacitance!(cell :: AbstractCell; bend=1.0, bstep=0.05, tend=1e-1)
	function dcap_stepdown!(cell :: AbstractCell; dbias=1e-2, tend=tend)
    	bias = cell.system.physics.data.bias
		measurement(cell) = get_stored_charge(cell)
    	m1 = measurement(cell)
    	stationary_update!(cell, Dict(:bias => bias-dbias), tend=tend)
    	m2 = measurement(cell)
    	return (m1-m2)/dbias
		#return m1
	end #dcap_stepdown!
    cell.U = inival(cell)
    cell.Uold = copy(cell.U)
    #results = DataFrame(bias=Float64[], U=typeof(cell.U))
    results = DataFrame(bias=Float64[], U=Any[], capacitance=Float64[])
    for bias ∈ collect(0.0:bstep:bend)
        stationary_update!(cell , Dict(:bias => bias), tend=tend)
        push!(results, Dict(:bias => bias, :U => copy(cell.U), :capacitance => copy(dcap_stepdown!(cell))))
    end
    cell.U = copy(cell.Uold)
    stationary_update!(cell , Dict(:bias => 0.0), tend=tend)
    for bias ∈ collect(-bstep:-bstep:-bend)
        stationary_update!(cell , Dict(:bias => bias), tend=tend)
        push!(results, Dict(:bias => bias, :U => copy(cell.U), :capacitance => copy(dcap_stepdown!(cell))))
    end
    return sort!(results, [:bias])
end

# ╔═╡ 3f35fdfe-5ee7-4ff8-9593-34ef39f413da
begin
	update_parameters!(Bcell, Dict(:alpha => alphas[1], :alphas => alphas[2]))
	Bcell.U = inival(Bcell)
	Bdf = biascapacitance!(Bcell, tend=1e1);
end;

# ╔═╡ 061dd647-6396-4c60-abb6-39ed4579bf4a
stateplot(Bcell, Bdf.U[modvalbiasindex])

# ╔═╡ 65ba751d-b21a-46ae-bb34-a8f329f64daa
plot(Bdf.bias, Bdf.capacitance)

# ╔═╡ 053e25e7-fdbd-40a5-8867-368f7782f99c
modvalbias = Bdf.bias[modvalbiasindex]

# ╔═╡ e01b4c78-bd31-44a8-aaf4-14b20c9084e0
println(Bdf.capacitance)

# ╔═╡ e2076e94-8101-46f0-9ce4-92cc6370402b
begin
	update_parameters!(BNcell, Dict(:alpha => alphas[1], :alphas => alphas[2]))
	BNcell.U = inival(BNcell)
	BNdf = biascapacitance!(BNcell, tend=1.0e0);
end;

# ╔═╡ 16a57ed5-cdd5-4ded-b98a-f0ffcb2c4854
stateplot(BNcell, BNdf.U[LGbiasindex])

# ╔═╡ 01e85773-a85a-4ae9-bdb9-93c762359230
plot(BNdf.bias, BNdf.capacitance)
#println(BNdf.capacitance)

# ╔═╡ 2dce89ec-d39e-4a9c-8fb0-8fc7f4f722d5
compareplot(BNdf.U[🐶], Bdf.U[🐶])

# ╔═╡ 7dcbc944-c1dc-4c6b-b5ed-0accbd1ccc4a
plot(Bdf.bias, Bdf.capacitance .- BNdf.capacitance)

# ╔═╡ 0b2d3b1c-663b-4ecb-ac7e-81c26faeff80
begin
	capcell = AYALG1iBoltzmann()
	refcapdf = biascapacitance!(capcell)
	#update_parameters!(capcell, Dict(:GA => 0.0*e0, :SL => 1e0, :AYSZ => 0.0, :As => 0.0, :Ge => 0.1*e0))
	update_parameters!(capcell, BLG1icellpars)
	capdf = biascapacitance!(capcell)
end;

# ╔═╡ db3a45ba-970f-4a5b-922c-3b784aca061e
begin
	#cp = plot(capdf.bias, capdf.capacitance)
	#plot!(cp, refcapdf.bias, refcapdf.capacitance, linestyle=:dot)
	plot(refcapdf.bias, refcapdf.capacitance, linestyle=:dot)
	println(refcapdf.capacitance)
end

# ╔═╡ e1a9a96b-aada-4f64-b1b9-37c042ea8806
md"""
Observations
- `tend ~ 1e-1` in `stationary_update(tend)`
- the surface charge on the *left* ISR is 16 orders of magnitude larger than the bulk capacitance
- the surface charge on the *right* ISR is, independently of bias, 0.0
- capacitance is 10*ISR charge 
- setting `f[ipsi] = 0.0` in the `breaction()` reduces
  - ISR charge density by 2 orders of magnitude
  - ISR capacitance contribution by factor of 10
  - ...**why? It should effectively send it to 0.0**
  - unless there's something going on with the boundary integration
- alternative computation of the surface charge gives reasonable ISR capacitance, an order of magnitude larger than the bulk contribution

To-Do
- [OK] check the `bfacemasks!()` vs ISR regions (grids.jl -> AuYSZAu)
- [ ] check the dimensions of ISR lattice densities and its constants
- [fail] check the consistency of the `VoronoiFVM.integrate(boundary=true)`
- [ ] find out why is the capacitance now inverted...
"""

# ╔═╡ Cell order:
# ╠═0880ba1e-1f25-465e-b2e5-cc34acc7cea1
# ╠═3ad982f8-81a4-4f05-a2a2-edb374a5b347
# ╠═fd4f4ac3-f921-4ac9-a271-53a370950eec
# ╠═69c19db2-23cd-42dc-be00-f591ae44af94
# ╟─02041fc3-1078-458b-bfb2-8b5251948b0d
# ╠═2c2c845a-4a96-4d50-92e8-baa0bc54cd35
# ╠═3f35fdfe-5ee7-4ff8-9593-34ef39f413da
# ╠═061dd647-6396-4c60-abb6-39ed4579bf4a
# ╟─3530836f-a613-49ac-a0c0-bc3223e39037
# ╠═65ba751d-b21a-46ae-bb34-a8f329f64daa
# ╠═053e25e7-fdbd-40a5-8867-368f7782f99c
# ╟─58fafa95-7562-45cd-9a55-ec934615a5a9
# ╠═e01b4c78-bd31-44a8-aaf4-14b20c9084e0
# ╟─2fe99ebb-ef41-4d75-bde2-85867f7fccca
# ╠═b5a1c520-a787-444e-a087-9015e20b22c7
# ╠═e2076e94-8101-46f0-9ce4-92cc6370402b
# ╟─cb51fb7c-ec4e-4bbb-9722-b9d24151dd44
# ╠═16a57ed5-cdd5-4ded-b98a-f0ffcb2c4854
# ╟─89f67db7-78bb-4427-893e-5e5ebe25ed23
# ╠═01e85773-a85a-4ae9-bdb9-93c762359230
# ╟─e17cc3cb-ad06-4f58-a698-8ed590ee7d6e
# ╟─94b4126d-da32-4204-b550-4ab9cbbf7b3f
# ╠═2dce89ec-d39e-4a9c-8fb0-8fc7f4f722d5
# ╟─3c75a7eb-fc7d-4d3c-824c-7b8f6eb43068
# ╠═7dcbc944-c1dc-4c6b-b5ed-0accbd1ccc4a
# ╟─4942295d-7aaa-43ac-89bf-857655378809
# ╟─44f8cf91-ca35-4081-a28e-79ee679c3f1b
# ╟─13825791-9960-496c-ad4e-9cef09aeb7da
# ╠═2b75684e-3182-42a6-bb12-7b8cada8f4d9
# ╠═c58123ca-4637-4188-8509-d209f428fa54
# ╠═3ddeccf0-9382-442c-a987-c50447c96a9c
# ╠═7cb18a48-4012-4e52-88ec-a0ef10f1ecd2
# ╟─67b8e041-48f9-465b-a4e8-516ea4b377aa
# ╠═d389d5ec-923f-42c9-9903-f43d941b5a9b
# ╠═6aae8169-96cf-411c-ac28-a8e9164e415d
# ╟─3d601f9a-e718-441f-a676-7652bdd048bd
# ╠═db3a45ba-970f-4a5b-922c-3b784aca061e
# ╠═0b2d3b1c-663b-4ecb-ac7e-81c26faeff80
# ╠═b5e3ac3a-627c-4ee5-84a6-6f5391044739
# ╠═e929ca9d-b7d8-4868-b4bf-bbe44de041b4
# ╠═d6d5bdfe-d730-445c-80bd-6c5a4d57cbb6
# ╟─e1a9a96b-aada-4f64-b1b9-37c042ea8806
