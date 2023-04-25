function compose_dict(this)
    d = Dict()
    for name in fieldnames(typeof(this))
        d[name] = getfield(this,name)
    end
    return d
end

function convert_to_strings(old_syms)
    new_strings=Dict()
    for (key,value) in old_syms
        new_strings[string(key)] = value
    end
    return new_strings
end

function set_parameters!(parameters, d::Dict)
    for key in keys(d)
        key in fieldnames(typeof(parameters)) ? setfield!(parameters, Symbol(key), d[key]) : @warn key "not found in parameters"
    end
end

function ishless(a::Real,b::Real)
    c = a*b
    if c < 0.0
        if a < 0.0
            return false
        else
            return true
        end
    else
        return isless(abs(a), abs(b))
    end
end

function row2dict(row)
    return Dict(Symbol.(names(row)) .=> values(row))
end

function parameters2VoronoiData(model, name :: Symbol)
    # https://discourse.julialang.org/t/create-parametric-structure-with-metaprogramming/77184
    fields = [:( $(nm)::T ) for nm in model[:fieldname]]
    eval(quote
        mutable struct $name{T}
            $(fields...)
        end
        $name() = $name($model[:val]...)
    end
    )
end

function integrate_arrays_by_trapezoid(x, y)
	sum = 0
	for i in 1:length(x)-1
		sum += (x[i+1]-x[i])*(y[i]+y[i+1])/2
	end
	return sum
end