using ComponentArrays

"""
Equivalent to 
vortex_biot_savart!( D, u, γ, t, u; ...)

"""
function vortex_biot_savart!( D, u, γ, t; kwargs... )
    vortex_biot_savart!( D, u, γ, t, u; kwargs...)
end

"""
function vortex_biot_savart!( D, u, γ, t, w; core = 0.0)

Use Biot--Savart to compute the velocity field induced by vortices (with locations in w.x, w.y and circulations γ) onto locations u.

Shape of w.x, w.y, and γ is enforced to be the same.

Time t is the parameter required by the format of ODEProblem. It is ignored (for now).

core is the regularization parameter (added to the distance to avoid the velocity blowup).

If u ≡ w, e.g., when evolving the locations of vortices under their field, omit the w parameter.
"""
function vortex_biot_savart!( D, u, γ, t, w; core = 0.0)

    @assert size(γ) == size(w.x) == size(w.y)
    # perp of vectors connecting pairs of points in u, w
    P =  ComponentArray( 
        x =   -( u.y .- transpose(w.y) ),
        y =   u.x .- transpose(w.x),
    )


    # (! perp operation doesn't change the length)=
    dist2 = @. P.x^2 + P.y^2

    # components of the velocity field induced by each vortex
    for i = [:x,:y]
        local Pc = @view(P[i])
        Pc .= Pc ./ (dist2 .+ core^2)
        Pc .= Pc .* transpose(γ)

        Pc .*= (1/2/π)
        Pc[dist2 .== 0] .= 0 # ignore self-contribution
        @view(D[i]) .= sum(Pc,dims=2)
    end

    return nothing
end

"""
Allocating version.
"""
function vortex_biot_savart( u, args...; kwargs... )

    D = similar(u)
    vortex_biot_savart!( D, u, args...; kwargs...)

    return D
end


"""
function vorticity_model( w, γ, t; model=:burgers, ν = 1e-3 )

Returns a function p -> ω(p) that evaluates the chosen vorticity model at a specific gridpoint.

"""
function vorticity_burgers( γ, t, w; ν = 1e-3 )

    return p -> [  (γi) /(2 * π * ν) * exp( -(1/2/ν) * hypot(x-p.x, y-p.y).^2 ) 
                for (x, y, γi) in zip(w.x, w.y, γ) 
                ] |> sum

end

function velocity( γ, t, w; kwargs... )

    return p -> vortex_biot_savart( p, γ, t, w; kwargs...)

end