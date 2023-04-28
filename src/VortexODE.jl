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

        Pc .*= (-1/2/π)
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

function vortex_biot_savart( u::Tuple{LinRange,LinRange}, args...; kwargs... )
    z0 = zeros( length(u[1]) * length(u[2]))
    D = ComponentArray(x=similar(z0), y=similar(z0))
    vortex_biot_savart!( D, u, args...; kwargs...)

    return D
end


function vortex_biot_savart!( D, u::Tuple{LinRange,LinRange}, args...; kwargs...)


p = [ ComponentVector(x=x, y=y) for x in u[1], y in u[2]]
p = ComponentArray( 
    x= vec( map( pp -> pp.x, p) ),  
    y= vec( map( pp -> pp.y, p) ) )


vortex_biot_savart!(D, p, args...; kwargs...)

D.x = reshape( D.x, length.(u))
D.y = reshape( D.y, length.(u))

end