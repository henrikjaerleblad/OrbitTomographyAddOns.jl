"""
    add_noise(S,0.05)
    add_noise(-||-, k=0.15)

A function that adds noise to a signal in a consistent way. The function adds both background noise as well as
noise to the signal itself. The level of the background noise and the signal noise is determined by the b
and k variables respectively. The levels are input as a fraction of the mean of the original signal S strength.
"""
function add_noise(s, b; k=0.1)
    sn = max.(s,0.0) .+ k.*(mean(sqrt.(abs.(s)))).*rand.(Normal.(0.0, max.(sqrt.(max.(s,0)), sqrt.(b))))
    err = k.*mean(sqrt.(abs.(s))).*max.(sqrt.(max.(s,0)), sqrt.(b))
    return sn, err
end

"""
    get_orbel_volume(myOrbitGrid, os_equidistant)

Get the orbit element volume of this specific orbit-grid.
If os_equidistant=false, then return a 3D array with all the orbit volumes.
Otherwise, just a scalar.
"""
function get_orbel_volume(og::OrbitGrid, os_equidistant::Bool)

    if os_equidistant
        dRm = abs(og.r[2]-og.r[1]) # The difference between the first and second element should be representative for the whole grid, due to equidistancy
        dE = abs(og.energy[2]-og.energy[1])
        dpm = abs(og.pitch[2]-og.pitch[1])
        dO = dE*dpm*dRm
        return dO
    else
        dO = zeros(length(og.energy), length(og.pitch), length(og.r)) # If not equidistant, create a 3D array with zeros
        for Ei=1:length(og.energy) # For all energies...
            if Ei==length(og.energy) # If at the edge of the orbit grid...
                # Assume edge orbit-element volume to be same as next-to-edge
                dEi = abs(og.energy[end]-og.energy[end-1])
            else
                dEi = abs(og.energy[Ei+1]-og.energy[Ei])
            end
            for pmi=1:length(og.pitch)
                if pmi==length(og.pitch)
                    dpmi = abs(og.pitch[end]-og.pitch[end-1])
                else
                    dpmi = abs(og.pitch[pmi+1]-og.pitch[pmi])
                end
                for Rmi=1:length(og.r)
                    if Rmi==length(og.r)
                        dRmi = abs(og.r[end]-og.r[end-1])
                    else
                        dRmi = abs(og.r[Rmi+1]-og.r[Rmi])
                    end

                    dO[Ei, pmi, Rmi] = dEi*dpmi*dRmi # Replace the zero-element in the 3D array with the resulting orbit-element volume
                end
            end
        end

        return dO
    end
end