module Data

using StaticArrays


function inenout(data)
    
    part_capt = 0
    part_nocapt = 0

    for part in data
        final_pos = part.pos[end]

        if final_pos[2] < 0
            part_capt += 1
        else
            part_nocapt += 1

        end
    end
    return part_capt, part_nocapt
end

    
end