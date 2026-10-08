module Data

using StaticArrays
export inenout
function inenout(data, height, env)
    total_score = 0.0
    part_capt = 0   

    x_start = env.x_min_c
    x_exit = env.x_max_c
    mitad = height / 2.0
    total_length = x_exit - x_start

    for part in data
        final_pos = part.pos[end]
        
        # 1. CAPTURA PERFECTA (¡Premio mayor!)
        if final_pos[1] >= (x_exit - 1e-4) && final_pos[2] < mitad
            total_score += 1.0
            part_capt += 1
            
        # 2. SALIÓ POR ARRIBA (IA perezosa / Campo débil)
        elseif final_pos[1] >= (x_exit - 1e-4) && final_pos[2] >= mitad
            # 'd' va de 0.0 (rozando la meta) a 1.0 (pegado al techo)
            d = clamp((final_pos[2] - mitad) / mitad, 0.0, 1.0)
            
            # Caída exponencial severa: d=1 (techo) da ~0.005 pts. d=0 da 0.8 pts.
            score = 0.4 * exp(-5.0 * d)
            total_score += score
            
        # 3. CHOCÓ / SE ATORÓ (Campo fuerte)
        else
            avance = clamp(max(0.0, final_pos[1] - x_start) / total_length, 0.0, 1.0)
            
            # Crecimiento exponencial cúbico: avance=10% da ~0.001. avance=90% da ~0.58
            score = 0.4 * (avance ^ 3)
            total_score += score
        end
    end
    
    pseudo_efficiency = total_score / length(data)
    
    return pseudo_efficiency, part_capt
end
end