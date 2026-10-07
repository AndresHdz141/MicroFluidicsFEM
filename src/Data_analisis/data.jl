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
        
        # 1. LLEGARON AL FINAL (Recuperamos tu tolerancia de - 1e-4)
        if final_pos[1] >= (x_exit - 1e-4)
            
            if final_pos[2] < mitad
                # CAPTURA PERFECTA
                total_score += 1.0
                part_capt += 1
            else
                # SALIÓ POR ARRIBA (Protegemos que la penalización no pase de 1.0)
                dist_a_mitad = final_pos[2] - mitad
                penalizacion = clamp(max(0.0, dist_a_mitad) / mitad, 0.0, 1.0)
                total_score += 1.0 - (0.5 * penalizacion) 
            end
            
        # 2. CHOCÓ O SE ATORÓ EN EL CAMINO
        else
            # Previene puntajes negativos si la partícula fue empujada en reversa
            distancia_recorrida = max(0.0, final_pos[1] - x_start)
            
            # Otorga puntos según qué tanto avanzó (máximo 0.5)
            avance_porcentual = clamp(distancia_recorrida / total_length, 0.0, 1.0)
            total_score += 0.5 * avance_porcentual
        end
    end
    
    pseudo_efficiency = total_score / length(data)
    
    return pseudo_efficiency, part_capt
end

end