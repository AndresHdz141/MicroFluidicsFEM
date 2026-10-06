module ask_data
using DataFrames

export ask_array_type, ask_simulate_again,select_particle,ask_simulate_new_field

function ask_array_type()
    # A tuple where index 1 is 0, index 2 is 1, etc. (Julia is 1-indexed)
    names = ("Unipolar.jld2", "Alterno.jld2", "Halbach.jld2")
    if bytesavailable(stdin) > 0
        read(stdin, Char) 
    end
    while true
        println("Selecciona el tipo de arreglo: ")
        println("0 - Unipolar | 1 - Alterno | 2 - Halbach")
        input = strip(readline())
        
        if isempty(input)
            println("Input Invalido")
            continue
        end
        parsed = tryparse(Int, input)        
        # Since Julia arrays start at 1, option 0 corresponds to index 1, option 1 to index 2, etc.
        if parsed !== nothing && parsed in (0, 1, 2)
            return names[parsed + 1] 
        end
        
        println("Input Invalido")
    end
end
function ask_simulate_new_field()
    if bytesavailable(stdin) > 0
        read(stdin, Char) 
    end
    while true
            println("¿Simular nuevos campos fem?")
            println("(y)/(n)")
            input = strip(readline())
            if isempty(input)
                println("Input Invalido")
                continue
            end
            if input !==nothing
                if input == "y"
                    return true
                elseif input == "n"
                    return false
                end
            end
        println("Input Invalido")
    end

end

function ask_simulate_again()
    if bytesavailable(stdin) > 0
        read(stdin, Char) 
    end
    while true
            println("¿Simular otra particula en el mismo canal?")
            println("(y)/(n)")
            input = strip(readline())
            if isempty(input)
                println("Input Invalido")
                continue
            end
            if input !==nothing
                if input == "y"
                    return true
                elseif input == "n"
                    return false
                end
            end
        println("Input Invalido")
    end

end


function select_particle(data)
    while true
    println("=== PARTICULAS DISPONIBLES POR SIMULAR ===")
    for r in eachrow(data)
        if r.tag!= 0 
            println("Tag [$(r.tag)] -> Nombre: $(r.Nombre)") 
        end
    end
    println("-----------------------------------------")
    println("Seleccione el 'Tag' de la particula para simular:")

    input = strip(readline())

    try
        target_tag = parse(Int, input)
        if target_tag == 0 
            println("Error:Tag $target_tag no existe en el dataset \n")
            continue
        end
        matching_rows = filter(row -> row.tag == target_tag, data)
        
        if !isempty(matching_rows)
            println("Exito: Se simulara la particula con el tag $target_tag:")
            extracted_vector = Vector(matching_rows[1, :])
            return extracted_vector
        else
        println("\n Error: Tag $target_tag no existe en el dataset \n")
        end
    catch e
        println("\n Error: Seleccione un valor valido \n")
    end
end


end


end