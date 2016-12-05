__precompile__()

module Ising_cuad_images
export simulacion_images

using PyPlot

function configuration_run_images(Spins::Array{Float64,1}, L::Int64, S::Array{Float64,2}, T::Float64, n::Array{Int64,1},
    doc_path::AbstractString)

    N = L*L
    Steps = n[1]*N

    figure()
    axis("off")

    # Informacion sobre las energias y magnetizaciones
    folder_path = doc_path*@sprintf("/T_%.4f",T)

    for step in 0:Steps-1
        # Se calcula el cambio de un spin aleatorio
        i = rand(1:N)
        sp = rand(Spins)
        # Se calcula el cambio de energia debido al cambio del espin
        DeltaE = S[div(i-1,L)*L+mod(i,L)+1] + S[mod(i-1+L,N)+1]
        DeltaE += S[div(i-1,L)*L+mod(i-2,L)+1] + S[mod(i-1-L,N)+1]
        DeltaE *= S[i]-sp

        # Se determina la probabilidad de que se acepte el cambio
        if rand() < min(1.0, exp(-DeltaE/T))
            # Se actualiza el estado del sistema
            S[i] = sp        # Se actualiza el espin del sitio i,j
        end


        if mod(step, N*n[2]) == 0
            imshow(S)
            title(@sprintf("Paso MC = %d", step/N))
            if step==0
                colorbar()
                try
                    mkdir(folder_path)
                    println("Se creo la carpeta "*folder_path)
                end
            end
            savefig(folder_path*@sprintf("/plot_%.4f_%d.jpg", T, div(step,N)))
        end

    end
end

function simulacion_images(spin::Float64,L::Int64,Ti::Float64,Tf::Float64,Tstep::Float64,n::Array{Int64,1}, doc_path::AbstractString)
    @printf "Temperaturas -> %.4f - %.4f en pasos de %.4f\n" Ti Tf Tstep
    @printf "Red de %dx%d\n" L L
    @printf "Pasos Montecarlo -> %d\n" n[1]

    # Se definen los espines posibles del Sistema y la lista de posibles combinaciones
    Spins = collect(-2*spin:2:2*spin)
    @show Spins

    # Se define el tamano de los arreglos y los arreglos
    Temperatures = collect(Ti:Tstep:Tf)   # La temperatura

    data_folder = doc_path*@sprintf("/Ising_IMAGES_cuadrado_L%d_MC%d_T%.4f_%.4f", L, n[1], Ti, Tf)
    try
        mkdir(data_folder)
    end
    println("Se creo la carpeta "*data_folder)

    # Se defina el path donde se van a guardar los datos
    doc_logs = data_folder*"/logs.csv"
    out_logs = open(doc_logs, "w")

    write(out_logs, @sprintf("Tiempo inicial %s \n", string(now())))
    write(out_logs, @sprintf("L=%d - Mc=%d\n", L, n[1]))

    flush(out_logs);
    close(out_logs);

    # Se crea una matriz de spines comun a todas las temperaturas
    S_ini = rand(Spins,L,L)
    A = circshift(S_ini,(1,0))+circshift(S_ini,(-1,0))+circshift(S_ini,(0,1))+circshift(S_ini,(0,-1))


    @sync @parallel for t in 1:size(Temperatures)[1]
        tic()
        S = copy(S_ini)
        configuration_run_images(Spins, L, S, Temperatures[t], n, data_folder)
        tot_time = toq()

        out_logs = open(doc_logs, "a")
        write(out_logs, @sprintf("%.4f\n", tot_time))
        flush(out_logs); close(out_logs)

        # Se imprime el progreso del programa
        @printf "%f\n" Temperatures[t]
    end
    out_logs = open(doc_logs, "a")
    write(out_logs, @sprintf("Tiempo final %s \n", string(now())))
    flush(out_logs)
    close(out_logs)
end


end
