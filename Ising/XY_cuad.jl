__precompile__()

module XY_cuad
export simulacion

function configuration_run(L::Int64, S::Array{Float64,2}, T::Float64, n::Array{Int64,1},
    E_sys::Float64, Mag::Array{Float64,1}, doc_path::AbstractString)

    N = L*L
    Steps = n[1]*N
    Steps_prom = n[2]*N

    Corr_array = zeros(Float64, div(L,2)+1)

    # Informacion sobre las energias y magnetizaciones
    datos = zeros(Float64, 3)
    mag_prom = [0.0, 0.0]
    folder_path = doc_path*@sprintf("/T_%.4f",T)
    for step in 1:Steps
        # Se calcula el cambio de un spin aleatorio
        i = rand(1:N)
        sp = 2.0*pi*rand()

        # Se calcula el cambio de energia debido al cambio del espin
        DeltaE = cos(S[i]-S[div(i-1,L)*L+mod(i,L)+1]) + cos(S[i]-S[mod(i-1+L,N)+1])
        DeltaE += cos(S[i]-S[div(i-1,L)*L+mod(i-2,L)+1]) + cos(S[i]-S[mod(i-1-L,N)+1])

        DeltaE -= cos(sp-S[div(i-1,L)*L+mod(i,L)+1]) + cos(sp-S[mod(i-1+L,N)+1])
        DeltaE -= cos(sp-S[div(i-1,L)*L+mod(i-2,L)+1]) + cos(sp-S[mod(i-1-L,N)+1])

        # Se determina la probabilidad de que se acepte el cambio
        if rand() < min(1.0, exp(-DeltaE/T))
            # Se actualiza el estado del sistema
            E_sys += DeltaE  # Se cambia la energia del sistema
            Mag += [cos(sp)-cos(S[i]), sin(sp)-sin(S[i])]   # Se cambia la magnetizacion del sistema
            S[i] = sp        # Se actualiza el espin del sitio i,j
        end

        if step > Steps-Steps_prom
            # Se calculan los promedios de la energia y la magnetizacion del sistema
            # Se calculan los promedios de los cuadrados la energia y la magnetizacion del sistema
            datos += [E_sys/Steps_prom, E_sys*E_sys/Steps_prom, norm(Mag)*norm(Mag)/Steps_prom]
            mag_prom += Mag/Steps_prom

            if mod(step,N) == 0
                for r in 0:div(L,2)
                    Corr = cos(circshift(S,(r,0))-S)+cos(circshift(S,(-r,0))-S)
                    Corr += cos(circshift(S,(0,r))-S)+cos(circshift(S,(0,-r))-S)
                    Corr_array[r+1] += 0.25*mean(Corr)
                end
            end
        end
    end

    # Se calculan las susceptibilidades
    Cv = (datos[2] - datos[1]*datos[1])/(T*T)
    Chi_m = (datos[3] - norm(mag_prom)*norm(mag_prom))/T

    Corr_array = Corr_array/n[2]

    return datos[1], mag_prom/N, Cv, Chi_m, Corr_array
end

function simulacion(L::Int64,Ti::Float64,Tf::Float64,Tstep::Float64,n::Array{Int64,1}, doc_path::AbstractString)
    @printf "Temperaturas -> %.4f - %.4f en pasos de %.4f\n" Ti Tf Tstep
    @printf "Red de %dx%d\n" L L
    @printf "Pasos Montecarlo -> %d\n" n[1]

    # Se define el tamano de los arreglos y los arreglos
    Temperatures = collect(Ti:Tstep:Tf)   # La temperatura

    data_folder = doc_path*@sprintf("/XY_DATA_cuadrado_L%d_MC%d_T%.4f_%.4f", L, n[1], Ti, Tf)
    try
        mkdir(data_folder)
    end
    println("Se creo la carpeta "*data_folder)

    # Se defina el path donde se van a guardar los datos
    doc_data = data_folder*"/data.csv"
    doc_logs = data_folder*"/logs.csv"
    doc_corr = data_folder*"/correlacion.csv"

    out_data = open(doc_data, "w")
    out_logs = open(doc_logs, "w")
    out_corr = open(doc_corr, "w")

    write(out_logs, @sprintf("Tiempo inicial %s \n", string(now())))
    write(out_logs, @sprintf("L=%d - Mc=%d\n", L, n[1]))
    write(out_data, @sprintf("T,E,M_x,M_y,Cv,Chim\n"))

    r_s = ""
    for r in 0:div(L,2)
        r_s *= @sprintf(",%d", r)
    end
    write(out_corr, "T"*r_s*"\n")

    flush(out_logs); flush(out_data); flush(out_corr)
    close(out_logs); close(out_data); close(out_corr)

    # Se crea una matriz de spines comun a todas las temperaturas
    S_ini = 2.0*pi*rand(L,L)
    A = cos(S_ini-circshift(S_ini,(1,0))) + cos(S_ini-circshift(S_ini,(-1,0)))
    A += cos(S_ini-circshift(S_ini,(0,1)))+ cos(S_ini-circshift(S_ini,(0,-1)))

    # Se divide el producto de S con A entre 2 ya que se cuenta la
    # interaccion 2 veces cuando se realiza el producto
    E_sys = sum(A)/2.0
    Mag = [sum(cos(S_ini)), sum(sin(S_ini))]

    @sync @parallel for t in 1:size(Temperatures)[1]
        tic()
        S = copy(S_ini)
        E, M, Cv, Chi_m, Correl_array = configuration_run(L, S, Temperatures[t], n, E_sys, Mag, data_folder)
        tot_time = toq()
        out_data = open(doc_data, "a")
        write(out_data, @sprintf("%.5f,%.10f,%.10f,%.10f,%.10f,%.10f\n", Temperatures[t], E, M[1], M[2], Cv, Chi_m))
        flush(out_data); close(out_data)

        corr_str = @sprintf("%.5f", Temperatures[t])
        for r in 1:div(L,2)+1
            corr_str *= @sprintf(",%.5f", Correl_array[r])
        end
        corr_str *= "\n"
        out_corr = open(doc_corr, "a")
        write(out_corr, corr_str)
        flush(out_corr); close(out_corr)

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
