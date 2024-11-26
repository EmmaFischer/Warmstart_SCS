using JuMP
import SCS
using LinearAlgebra
import Random
using CSV, DataFrames

function ws1(delta_A, delta_b, delta_C, A_bar, x, y, s, n, m, modus = 0) #Eingabe: Optimallösung von model_1, modus; Ausagbe: Warmstartpunkt für model_2 nach Methode 1 und je nach Modus
    #Umformen der Eingaben
    delta_A = tensor_to_matrix(delta_A, n, m)
    A_bar = tensor_to_matrix(A_bar, n, m)
    delta_c = matrix_to_svec(delta_C, n) 
    x = jumpvec_to_svec(x, n)
    s = jumpvec_to_svec(s, n)

    #Definieren der Dimension von x
    d = Int(n*(n+1)/2)

    #Least-SquareAdjustment
    #Berchnen von delta_x, delta_y, delta_s je nach Modus
    if(modus == 0) #Sigma = Lambda = I_d
        B = inv(A_bar*transpose(A_bar))
        tr_delta_A_y = transpose(delta_A)*y

        delta_x = transpose(A_bar)*B*(delta_b - delta_A*x)
        delta_y = B*A_bar*(delta_c - tr_delta_A_y)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
    elseif(modus == 1) #Sigma = X^-1, Lambda = S^-1
        Sigma_inv_2 = diagm(x)*diagm(x)
        Lambda_2 = inv(diagm(s)*diagm(s))
        tr_delta_A_y = transpose(delta_A)*y

        delta_x = Sigma_inv_2*transpose(A_bar)*inv(A_bar*Sigma_inv_2*transpose(A_bar))*(delta_b - delta_A*x)
        delta_y = inv(A_bar*Lambda_2*transpose(A_bar))*A_bar*Lambda_2*(delta_c - tr_delta_A_y)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
    elseif(modus == 2) #Sigma = X^-1/2S^1/2, Lambda = X^1/2S^-1/2
        Sigma_inv_2 = diagm(x)*inv(diagm(s))
        tr_delta_A_y = transpose(delta_A)*y
        B = inv(A_bar*Sigma_inv_2*transpose(A_bar))

        delta_x = Sigma_inv_2*transpose(A_bar)*B*(delta_b - delta_A*x)
        delta_y = B*A_bar*Sigma_inv_2*(delta_c - tr_delta_A_y)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
    elseif(modus == 3)#Sigma Lambda wie in meinem Vorschlag
        Sigma = sigmalambda(x, n)
        Lambda = sigmalambda(s, n)
        Sigma_inv_2 = Sigma*inv(Lambda)
        tr_delta_A_y = transpose(delta_A)*y
        B = inv(A_bar*Sigma_inv_2*transpose(A_bar))

        delta_x = Sigma_inv_2*transpose(A_bar)*B*(delta_b - delta_A*x)
        delta_y = B*A_bar*Sigma_inv_2*(delta_c - tr_delta_A_y)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
    elseif(modus == 4) #Sigma Lambda wie in meinem Vorschlag
        Sigma = sigmalambda(x, n)
        Lambda = sigmalambda(s, n)
        Sigma_inv_2 = Sigma*Sigma
        Lambda_2 = inv(Lambda*Lambda)
        tr_delta_A_y = transpose(delta_A)*y

        delta_x = Sigma_inv_2*transpose(A_bar)*inv(A_bar*Sigma_inv_2*transpose(A_bar))*(delta_b - delta_A*x)
        delta_y = inv(A_bar*Lambda_2*transpose(A_bar))*A_bar*Lambda_2*(delta_c - tr_delta_A_y)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
    #Newtonschritt-Adjustment
    else
        D = Sigma_inv_2 = diagm(x)*inv(diagm(s))
        tr_delta_A_y = transpose(delta_A)*y

        delta_y = inv(A_bar*D*transpose(A_bar))*(A_bar*D*(delta_c - tr_delta_A_y) + delta_b - delta_A*x)
        delta_s = delta_c - tr_delta_A_y - transpose(A_bar)*delta_y
        delta_x = -D*delta_s
    end
    w_x = x + delta_x
    w_y = y + delta_y
    w_s = s + delta_s
    #Umformen der Ausgabe
    w_x = svec_to_jumpvec(w_x, n)
    w_s = svec_to_jumpvec(w_s, n)
    return w_x, w_y, w_s
end

function sigmalambda(x, n) #Eingabe: sVec(X); Ausagbe: (Sigma)^(-1) bzw. (Lambda)^(-1) wie in BA
    k = Int(n*(n+1)/2)
    y = zeros(k)
    for l_1 = 1:n
        d_l_1 = d(l_1, n)
        y[d_l_1] = x[d_l_1]
        for l_2 = 1:(n-l_1)
            d_l_2 = d(l_1 + l_2, n)
            y[d_l_1 + l_2] = minimum([x[d_l_1], x[d_l_2]]) 
        end
    end
    return diagm(y)
end

function d(l, n) #Eingabe: l; Ausgabe: Index des lten Diagonalelements
    d = 1
    if(l == 1)
        return d
    end
    for k = 0:(l-2)
        d = d + (n - k)
    end
    return d
end

function ws2(x, y, s, n, lambda, modus = 0) #Eingabe: Optimallösung von model_1; Ausagbe: Warmstartpunkt für model_2 nach Methode 2
    e = e_to_jumpvec(n)
    w_x = lambda*x + (1 - lambda)*e
    w_y = lambda*y
    if(modus == 0)
        w_s = lambda*s + (1 - lambda)*e
    else
        w_s = e
    end 
    return w_x, w_y, w_s
end

function e_to_jumpvec(n) #Eingabe: n; Ausgabe: Darstellung von sMat(e) wie in JuMP
    k = Int(n*(n+1)/2)
    r2 = 2^0.5
    y = 1/r2 * ones(k)
    j = 0
    for i = 1:n
         j = j + i
        y[j] = 1
    end
    return y
end

function jumpvec_to_svec(x, n) #Eingabe: Darstellung von X wie in JuMP; Ausgabe: sVec(X)
    y = Float64[]
    r2 = 2^0.5
    for i = 1:n
        s = 0
        for j = i:n
            k = Int(i*(i+1)/2)
            if(j == i)
                push!(y, x[k+s])
            else
                push!(y, x[k+s]*r2)
            end
            s += j
        end
    end
    return y
end

function svec_to_jumpvec(x, n) #Eingabe: sVec(X); Ausgabe: Darstellung von X wie in JuMP
    y = Float64[]
    r2 = 2^0.5
    for i = 1:n
        s = 0
        for j = 1:i
            if(j == i)
                push!(y, x[i+s])
            else
                push!(y, x[i+s]/r2)
            end
            s += n-j
        end
    end
    return y
end

function matrix_to_svec(c, n) #Eingabe: Matrix C; Ausgabe: sVec(C)
    y = Float64[]
    r2 = 2^0.5
    for i = 1:n
        for j = i:n
            if(j == i)
                push!(y, c[i, j])
            else
                push!(y, c[i, j]*r2)
            end
        end
    end
    return y
end

function tensor_to_matrix(A, n, m) #Eingabe: Tensor A in IR^(mxnxn); Ausgabe: Matrix in IR^(mxd) Zeile i ist sVec(A_i) (siehe BA)
    B = zeros((m, Int(n*(n+1)/2)))
    for i = 1:m
        B[i, :] = matrix_to_svec(A[i, :, :], n)
    end
    return B
end

function generate_sym_matrix(n) #Eingabe: n; Ausgabe: zufällige symmetrische Matrix mit Dimension nxn
    y = Symmetric(2*rand(n, n) - ones(n, n))
    return y
end

function define_model(A, b, C, m, n) #Eingabe: A, b, C; Ausgabe: Model, Nebenbedingungen
    model = Model(SCS.Optimizer)
    @variable(model, x[1:n, 1:n], PSD)
    con = []
    for i = 1:m
        push!(con, @constraint(model, LinearAlgebra.dot(x, A[i, :, :]) == b[i]))
    end
    @objective(model, Min, LinearAlgebra.dot(x, C))
    return model, con
end

function extract_vars(model, con, m) #Eingabe: Model; Ausgabe: Optimallösung
    x = Float64[]
    for x_i in all_variables(model)
        push!(x, value(x_i))
    end
    y = Float64[]
    for i = 1:m
        push!(y, dual(con[i]))
    end
    s = dual(all_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)[1])
    return x, y, s
end

function warmstart(x, y, s, model, con, m) #Eingabe: Warmstartpunkt und Model; Übergabe des Warmstartpunktes an das Model
    #x
    x_2 = all_variables(model) #Liste aller primalen Variablen ([x[1, 1], x[2, 1], ..., x[n, n]])
    for i = 1:length(x_2)
        set_start_value(x_2[i],x[i])
    end
    #y
    for i = 1:m
        set_dual_start_value(con[i], y[i])
    end
    #s
    con_x = all_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle)[1]
    set_dual_start_value(con_x, s)
end

N = [100, 200]
Alphas = [0.01, 0.1, 1]
Iter = [100, 50]

file_names = ["times_c_100_001.csv", "times_c_100_010.csv", "times_c_100_100.csv", "times_c_200_001.csv", "times_c_200_010.csv", "times_c_200_100.csv"]
for (k_1, n) in enumerate(N)
    m = n
    for (k_2, alpha) in enumerate(Alphas)
        iter = Iter[k_1]
        ts = zeros(iter, 19)
        for k_3 = 1:iter
            Random.seed!(k_3)
            A_1 = zeros((m, n, n))
            for i = 1:m
                B = zeros(n, n)
                B[i, i] = 1
                A_1[i, :, :] = copy(B)
            end
            b_1 = ones(m)
            C_1 = generate_sym_matrix(n)

            A_2 = zeros((m, n, n))
            for i = 1:m
                A_2[i, :, :] = copy(A_1[i, :, :])
            end
            b_2 = copy(b_1)
            C_2 = copy(C_1) + generate_sym_matrix(n)*alpha


            #Definieren der Modelle
            model_1, con_1 = define_model(A_1, b_1, C_1, m, n)
            model_2, con_2 = define_model(A_2, b_2, C_2, m, n)
                
            #Optimieren der Modelle (Kaltstart)
            optimize!(model_1)
            t_c = @elapsed begin
                optimize!(model_2)
            end
            ts[k_3, 1] = t_c
            
            #Optimalen Punkt des ersten Problems extrahieren (x, y, s)
            x, y, s = extract_vars(model_1, con_1, m)

            #Berechnung der Warmstartpunkte
            for i = 0:5
                ws1_x, ws1_y, ws1_s = ws1(A_2-A_1, b_2-b_1, C_2-C_1, A_2, x, y, s, n, m, i)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws1_x, ws1_y, ws1_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 2] = t_w
            end

            for i = 0:5
                ws2_x, ws2_y, ws2_s = ws2(x, y, s, n, 0.2*i, 0)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws2_x, ws2_y, ws2_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 8] = t_w
            end

            for i = 0:5
                ws2_x, ws2_y, ws2_s = ws2(x, y, s, n, 0.2*i, 1)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws2_x, ws2_y, ws2_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 14] = t_w
            end

            CSV.write(file_names[3*(k_1-1) + k_2], DataFrame(ts, :auto), header = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])
        end
    end
end

file_names = ["times_b_100.csv", "times_b_200.csv", "times_bc_100.csv", "times_bc_200.csv", "times_abc_100.csv", "times_abc_200.csv"]
for k_0 = 1:3 
    for (k_1, n) in enumerate(N)
        m = n
        eps = 0.01
        iter = Iter[k_1]
        ts = zeros(iter, 19)
        for k_3 = 1:iter
            Random.seed!(k_3)
            A_1 = zeros((m, n, n))
            for i = 1:m
                B = zeros(n, n)
                B[i, i] = 1
                A_1[i, :, :] = copy(B)
            end
            b_1 = ones(m)
            C_1 = generate_sym_matrix(n)
    
            A_2 = zeros((m, n, n))
            for i = 1:m
                if(k_0 == 3)
                    A_2[i, :, :] = copy(A_1[i, :, :]) + generate_sym_matrix(n)*eps
                else
                    A_2[i, :, :] = copy(A_1[i, :, :])
                end
            end

            b_2 = copy(b_1) + eps*(2*rand(m) - ones(m))

            if(k_0 == 2 || k_0 == 3)
                C_2 = copy(C_1) + generate_sym_matrix(n)*eps
            else
                C_2 = copy(C_1)
            end
    
            #Definieren der Modelle
            model_1, con_1 = define_model(A_1, b_1, C_1, m, n)
            model_2, con_2 = define_model(A_2, b_2, C_2, m, n)
                    
            #Optimieren der Modelle (Kaltstart)
            optimize!(model_1)
            t_c = @elapsed begin
                optimize!(model_2)
            end
            ts[k_3, 1] = t_c
                
            #Optimalen Punkt des ersten Problems extrahieren (x, y, s)
            x, y, s = extract_vars(model_1, con_1, m)
   
            #Berechnung der Warmstartpunkte
            for i = 0:5
                ws1_x, ws1_y, ws1_s = ws1(A_2-A_1, b_2-b_1, C_2-C_1, A_2, x, y, s, n, m, i)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws1_x, ws1_y, ws1_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 2] = t_w
            end
    
            for i = 0:5
                ws2_x, ws2_y, ws2_s = ws2(x, y, s, n, 0.2*i, 0)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws2_x, ws2_y, ws2_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 8] = t_w
            end

            for i = 0:5
                ws2_x, ws2_y, ws2_s = ws2(x, y, s, n, 0.2*i, 1)
                #Übergeben des Warmstartpunktes an das zweite Problem
                warmstart(ws2_x, ws2_y, ws2_s, model_2, con_2, m)
                #Optimieren des zweiten Modells (Warmstart)
                t_w = @elapsed begin
                    optimize!(model_2)
                end
                ts[k_3, i + 14] = t_w
            end
            CSV.write(file_names[2*(k_0-1) + k_1], DataFrame(ts, :auto), header = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])
        end
    end
end