using CSV, DataFrames

function f(t, r, len) #Ausgabe: prozentualer Anteil, welcher <= r ist
    sum = 0
    for k = 1:len
        if(t[k] <= r)
            sum += 1
        end
    end
    return (sum/len)*100
end


#Ratio zu bestem Startpunkt
file_names = ["data/times_c_100_001.csv", "data/times_c_100_010.csv", "data/times_c_100_100.csv", "data/times_c_200_001.csv", "data/times_c_200_010.csv", "data/times_c_200_100.csv", "data/times_b_100.csv", "data/times_b_200.csv", "data/times_bc_100.csv", "data/times_bc_200.csv", "data/times_abc_100.csv", "data/times_abc_200.csv"]
ana_file_names = ["analysed_data/ratio_c_100_001.csv", "analysed_data/ratio_c_100_010.csv", "analysed_data/ratio_c_100_100.csv", "analysed_data/ratio_c_200_001.csv", "analysed_data/ratio_c_200_010.csv", "analysed_data/ratio_c_200_100.csv", "analysed_data/ratio_b_100.csv", "analysed_data/ratio_b_200.csv", "analysed_data/ratio_bc_100.csv", "analysed_data/ratio_bc_200.csv", "analysed_data/ratio_abc_100.csv", "analysed_data/ratio_abc_200.csv"]
plot_file_names = ["analysed_data/ratio_plot_c_100_001.csv", "analysed_data/ratio_plot_c_100_010.csv", "analysed_data/ratio_plot_c_100_100.csv", "analysed_data/ratio_plot_c_200_001.csv", "analysed_data/ratio_plot_c_200_010.csv", "analysed_data/ratio_plot_c_200_100.csv", "analysed_data/ratio_plot_b_100.csv", "analysed_data/ratio_plot_b_200.csv", "analysed_data/ratio_plot_bc_100.csv", "analysed_data/ratio_plot_bc_200.csv", "analysed_data/ratio_plot_abc_100.csv", "analysed_data/ratio_plot_abc_200.csv"]

N = 200 #Auflösung des Plots
ratios = range(1, 1.6, length=N)
for (k, file_name) in enumerate(file_names)
    file = CSV.read(file_name, DataFrame)
    I = length(file[:, 1])
    J = length(file[1, :])
    for i = 1:I
        min_i = minimum(file[i, :])
        for j = 1:J
            file[i, j] = file[i, j] / min_i
        end
    end
    CSV.write(ana_file_names[k], file, header = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])

    data = zeros(N, J + 1)
    for i = 1:N
        data[i, 1] = ratios[i]
    end
    for n = 1:N 
        for j = 1:J
            data[n, j + 1] = f(file[:, j], ratios[n], I)
        end
    end
    CSV.write(plot_file_names[k], DataFrame(data, :auto), header = ["x", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])
end


#Ratio zu Kaltstart n = 100
file_names = ["data/times_c_100_001.csv", "data/times_c_100_010.csv", "data/times_c_100_100.csv", "data/times_b_100.csv", "data/times_bc_100.csv", "data/times_abc_100.csv"]
ana_file_names = ["analysed_data/ratio_cold_c_100_001.csv", "analysed_data/ratio_cold_c_100_010.csv", "analysed_data/ratio_cold_c_100_100.csv", "analysed_data/ratio_cold_b_100.csv", "analysed_data/ratio_cold_bc_100.csv", "analysed_data/ratio_cold_abc_100.csv"]

ratios = [0, 0.5, 1, 1.5]
data = zeros(18*4, 6)
for (k, file_name) in enumerate(file_names)
    file = CSV.read(file_name, DataFrame)
    I = length(file[:, 1])
    J = length(file[1, :])
    for i = 1:I
        t_cold = file[i, 1]
        for j = 1:J
            file[i, j] = file[i, j] / t_cold
        end
    end
    CSV.write(ana_file_names[k], file, header = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])

    for n = 1:4
        for j = 2:J
            if (n <= 3)
                data[(j-2)*4 + n, k] = f(file[:, j], ratios[n+1], I) - f(file[:, j], ratios[n], I)
            else
                data[(j-2)*4 + n, k] = 100 - data[(j-2)*4 + 1, k] - data[(j-2)*4 + 2, k] - data[(j-2)*4 + 3, k]
            end
        end
    end
end
CSV.write("analysed_data/ratios_100.csv", DataFrame(data, :auto), header = ["c001", "c010", "c100", "b001", "bc001", "abc001"])


#Ratio zu Kaltstart n = 200
file_names = ["data/times_c_200_001.csv", "data/times_c_200_010.csv", "data/times_c_200_100.csv", "data/times_b_200.csv", "data/times_bc_200.csv", "data/times_abc_200.csv"]
ana_file_names = ["analysed_data/ratio_cold_c_200_001.csv", "analysed_data/ratio_cold_c_200_010.csv", "analysed_data/ratio_cold_c_200_200.csv", "analysed_data/ratio_cold_b_200.csv", "analysed_data/ratio_cold_bc_200.csv", "analysed_data/ratio_cold_abc_200.csv"]

ratios = [0, 0.5, 1, 1.5]
data = zeros(18*4, 6)
for (k, file_name) in enumerate(file_names)
    file = CSV.read(file_name, DataFrame)
    I = length(file[:, 1])
    J = length(file[1, :])
    for i = 1:I
        t_cold = file[i, 1]
        for j = 1:J
            file[i, j] = file[i, j] / t_cold
        end
    end
    CSV.write(ana_file_names[k], file, header = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"])

    for n = 1:4
        for j = 2:J
            if (n <= 3)
                data[(j-2)*4 + n, k] = f(file[:, j], ratios[n+1], I) - f(file[:, j], ratios[n], I)
            else
                data[(j-2)*4 + n, k] = 100 - data[(j-2)*4 + 1, k] - data[(j-2)*4 + 2, k] - data[(j-2)*4 + 3, k]
            end
        end
    end
end
CSV.write("analysed_data/ratios_200.csv", DataFrame(data, :auto), header = ["c001", "c010", "c100", "b001", "bc001", "abc001"])

#Vergleich zwischen w_p und w_dp
file_names_100 = ["data/times_c_100_001.csv", "data/times_c_100_010.csv", "data/times_c_100_100.csv", "data/times_b_100.csv", "data/times_bc_100.csv", "data/times_abc_100.csv"]
file_names_200 = ["data/times_c_200_001.csv", "data/times_c_200_010.csv", "data/times_c_200_100.csv", "data/times_b_200.csv", "data/times_bc_200.csv", "data/times_abc_200.csv"]

sum_100 = [0.0, 0.0, 0.0, 0.0, 0.0]
for (k, file_name) in enumerate(file_names_100)
    file = CSV.read(file_name, DataFrame)
    I = length(file[:, 1])
    J = length(file[1, :])
    for i = 1:I
        for j = 1:5
            sum_100[j] += file[i, j+8+6]/file[i, j+8]
        end
    end
end
for j = 1:5
    sum_100[j] /= 600
end
sum_200 = [0.0, 0.0, 0.0, 0.0, 0.0]
for (k, file_name) in enumerate(file_names_200)
    file = CSV.read(file_name, DataFrame)
    I = length(file[:, 1])
    J = length(file[1, :])
    for i = 1:I
        for j = 1:5
            sum_200[j] += file[i, j+8+6]/file[i, j+8]
        end
    end
end
for j = 1:5
    sum_200[j] /= 300
end

println((sum_100 + sum_200)/2)