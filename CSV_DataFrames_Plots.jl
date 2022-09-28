using DataFrames
using CSV
using Plots
plotlyjs()

df = DataFrame(CSV.File("new_S210_220923_A00690_0267_BHL7VNDRX2_demux_report.csv", header = 16))

#remove the rows with "missing"
df = df[completecases(df), :]
#vscodedisplay(df)

df[:, :x] = parse.(Float64, df[:, :"%_of_the_lane"])
df[:, :y] = parse.(Float64, df[:, :"%_Perfect_barcode"])
df[:, :z] = df[:, :"%>=Q30bases"]

scatter(df.x, df.y, group = df.iLabID, legend=:bottomright)
scatter(df.x, df.y, group = df.Lane, legend=:bottomright)

scatter(df.x, df.z, group = df.iLabID, legend=:top)
scatter(df.x, df.z, group = df.Lane, legend=:top)