## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))


using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV,BlackBoxOptim
using Revise
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
# pyplot()
# plotlyjs()
gr()




## Load data and analysis pipeline
@load("data/case_data_with_pos_neg_21feb_to_30sept.jld2")
@load("data/serologydata_21feb_6aug.jld2")

#Zero fill serodata to be same length as case data
zerofill = zeros(Int64,size(case_data_with_pos_neg.cases,1)-size(sero_data.serodata,1),30,2)
zerofilledserodata = vcat(sero_data.serodata,zerofill)
filleddates = vcat(sero_data.dates,[sero_data.dates[end] + Day(k+1) for k = 1:size(zerofill,1)])

sero_data = (serodata =zerofilledserodata, areas = sero_data.areas,dates = filleddates)

@load("data/relative_testing_rate01102020.jld2")
@load("data/projected_contact_data_09102020.jld2")
@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")

trans_neg_PCR = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,χ = asℝ₊,M_PCR = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))
dir_name = "modelfits_pos_neg_opt"
model_paths = [ joinpath(dir_name,filename) for filename in readdir(dir_name)]

model_nairobi = load("modelfits_pos_neg_opt/Nairobi_model_opt.jld2")["model_opt"]

ct = KenyaSerology.optimise_ct_model(model_nairobi,projected_contactrate_nairobi,6.8,0.2,0.0)

plot(ct)
function run_optimisations(model_paths,projected_contactrate,trans;num_steps = 3,varname = "model")
        for filename in model_paths
                model_opt = KenyaSerology.EM_optimise(filename,projected_contactrate,trans,6.8,0.2,0.05;num_steps=num_steps,varname=varname)
                countyname = model_opt.areaname
                @save("modelfits_pos_neg_opt_vs2/$(countyname)_model_opt.jld2",model_opt)
        end
end

nai_path = "modelfits_pos_neg_opt/Nairobi_model_opt.jld2"
mom_path = "modelfits_pos_neg_opt/Mombasa_model_opt.jld2"
other_paths = setdiff(model_paths,[nai_path,mom_path])

#Optimise on Nairobi
run_optimisations([nai_path],projected_contactrate_nairobi,trans_neg_PCR;num_steps = 0,varname = "model_opt")
#Optimise on Mombasa
run_optimisations([mom_path],projected_contactrate_mombasa,trans_neg_PCR;num_steps = 0,varname = "model_opt")


#Optimise on group 1
run_optimisations(other_paths[1:23],projected_contactrate_kenya,trans_neg_PCR;num_steps = 0,varname = "model_opt")
#Optimise on group 2
run_optimisations(other_paths[23:end],projected_contactrate_kenya,trans_neg_PCR;num_steps = 0,varname = "model_opt")
