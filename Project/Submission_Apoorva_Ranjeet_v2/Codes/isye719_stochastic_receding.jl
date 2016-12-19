using JuMP
using PyPlot
using Gurobi

include("helper_functions.jl")
#############################################################################################################
# Read load and price data from csv files
loadExceldata = readcsv("data_loads_past_1year.csv");
rtmpriceExceldata = readcsv("data_RTM_prices_2015.csv")
dampriceExceldata = readcsv("data_DAM_prices_2015.csv")

generate_new_scenarios = 0; # 0: don't generate new scenarios for load profiles; 1: generate new scenarios
generate_new_sample_paths = 0; # 0: don't generate new sample paths for loads; 1: generate new sample paths
participate_rtm = 1; # 0: Don't participate in RTM; 1: participate in RTM
participate_dam = 1; # 0: Don't participate in DAM; 1: participate in DAM
makeplots = 0; # 0: Don't generate plots, 1: Generate plots

# Defining time parameters for the model and planning schedule
dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day-ahead market (hourly) intervals in a day
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtdam/dtrtm);		#No. of real time intervals in each hour
nrtmpoints = ndam*nrtm;		    #Total number of points in RTM data in one day
ndampoints = ndam;			  #Total number of points in hourly data in one day
weekly_ndays = 7;         #Number of days in every week
ndays = 365;              #Number of days data is available for
ndays_planning = 7;       #Number of days you want to plan for
nweeks_planning = Int64(ceil((ndays_planning/weekly_ndays)));
nhours_planning = ndays_planning*ndam;  #Number of hours we will plan the policy for
nrtm_planning = nhours_planning*nrtm;   #Number of rtm intervals we will plan the policy for
ndays_horizon = 1; # Number of days in horizon at every step of receding horizon scheme
nhours_horizon = ndays_horizon*ndam; # Number of hours in horizon at every step of receding horizon scheme
nrtm_horizon = nhours_horizon*nrtm; # Number of real-time intervals in horizon at every step of receding horizon scheme

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
ebat0 = ebat_max;		   #Initial State of charge
rampmin = -0.5*P_max;	          #Lower bound for ramp discharge, MW/5min
rampmax = 0.5*P_max;  	  #Upper bound for ramp discharge, MW/5min
NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;

# Price data
eprrtm = rtmpriceExceldata[1:nrtmpoints*ndays,4];	 #Real Time Market price, $/MWh
eprdam = dampriceExceldata[1:ndampoints*ndays,4];	 #Day Ahead Market price, $/MWh

# Generate new scenarios for loads
loadNSdatafile = "loads_scenarios_month.csv";
if generate_new_scenarios == 1
  load = Matrix{Float64}(loadExceldata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW
  loadNSdata = generate_weekly_loads_scenarios(load,1:52,ndays_planning,NS,loadNSdatafile);
end
loadNSdata = readcsv(loadNSdatafile)
ndays_data = (nweeks_planning+1)*weekly_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

#Reshape the data to matrices to be used in the model
rtmepr = reshape(eprrtm,nrtm,ndam,ndays);
damepr = reshape(eprdam,ndam,ndays);
load = loadNSplanningdata[:,:,1:ndays_planning,:]/1000; #MW

#Define sets to be used in the JuMP model
rtm = 1:nrtm; # {1,2,...,12}
dam = 1:ndam; # {1,2,...,24}
day = 1:ndays_planning; # {1,2,...,7}

if generate_new_sample_paths == 1
# Generate NS sample paths for realizations of loads in 7 days at hourly intervals
  (paths,loadperm) = generate_sample_paths(load,NS,"samplepaths.csv","sampleloadperm.csv");
end
# Take the NS sample paths for loads generated earlier
paths = Matrix{Int64}(readcsv("samplepaths.csv"));
loadperm = zeros(nrtm,ndam,ndays_planning,NS);
for s in S
  j = 1;
  for l in day
    for k in dam
        loadperm[:,k,l,s] = load[:,k,l,paths[j,s]];  #MW
        j = j+1;
    end
  end
end


load1 = loadNSplanningdata/1000; #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);





dailyprofit = zeros(ndays);
loadplot = zeros(ndays*nrtmpoints);
rtmeprplot = zeros(ndays*nrtmpoints);
dameprplot = zeros(ndays*ndampoints);
ebatplot = zeros(ndays*nrtmpoints);
socplot = zeros(ndays*nrtmpoints);
Prtmplot = zeros(ndays*nrtmpoints);
Pdamplot = zeros(ndays*ndampoints);
rtmeprofitplot = zeros(ndays*nrtmpoints);
dameprofitplot = zeros(ndays*ndampoints);

totalpower = zeros(nrtm,ndam,ndays);
totalpowerplot = zeros(nrtm,ndam,ndays);


Prtm_realized = zeros(nrtm,nhours_planning);
unmetload_realized = zeros(nrtm,nhours_planning)
unmetcost_realized = zeros(nrtm,nhours_planning);
profitErtm_realized = zeros(nrtm,nhours_planning);
profitEdam_realized = zeros(nhours_planning);
profitE_realized = zeros(nhours_planning);
profittotal_realized = zeros(nhours_planning);
netobjective_realized = zeros(nhours_planning);


m_rol = nothing;

obj_st_rh_NS = Vector()

for k in S # Loop to evaluate cost along each scenario

ebat0 = ebat_max;		  #Initial State of charge, 100 means fully charged


#realized_sequence = Vector{Int64}(k*ones(nhours_planning));
realized_sequence = paths[:,k];

Prtm_realized = zeros(nrtm,nhours_planning);
unmetload_realized = zeros(nrtm,nhours_planning)
unmetcost_realized = zeros(nrtm,nhours_planning);
profitErtm_realized = zeros(nrtm,nhours_planning);
profitEdam_realized = zeros(nhours_planning);
profitE_realized = zeros(nhours_planning);
profittotal_realized = zeros(nhours_planning);
netobjective_realized = zeros(nhours_planning);


j=1;
for p in 1:nhours_planning


    #Load and price data
    load = loaddata1[(p-1)*nrtm+(1:nrtm_horizon),:];	#Load, MW

    eprrtm = rtmpriceExceldata[(p-1)*nrtm+(1:nrtm_horizon),4];	    	#Real Time Market price, $/MWh
    eprdam = dampriceExceldata[(p-1)+(1:nhours_horizon),4];	    	#Day Ahead Market Selling price, $/MWh

    #Reshape the data to matrices
    rtmepr = reshape(eprrtm,nrtm,nhours_horizon);
    damepr = reshape(eprdam,nhours_horizon);
    load = reshape(load,nrtm,nhours_horizon,NS);

    #Define sets to be used in the model defVar and addConstraint
    rtm = 1:nrtm;
    dam = 1:nhours_horizon;

    ################ Model ##################

    m_rol = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m_rol, -P_max <= Prtm[rtm,dam,S] <= P_max)	                #Power sold to the real time market, kW
    @variable(m_rol, -P_max <= Pdam[dam,S] <= P_max)    	                #Power sold to the day ahead market, kW
    @variable(m_rol, 0 <= ebat[rtm,dam,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
    @variable(m_rol, suppliedload[rtm,dam,S] >= 0)
    @variable(m_rol, unmetload[rtm,dam,S] >= 0)
		@expression(m_rol, Pnet[i in rtm,k in dam,s in S], Prtm[i,k,s] + Pdam[k,s] + suppliedload[i,k,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m_rol, profitErtm[rtm,dam,S])				        #Profit from the real time market, USD
    @variable(m_rol, profitEdam[dam,S])	        			#Profit from the day ahead market, USD
    @variable(m_rol, profittotal[S])		        	#Total profit in the day, USD
    @variable(m_rol, unmetcost[rtm,dam,S])

    @constraint(m_rol, InitialEnergy[s in S], ebat[1,1,s] == ebat0 - Pnet[1,1,s]*dtrtm)	#Inital energy in the battery
    @constraint(m_rol, rtmEBalance[i in rtm[2:end],k in dam,s in S], ebat[i,k,s] == ebat[i-1,k,s] - Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rol, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],s in S], ebat[i,k,s] == ebat[iend,k-1,s] - Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rol, UnmetLoad[i in rtm,k in dam,s in S], suppliedload[i,k,s] + unmetload[i,k,s] >=  load[i,k,s])
    @constraint(m_rol, BoundSupplied[i in rtm,k in dam,s in S], suppliedload[i,k,s] <= load[i,k,s])
    @constraint(m_rol, BoundUnmet[i in rtm,k in dam,s in S], unmetload[i,k,s] <= load[i,k,s])
		@constraint(m_rol, RTMEProfits[i in rtm,k in dam,s in S], profitErtm[i,k,s] == rtmepr[i,k]*Prtm[i,k,s]*dtrtm)	#Economic calculation
    @constraint(m_rol, DAMEProfits[k in dam,s in S], profitEdam[k,s] == damepr[k]*Pdam[k,s]*dtdam)	#Economic calculation
    @constraint(m_rol, TotalProfit[s in S], profittotal[s] ==
                        sum{profitErtm[i,k,s], i in rtm, k in dam} + sum{profitEdam[k,s], k in dam})
    @constraint(m_rol, UnmetCost[i in rtm, k in dam, s in S], unmetcost[i,k,s] == rtmepr[i,k]*unmetload[i,k,s])
    @constraint(m_rol, NetDischarge1[i in rtm,k in dam,s in S], Pnet[i,k,s] <= P_max)
    @constraint(m_rol, NetDischarge2[i in rtm,k in dam,s in S], Pnet[i,k,s] >= -P_max)
    # Non-anticipativity constraints for first stage variables
    @constraint(m_rol, Nonant_PDAM[k in dam,s in S], Pdam[k,s] == (1/NS)*sum{Pdam[k,s], s in S})
    @objective(m_rol, Min, (1/NS)*sum{-profittotal[s] + sum{unmetcost[i,k,s],i in rtm, k in dam}, s in S})
    tic()
    status = solve(m_rol)

##########################################################################

    ebat0 = getvalue(getvariable(m_rol,:ebat))[rtm[end],dam[1],realized_sequence[p]];
    Prtm_realized[:,p] = getvalue(getvariable(m_rol,:Prtm))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetload_realized[:,p] = getvalue(getvariable(m_rol,:unmetload))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetcost_realized[:,p] = getvalue(getvariable(m_rol,:unmetcost))[1:rtm[end],dam[1],realized_sequence[p]];
    profitErtm_realized[:,p] = getvalue(getvariable(m_rol,:profitErtm))[1:rtm[end],dam[1],realized_sequence[p]];
    profitEdam_realized[p] = getvalue(getvariable(m_rol,:profitEdam))[dam[1],realized_sequence[p]];
    profittotal_realized[p] = sum(profitErtm_realized[:,p]) + profitEdam_realized[p];
    netobjective_realized[p] = sum(unmetcost_realized[:,p]) - profittotal_realized[p];
    time_taken_st_rolling = toc();
    println("Scenario $k, Step $p, netobjective_realized = $(netobjective_realized[p]), time = $(time_taken_st_rolling)")

j = j+1;

end # End rolling horizon


totalcost_after_rolling_st = sum(netobjective_realized);

obj_st_rh = totalcost_after_rolling_st;

push!(obj_st_rh_NS,obj_st_rh)

end


expected_obj_st_rh = mean(obj_st_rh_NS);
