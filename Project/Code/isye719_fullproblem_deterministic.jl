
using JuMP
using PyPlot
using Gurobi

include("helper_functions.jl")
#############################################################################################################
# Read load and price data from csv files
loadExceldata = readcsv("DATARCoast.csv");
rtmpriceExceldata = readcsv("AggregatedData_RTM_ALTA31GT_7_B1.csv")
dampriceExceldata = readcsv("AggregatedData_DAM_ALTA31GT_7_B1.csv")

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
        loadperm[:,k,l,s] = load[:,k,l,paths[j,s]];
        j = j+1;
    end
  end
end
loadperm_mv = mean(loadperm,4); # Mean loads for mean-value problem

################ Mean Value Model ##################
tic()
mv = Model(solver = GurobiSolver(Threads=2))
		@variable(mv, -P_max <= Prtm[rtm,dam,day] <= P_max)	                #Power sold to the real time market, kW
		@variable(mv, -P_max <= Pdam[dam,day] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(mv, 0 <= ebat[rtm,dam,day] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(mv, suppliedload[rtm,dam,day] >= 0)
		@variable(mv, unmetload[rtm,dam,day] >= 0)
		@expression(mv, Pnet[i in rtm,k in dam,l in day], Prtm[i,k,l] + Pdam[k,l] + suppliedload[i,k,l])    #Net power discharged from battery in all 5-min interval, kW
		@variable(mv, profitErtm[rtm,dam,day])				        #Profit from the real time market, USD
		@variable(mv, profitEdam[dam,day])	        			#Profit from the day ahead market, USD
		@variable(mv, profittotal)	        	#Total profit in the day, USD
		@variable(mv, unmetcost)

		@constraint(mv, InitialEnergy[s in S], ebat[1,1,1] == ebat0 - Pnet[1,1,1]*dtrtm)	#Inital energy in the battery
		@constraint(mv, rtmEBalance[i in rtm[2:end],k in dam,l in day], ebat[i,k,l] == ebat[i-1,k,l] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day], ebat[i,k,l] == ebat[iend,k-1,l] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], ebat[i,k,l] == ebat[iend,kend,l-1] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, UnmetLoad[i in rtm,k in dam,l in day], suppliedload[i,k,l] + unmetload[i,k,l] >=  loadperm_mv[i,k,l])
		@constraint(mv, BoundSupplied[i in rtm,k in dam,l in day], suppliedload[i,k,l] <= loadperm_mv[i,k,l])
		@constraint(mv, BoundUnmet[i in rtm,k in dam,l in day], unmetload[i,k,l] <= loadperm_mv[i,k,l])
    @constraint(mv, NetDischarge1[i in rtm,k in dam,l in day], Pnet[i,k,l] <= P_max)
    @constraint(mv, NetDischarge2[i in rtm,k in dam,l in day], Pnet[i,k,l] >= -P_max)
		#=    @constraint(mv, RTMRamp1[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, RTMRamp2[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(mv, RTMEProfits[i in rtm,k in dam,l in day], profitErtm[i,k,l] == rtmepr[i,k,l]*Prtm[i,k,l]*dtrtm)	#Economic calculation
		@constraint(mv, DAMEProfits[k in dam,l in day], profitEdam[k,l] == damepr[k,l]*Pdam[k,l]*dtdam)	#Economic calculation
		@constraint(mv, TotalProfit, profittotal ==
												sum{profitErtm[i,k,l], i in rtm, k in dam, l in day} + sum{profitEdam[k,l], k in dam, l in day})
		@constraint(mv, UnmetCost, unmetcost == sum{rtmepr[i,k,l]*unmetload[i,k,l], i in rtm, k in dam, l in day})
		@objective(mv, Min, -profittotal + unmetcost)
		if participate_dam == 0
      @constraint(mv, NoDAM[k in dam,l in day, s in S], Pdam[k,l,s] == 0)
    end
    if participate_rtm == 0
      @constraint(mv, NoRTM[i in rtm,k in dam,l in day, s in S], Prtm[i,k,l,s] == 0)
    end
status = solve(mv);
Pdam_mv = getvalue(getvariable(mv,:Pdam));

######################### Mean value Model End ###########################

################ Evaluating costs in all scenarios with mv-solution ##################

m_d = Model(solver = GurobiSolver(Threads=2))
		@variable(m_d, -P_max <= Prtm[rtm,dam,day,S] <= P_max)	                #Power sold to the real time market, kW
		@variable(m_d, -P_max <= Pdam[dam,day,S] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(m_d, 0 <= ebat[rtm,dam,day,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(m_d, suppliedload[rtm,dam,day,S] >= 0)
		@variable(m_d, unmetload[rtm,dam,day,S] >= 0)
		@expression(m_d, Pnet[i in rtm,k in dam,l in day,s in S], Prtm[i,k,l,s] + Pdam[k,l,s] + suppliedload[i,k,l,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m_d, profitErtm[rtm,dam,day,S])				        #Profit from the real time market, USD
		@variable(m_d, profitEdam[dam,day,S])	        			#Profit from the day ahead market, USD
		@variable(m_d, profittotal[S])		        	#Total profit in the day, USD
		@variable(m_d, unmetcost[S])

		@constraint(m_d, InitialEnergy[s in S], ebat[1,1,1,s] == ebat0 - Pnet[1,1,1,s]*dtrtm)	#Inital energy in the battery
		@constraint(m_d, rtmEBalance[i in rtm[2:end],k in dam,l in day,s in S], ebat[i,k,l,s] == ebat[i-1,k,l,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], ebat[i,k,l,s] == ebat[iend,k-1,l,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], ebat[i,k,l,s] == ebat[iend,kend,l-1,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, UnmetLoad[i in rtm,k in dam,l in day, s in S], suppliedload[i,k,l,s] + unmetload[i,k,l,s] >=  loadperm[i,k,l,s])
		@constraint(m_d, BoundSupplied[i in rtm,k in dam,l in day,s in S], suppliedload[i,k,l,s] <= loadperm[i,k,l,s])
		@constraint(m_d, BoundUnmet[i in rtm,k in dam,l in day,s in S], unmetload[i,k,l,s] <= loadperm[i,k,l,s])
    @constraint(m_d, NetDischarge1[i in rtm,k in dam,l in day,s in S], Pnet[i,k,l,s] <= P_max)
    @constraint(m_d, NetDischarge2[i in rtm,k in dam,l in day,s in S], Pnet[i,k,l,s] >= -P_max)
		#=    @constraint(m_d, RTMRamp1[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, RTMRamp2[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(m_d, RTMEProfits[i in rtm,k in dam,l in day,s in S], profitErtm[i,k,l,s] == rtmepr[i,k,l]*Prtm[i,k,l,s]*dtrtm)	#Economic calculation
		@constraint(m_d, DAMEProfits[k in dam,l in day,s in S], profitEdam[k,l,s] == damepr[k,l]*Pdam[k,l,s]*dtdam)	#Economic calculation
		@constraint(m_d, TotalProfit[s in S], profittotal[s] ==
												sum{profitErtm[i,k,l,s], i in rtm, k in dam, l in day} + sum{profitEdam[k,l,s], k in dam, l in day})
		@constraint(m_d, UnmetCost[s in S], unmetcost[s] == sum{rtmepr[i,k,l]*unmetload[i,k,l,s], i in rtm, k in dam, l in day})
    # Fixing first stage variables at mean-value solution
		@constraint(m_d, Fix_PDAM[k in dam,l in day,s in S], Pdam[k,l,s] == Pdam_mv[k,l])
    @objective(m_d, Min, (1/NS)*sum{-profittotal[s] + unmetcost[s], s in S})
		if participate_dam == 0
      @constraint(m_d, NoDAM[k in dam,l in day, s in S], Pdam[k,l,s] == 0)
    end
    if participate_rtm == 0
      @constraint(m_d, NoRTM[i in rtm,k in dam,l in day, s in S], Prtm[i,k,l,s] == 0)
    end
		status = solve(m_d)

time_taken_dt_fullproblem = toc();
###############################################################

println("\nExpected Objective with Mean Value first stage solution ", getobjectivevalue(m_d),"\n" )

obj_dt_fp = getobjectivevalue(m_d);


if makeplots == 1
################# PLOTTING #################
#=
# Plot of Scenarios of loads
xplot = 0:dtrtm:dtrtm*nrtm_planning
loadNSplot = loadNSdata[1:nrtm_planning,:];
loadNSplot = [loadNSplot;loadNSplot[end,:]];
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S
plot(xplot,loadNSplot[:,s], color="grey", drawstyle="steps-post");
end
plot(xplot,mean(loadNSplot,2), color="blue", drawstyle="steps-post",label="Mean scenario");
grid()
xlim(0,nhours_planning)
ylabel("Loads (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/loads_scenarios.pdf"))
close("all")
=#

r = 5;
# Plot of Pdam
n3 = [ndam,ndays_planning,NS];
eprdamplot = eprdam[1:nhours_planning];
push!(eprdamplot,eprdamplot[end]);
Pdamarray = convertToArray3(getvalue(getvariable(m_d,:Pdam)),n3);
Pdamplot = reshape(Pdamarray[:,:,1],nhours_planning);
push!(Pdamplot,Pdamplot[end]);
xplot = 0:dtdam:dtdam*nhours_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
subplot(2,1,1)
hold(true)
plot(xplot,eprdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Energy price (\$/MWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
plot(xplot,Pdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("DAM Power (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Pdam_fp_dt.pdf"))
close("all")

# Plot of Prtm
n4 = [nrtm,ndam,ndays_planning,NS];
eprrtmplot = eprdam[1:nrtm_planning];
push!(eprrtmplot,eprrtmplot[end]);
Prtmarray = convertToArray4(getvalue(getvariable(m_d,:Prtm)),n4);
Prtmplot = reshape(Prtmarray,nrtm_planning,NS);
Prtmplot = [Prtmplot;Prtmplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
subplot(2,1,1)
hold(true)
plot(xplot,eprrtmplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Energy price (\$/MWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
for s in S[r]
	plot(xplot,Prtmplot[:,s],  drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("RTM Power (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Prtm_fp_dt.pdf"))
close("all")

# Plot of supplied load and unmet load
n4 = [nrtm,ndam,ndays_planning,NS];
suppliedloadarray = convertToArray4(getvalue(getvariable(m_d,:suppliedload)),n4);
suppliedloadplot = reshape(suppliedloadarray,nrtm_planning,NS);
suppliedloadplot = [suppliedloadplot;suppliedloadplot[end,:]];
unmetloadarray = convertToArray4(getvalue(getvariable(m_d,:unmetload)),n4);
unmetloadplot = reshape(unmetloadarray,nrtm_planning,NS);
unmetloadplot = [unmetloadplot;suppliedloadplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,suppliedloadplot[:,s], drawstyle="steps-post", color = "green", label="Supplied");
	plot(xplot,unmetloadplot[:,s], drawstyle="steps-post", color = "red", label="Unmet");
end
grid()
xlim(0,nhours_planning)
ylabel("Supplied & unmet loads (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/supp_unmet_fp_dt.pdf"))
close("all")


# Plot of SOC
n4 = [nrtm,ndam,ndays_planning,NS];
socarray = convertToArray4(getvalue(getvariable(m_d,:ebat)),n4)/ebat_max*100;
socplot = reshape(socarray,nrtm_planning,NS);
socplot = [socplot;socplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,socplot[:,s], drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylim(-0.02,100.02)
ylabel("State of charge (%)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/soc_fp_dt.pdf"))
close("all")

# Pnet, regulation bands calculation
netpower = zeros(nrtm,ndam,ndays_planning,NS);
for i in rtm
        for k in dam
            for l in day
							for s in S
                netpower[i,k,l,s] = Prtmarray[i,k,l,s] + Pdamarray[k,l,s] + suppliedloadarray[i,k,l,s];
							end
            end
        end
end
netpowerplot = reshape(netpower,nrtm_planning,NS);
netpowerplot = [netpowerplot;netpowerplot[end,:]];

# Plot of netpower
n4 = [nrtm,ndam,ndays_planning,NS];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,netpowerplot[:,s], drawstyle="steps-post", color = "blue", label="Net discharge");
end
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("Net discharge (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/netpower_fp_dt.pdf"))
close("all")

# Plot of Profits from Pdam and Prtm and Total
n3 = [ndam,ndays_planning,NS];
profitEdamarray = convertToArray3(getvalue(getvariable(m_d,:profitEdam)),n3);
profitEdamplot = reshape(profitEdamarray[:,:,1],nhours_planning);
cumulprofitEdamplot = cumul(profitEdamplot);
cumulprofitEdamtortmplot = damtortm(cumulprofitEdamplot);

n4 = [nrtm,ndam,ndays_planning,NS];
profitErtmarray = convertToArray4(getvalue(getvariable(m_d,:profitErtm)),n4);
profitErtmplot = reshape(profitErtmarray,nrtm_planning,NS);
profitErtmplotmean = mean(profitErtmplot,2);
cumulprofitErtmplot = zeros(nrtm_planning,NS);
for s in S
  cumulprofitErtmplot[:,s] = cumul(profitErtmplot[:,s]);
end
cumulprofitErtmplotmean = cumul(profitErtmplotmean);


n4 = [nrtm,ndam,ndays_planning,NS];
eprrtmplot = eprrtm[1:nrtm_planning];
unmetloadarray = convertToArray4(getvalue(getvariable(m_d,:unmetload)),n4);
unmetloadplot = reshape(unmetloadarray,nrtm_planning,NS);
costunmetloadplot = zeros(nrtm_planning,NS);
cumulcostunmetloadplot = zeros(nrtm_planning,NS);
for s in S
  costunmetloadplot[:,s] = eprrtmplot.*unmetloadplot[:,s];
  cumulcostunmetloadplot[:,s] = cumul(costunmetloadplot[:,s]);
end
cumulcostunmetloadplotmean = mean(cumulcostunmetloadplot,2);
cumultotalprofitplot = repmat(cumulprofitEdamtortmplot,1,NS) + cumulprofitErtmplot - cumulcostunmetloadplot;
cumultotalprofitplotmean = mean(cumultotalprofitplot,2);

cumulprofitEdamtortmplot = [cumulprofitEdamtortmplot;cumulprofitEdamtortmplot[end]];
cumulprofitErtmplotmean = [cumulprofitErtmplotmean;cumulprofitErtmplotmean[end]];
cumulcostunmetloadplotmean = [cumulcostunmetloadplotmean;cumulcostunmetloadplotmean[end]];
cumultotalprofitplotmean = [cumultotalprofitplotmean;cumultotalprofitplotmean[end]];
cumulprofitErtmplot = [cumulprofitErtmplot;cumulprofitErtmplot[end,:]];
cumulcostunmetloadplot = [cumulcostunmetloadplot;cumulcostunmetloadplot[end,:]];
cumultotalprofitplot = [cumultotalprofitplot;cumultotalprofitplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
plot(xplot,cumulprofitEdamtortmplot, drawstyle="steps-post", color = "blue", label="Day-ahead market");
plot(xplot,cumulprofitErtmplotmean, drawstyle="steps-post", color = "green", label="Real-time market");
plot(xplot,cumulcostunmetloadplotmean, drawstyle="steps-post", color = "red", label="Unmet load");
plot(xplot,cumultotalprofitplotmean, drawstyle="steps-post", color = "indigo", label="Total revenue");
grid()
xlim(0,nhours_planning)
#ylim(-1.02,1.02)
ylabel("Expected Cumulative Revenue (\$)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=20)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/cumulative_rev_fp_dt.pdf"))
close("all")

end #END IF Makeplots
