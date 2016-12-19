using JuMP
using PyPlot
using Gurobi
using StatsFuns

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
ndays_horizon = 1;    # Number of days in horizon at every step of receding horizon scheme
nhours_horizon = ndays_horizon*ndam;  # Number of hoours in horizon at every step of receding horizon scheme
nrtm_horizon = nhours_horizon*nrtm;   # Number of real-time intervals in horizon at every step of receding horizon scheme

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
ebat0 = ebat_max;		   #Initial State of charge
rampmin = -0.5*P_max;	          #Lower bound for ramp discharge, MW/5min
rampmax = 0.5*P_max;  	  #Upper bound for ramp discharge, MW/5min
NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;

# Price data
eprrtm = rtmpriceExceldata[1:nrtm_planning,4];	 #Real Time Market price, $/MWh
eprdam = dampriceExceldata[1:nhours_planning,4];	 #Day Ahead Market price, $/MWh

# Generate new scenarios for loads
loadNSdatafile = "loads_scenarios_month.csv";
if generate_new_scenarios == 1
  load = Matrix{Float64}(loadExceldata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW
  loadNSdata = generate_weekly_loads_scenarios(load,1:52,ndays_planning,NS,loadNSdatafile);
end
loadNSdata = readcsv(loadNSdatafile);
ndays_data = (nweeks_planning+1)*weekly_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW
load1 = loadNSplanningdata/1000; #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);


######################################################
# STOCHASTIC DUAL DYNAMIC PROGRAMMING

edamprofitlimit = P_max*(sum(eprdam[find(eprdam.>0)])-sum(eprdam[find(eprdam.<0)]));
ertmprofitlimit = P_max*(sum(eprdam[find(eprdam.>0)])-sum(eprdam[find(eprdam.<0)]));
thetalimit = -edamprofitlimit - ertmprofitlimit;
function forwardmodel(ebat0,L)
    m = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m, -P_max <= Prtm[rtm] <= P_max)	                #Net Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam <= P_max)    	                #Net Power sold to the day-ahead market, kW
    @variable(m, 0 <= ebat[rtm] <= ebat_max)                	#Energy stored in the battery at the end of each real time interval
    @variable(m, suppliedload[rtm] >= 0)
    @variable(m, unmetload[rtm] >= 0)
    @expression(m, Pnet[i in rtm], Prtm[i] + Pdam + suppliedload[i])              #Net power discharged from battery in all 5-min interval, kW
    @variable(m, profitErtm[rtm])				#Profit from the real time market, USD
    @variable(m, profitEdam)	        		#Profit from the day ahead market, USD
    @variable(m, profittotal)		                	#Total profit in the day, USD
    @variable(m, unmetcost)
    @variable(m, theta >= thetalimit)
    @constraint(m, InitialEnergy, ebat[1] == ebat0 - Pnet[1]*dtrtm)	#Inital energy in the battery
    @constraint(m, rtmEBalance[i in rtm[2:end]], ebat[i] == ebat[i-1] - Pnet[i]*dtrtm)	#Dynamics constraint
    @constraint(m, NetDischarge1[i in rtm], Pnet[i] <= P_max)
    @constraint(m, NetDischarge2[i in rtm], Pnet[i] >= -P_max)
    @constraint(m, UnmetLoad[i in rtm], suppliedload[i] + unmetload[i] >=  L[i][1])
    @constraint(m, BoundSupplied[i in rtm], suppliedload[i] <= L[i][1])
    @constraint(m, BoundUnmet[i in rtm], unmetload[i] <= L[i][1])
    @constraint(m, RTMEProfits[i in rtm], profitErtm[i] == rtmepr[i]*Prtm[i]*dtrtm)	#Economic calculation
    @constraint(m, DAMEProfits, profitEdam == damepr*Pdam*dtdam)        	#Economic calculation
    @constraint(m, TotalProfit, profittotal == sum{profitErtm[i], i in rtm} + profitEdam)
    @constraint(m, UnmetCost, unmetcost == sum{rtmepr[i]*unmetload[i], i in rtm})
    @objective(m, Min, -profittotal + unmetcost + theta)
    return m;
end

function scenariosubproblem(ebat0,s)
    m = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m, -P_max <= Prtm[rtm] <= P_max)	                #Net Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam <= P_max)    	                #Net Power sold to the day-ahead market, kW
    @variable(m, 0 <= ebat[rtm] <= ebat_max)                	#Energy stored in the battery at the end of each real time interval
    @variable(m, suppliedload[rtm] >= 0)
    @variable(m, unmetload[rtm] >= 0)
    @expression(m, Pnet[i in rtm], Prtm[i] + Pdam + suppliedload[i])              #Net power discharged from battery in all 5-min interval, kW
    @variable(m, profitErtm[rtm])				#Profit from the real time market, USD
    @variable(m, profitEdam)	        		#Profit from the day ahead market, USD
    @variable(m, profittotal)		                	#Total profit in the day, USD
    @variable(m, unmetcost)
    @variable(m, theta >= thetalimit)
    @constraint(m, InitialEnergy, ebat[1] == ebat0 - Pnet[1]*dtrtm)	#Inital energy in the battery
    @constraint(m, rtmEBalance[i in rtm[2:end]], ebat[i] == ebat[i-1] - Pnet[i]*dtrtm)	#Dynamics constraint
    @constraint(m, NetDischarge1[i in rtm], Pnet[i] <= P_max)
    @constraint(m, NetDischarge2[i in rtm], Pnet[i] >= -P_max)
    @constraint(m, UnmetLoad[i in rtm], suppliedload[i] + unmetload[i] >=  load[i,s][1])
    @constraint(m, BoundSupplied[i in rtm], suppliedload[i] <= load[i,s][1])
    @constraint(m, BoundUnmet[i in rtm], unmetload[i] <= load[i,s][1])
    @constraint(m, RTMEProfits[i in rtm], profitErtm[i] == rtmepr[i]*Prtm[i]*dtrtm)	#Economic calculation
    @constraint(m, DAMEProfits, profitEdam == damepr*Pdam*dtdam)        	#Economic calculation
    @constraint(m, TotalProfit, profittotal == sum{profitErtm[i], i in rtm} + profitEdam)
    @constraint(m, UnmetCost, unmetcost == sum{rtmepr[i]*unmetload[i], i in rtm})
    @objective(m, Min, -profittotal + unmetcost + theta)
    return m;
end

#Define sets to be used in the model
rtm = 1:nrtm;
dam = 1:ndam;
rtmepr = Vector(nrtm);	    	#Real Time Market price, $/MWh
damepr = nothing;
load = nothing;


lowerbound_Vector = Vector{Float64}();
upperbound_Vector = Vector{Float64}();
upperbound_up_Vector = Vector{Float64}();
upperbound_down_Vector = Vector{Float64}();


m_f = Array{JuMP.Model}(nhours_planning);
m_b = Array{JuMP.Model}(NS,nhours_planning);
dual = Vector(NS);
obj = Vector(NS);
v = Vector(NS);
lowerbound = -1e10; upperbound = 1e10;
policy_cost = Vector();

# Stochastic Dual Dynamic Programming Algorithm Starts
tic()
# Solving root node to begin forward pass for first iteration
ebat0 = ebat_max;		  #Initial State of charge, 100 means fully charged
node0problem = Model(solver = GurobiSolver(OutputFlag=0,Threads=2));
@variable(node0problem, theta >= thetalimit)
@objective(node0problem, Min, theta)
solve(node0problem)
lowerbound = getobjectivevalue(node0problem);
push!(lowerbound_Vector,lowerbound)

ebat0 = ebat_max;
# First iteration of Forward and backward passes
j=1; # Iteration number
println("Iteration Number : $(j)")
# Forward pass for first iteration starting
for p in 1:nhours_planning # Forward pass
    #Load and price data
    load = loaddata1[(p-1)*nrtm+(1:nrtm),:];	#Load, MW
    rtmepr = eprrtm[(p-1)*nrtm+(1:nrtm)];	    	#Real Time Market price, $/MWh
    damepr = eprdam[p];	    	#Day Ahead Market Selling price, $/MWh
    s_realized = rand(S);
    m_f[p] = forwardmodel(ebat0,load[:,s_realized]);
    status = solve(m_f[p]);
    ebat0 = getvalue(getvariable(m_f[p],:ebat))[rtm[end]];
end # End forward pass

cost_hourly = Vector(nhours_planning);
for i in 1:nhours_planning
  cost_hourly[i] = getobjectivevalue(m_f[i])-getvalue(getvariable(m_f[i],:theta));
end
push!(policy_cost, sum(cost_hourly));
upperbound = mean(policy_cost);
println("Current upperbound = $(upperbound), Current lowerbound = $(lowerbound)")
push!(upperbound_Vector,upperbound);
upperbound_up = upperbound;
upperbound_down = upperbound;
push!(upperbound_up_Vector,upperbound_up);
push!(upperbound_down_Vector,upperbound_down);

# Forward pass for first iteration ended and got upperbound from first iteration
# Backward pass starting for first iteration
v_avg = 0; pi_avg = 0;
ebat0_b = nothing;
for p in nhours_planning:-1:1 # Backward pass
  load = loaddata1[(p-1)*nrtm+(1:nrtm),:];	#Load, MW
  rtmepr = eprrtm[(p-1)*nrtm+(1:nrtm)];	    	#Real Time Market price, $/MWh
  damepr = eprdam[p];	    	#Day Ahead Market Selling price, $/MWh
  if p == 1
    ebat0_b = ebat_max;
  else
    ebat0_b = getvalue(getvariable(m_f[p-1],:ebat))[rtm[end]];
  end
    for s in S
        m_b[s,p] = scenariosubproblem(ebat0_b,s);
        theta = getvariable(m_b[s,p],:theta)
        ebat = getvariable(m_b[s,p],:ebat)
        @constraint(m_b[s,p], ThetaConst, theta >= v_avg + pi_avg*ebat[rtm[end]])
        status = solve(m_b[s,p]);
        dual[s] = getdual(getconstraint(m_b[s,p],:InitialEnergy))+sum(getdual(getconstraint(m_b[s,p],:rtmEBalance)[:]));
        obj[s] = getobjectivevalue(m_b[s,p]);
        v[s] = obj[s] - dual[s]*ebat0_b;
    end
    v_avg = (1/NS)*(sum(v));
    pi_avg = (1/NS)*(sum(dual));
end # Backward pass end

theta = getvariable(node0problem,:theta);
@constraint(node0problem, ThetaConst, theta >= v_avg + pi_avg*ebat_max)
status0 = solve(node0problem);
lowerbound = getobjectivevalue(node0problem);
push!(lowerbound_Vector,lowerbound)
# Backward pass for first iteration ended at root node
# Iterations 2,3,... start
while abs(upperbound - lowerbound) >= 0.1 # Iterations
j = j+1; # Iteration number
println("Iteration Number : $(j)")
# Root node for this forward pass was already solved at the end of backward pass of previous iteration
ebat0 = ebat_max;
s_realized_for = Vector()
for p in 1:nhours_planning # Forward pass
    s_realized = rand(S);
    push!(s_realized_for,s_realized)
    JuMP.setRHS(getconstraint(m_b[s_realized,p],:InitialEnergy), ebat0)
    status = solve(m_b[s_realized,p]);
    ebat0 = getvalue(getvariable(m_b[s_realized,p],:ebat))[rtm[end]];
end # End forward pass

cost_hourly = Vector(nhours_planning);
for i in 1:nhours_planning
  cost_hourly[i] = getobjectivevalue(m_b[s_realized_for[i],i])-getvalue(getvariable(m_b[s_realized_for[i],i],:theta));
end
push!(policy_cost, sum(cost_hourly));
upperbound = mean(policy_cost);
push!(upperbound_Vector,upperbound);
# 95% 2-sided confidence interval on upperbound
upperbound_up = upperbound + tdistinvcdf(length(upperbound_Vector)-1, 0.975)*std(upperbound_Vector)/sqrt(length(upperbound_Vector));
upperbound_down = upperbound -  tdistinvcdf(length(upperbound_Vector)-1, 0.975)*std(upperbound_Vector)/sqrt(length(upperbound_Vector));;
push!(upperbound_up_Vector,upperbound_up);
push!(upperbound_down_Vector,upperbound_down);

# Forward pass for iteration j ended and got upperbound from iteration j
# Backward pass starting for iteration j
v_avg = 0; pi_avg = 0;
ebat0_b = nothing;
for p in nhours_planning:-1:1 # Backward pass
  if p ==1
    ebat0_b = ebat_max;
  else
    ebat0_b = getvalue(getvariable(m_b[s_realized_for[p-1],p-1],:ebat))[rtm[end]];
  end
    for s in S
        JuMP.setRHS(getconstraint(m_b[s,p],:InitialEnergy), ebat0_b)
        theta = getvariable(m_b[s,p],:theta)
        ebat = getvariable(m_b[s,p],:ebat)
        @constraint(m_b[s,p], ThetaConst, theta >= v_avg + pi_avg*ebat[rtm[end]])
        status = solve(m_b[s,p]);
        dual[s] = getdual(getconstraint(m_b[s,p],:InitialEnergy))+sum(getdual(getconstraint(m_b[s,p],:rtmEBalance)[:]));
        obj[s] = getobjectivevalue(m_b[s,p]);
        v[s] = obj[s] - dual[s]*ebat0_b;
    end
    v_avg = (1/NS)*(sum(v));
    pi_avg = (1/NS)*(sum(dual));
end # End backward pass
theta = getvariable(node0problem,:theta);
@constraint(node0problem, ThetaConst, theta >= v_avg + pi_avg*ebat_max)
status0 = solve(node0problem);
lowerbound = getobjectivevalue(node0problem);
push!(lowerbound_Vector,lowerbound)
# Backward pass for iteration j ended at root node
println("Current upperbound = $(upperbound), Current lowerbound = $(lowerbound)")
end # Iteration j ends
time_taken_dual_dynamic = toc();

if makeplots == 1
  xplot = 0:length(lowerbound_Vector)-1;
  figure()
  plt[:get_current_fig_manager]()[:full_screen_toggle]()
  hold(true)
  plot(xplot,lowerbound_Vector, color = "blue", label="Lower bound",LineWidth=2);
  plot(xplot[1:end-1],upperbound_Vector, color = "red", label="Upper bound",LineWidth=2);
  plot(xplot[1:end-1],upperbound_up_Vector, color = "grey",label="Confidence interval on upper bound",LineWidth=2);
  plot(xplot[1:end-1],upperbound_down_Vector, color = "grey",LineWidth=2);
  grid()
  ylabel("Expected cost",size = 24)
  xlabel("Iteration",size = 24)
  xlim(0,220)
  ylim(-565,-540)
  tick_params(labelsize=20)
  legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
  savefig(string("cs719figures/dual/dual_dynamic_bounds.pdf"))
  close("all")
  figure()
  xplot = 0:length(lowerbound_Vector)-1;
  figure()
  plt[:get_current_fig_manager]()[:full_screen_toggle]()
  hold(true)
  plot(xplot[5:end],lowerbound_Vector[5:end], color = "blue", label="Lower bound");
  plot(xplot[6:end],upperbound_Vector[5:end], color = "red", label="Upper bound");
  plot(xplot[6:end],upperbound_up_Vector[5:end], color = "grey",label="Confidence interval on upper bound");
  plot(xplot[6:end],upperbound_down_Vector[5:end], color = "grey");
  grid()
  xlim(5,xplot[end])
  ylabel("Upper & lower bounds",size = 24)
  xlabel("Iteration",size = 24)
  tick_params(labelsize=20)
  legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
  savefig(string("cs719figures/dual/dual_dynamic_bounds_zoomed.pdf"))
  close("all")

end # End makeplots
