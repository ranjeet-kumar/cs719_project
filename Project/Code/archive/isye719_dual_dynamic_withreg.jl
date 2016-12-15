using JuMP
using Gurobi
using PyPlot

############# Funtions to convert JuMP returned variables to arrays ################
function convertToArray(x)
    y = getvalue(x)
    n = length(y)
    a = zeros(n)
    for i = 1:n
	a[i] = y[i]
    end
    return a
end

function convertToArray2(A,n)
    AA = getvalue(A)
    m = (n[1],n[2])
    B = zeros(m)
    for i in 1:n[1]
	for j in 1:n[2]
	    B[i,j] = AA[i,j]
	end
    end
    return B
end

function convertToArray3(A,n)
    AA = getvalue(A)
    m = (n[1],n[2],n[3])
    B = zeros(m)
    for i in 1:n[1]
	for j in 1:n[2]
	    for k in 1:n[3]
		B[i,j,k] = AA[i,j,k]
	    end
	end
    end
    return B
end

#############################################################################################################

loaddata = readcsv("DATARCoast.csv");
rtmpricedata = readcsv("AggregatedData_RTM_ALTA31GT_7_B1.csv")
dampricedata = readcsv("AggregatedData_DAM_ALTA31GT_7_B1.csv")
rtmepricedata = rtmpricedata[:,4];
damepricedata = dampricedata[:,4];
damreguppricedata = dampricedata[:,5];
damregdownpricedata = dampricedata[:,6];


dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day ahead market (hourly) intervals in a day
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtdam/dtrtm);		#No. of real time intervals in each quarter-hour [hour]

nrtmpoints = ndam*nrtm;		        #Total number of points in RTM data
ndampoints = ndam;			#Total number of points in hourly data

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
regup_max = 0.5*P_max;    #Regulation Up Capacity, MW
regdown_max = 0.5*P_max;  #Regulation Up Capacity, MW
rampmin = -100;	          #Lower bound for ramp discharge, MW/s
rampmax = 100;  	  #Upper bound for ramp discharge, MW/s
eff = 1;                #Discharging Efficiency of battery
socend = 100;		  #State of charge at the end of the day
ndays = 365;              #Number of days data is available for
ndays_planning = 7;       #Number of days you want to plan for
ndays_horizon = 1;        #Number of days in the horizon at every step

NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;


load = Matrix{Float64}(loaddata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW
load = reshape(load,nrtm,ndam,ndays);
loadvec = vec(load);

# Reshaping the data as weekly profiles
reshape_ndays = 7;
reshape_nrows = reshape_ndays*nrtmpoints;
reshape_ncolumns = Int64(floor(length(load)/reshape_nrows));
load_estimationdata = loadvec[1:reshape_nrows*reshape_ncolumns];
load_weekly = reshape(load_estimationdata,reshape_nrows,reshape_ncolumns);


#=
# Using Ledoit-Wolfe Sample Covariance Estimator
(p,n) = size(load_weekly);
load_weeklymean = mean(load_weekly,2);
X = load_weekly-repmat(load_weeklymean,1,n); # Columns are 168*1 random vectors with mean 0 and covariance Sigma
Sn = X*X'/n;
mn = trace(Sn*eye(size(Sn)[1])')/p;
dn2 = trace((Sn-mn*eye(size(Sn)[1]))*(Sn-mn*eye(size(Sn)[1]))')/p;
bnbar2 = 0;
for k = 1:n
    bnbar2 = bnbar2 + trace((X[:,k]*X[:,k]'-Sn)*(X[:,k]*X[:,k]'-Sn)')/p;
end
bnbar2 = bnbar2/n^2;
bn2 = min(bnbar2,dn2);
an2 = dn2 - bn2;
Snstar = bn2/dn2*eye(p) + an2/dn2*Sn; # Estimator of variance Sigma


# Generating NS scenarios for weekly load profiles in kW
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));
loadNSdata = zeros(reshape_nrows*(nweeks_planning+1),NS);
R = chol(Snstar);
for j in 1:nweeks_planning+1
    loadNSdata[(j-1)*reshape_nrows+(1:reshape_nrows),:]  = repmat(load_weeklymean,1,NS) + R'*randn(p,NS);
end
loadNSdata[loadNSdata.<=0] = minimum(loadvec);
writecsv("loads_scenarios_month.csv",loadNSdata)


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

######################################################

load1 = loadNSplanningdata/1000;   #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);

=#


# Loading the NS scenarios for weekly load profiles in kW generated from the fullproblem_stochastic code
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));

loadNSdata = readcsv("loads_scenarios_month.csv")


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

################################################################


load1 = loadNSplanningdata/1000; #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);


nhours_planning = ndays_planning*ndam;
nrtm_planning = nhours_planning*nrtm;
nhours_horizon = ndays_horizon*ndam;
nrtm_horizon = nhours_horizon*nrtm;


#=
realized_sequence = rand(S,nhours_planning);
writecsv("realized_sequence.csv",realized_sequence)
=#

# realized_sequence = readcsv("realized_sequence.csv")
# realized_sequence = Vector{Int64}(ones(nhours_planning));



######################################################
# TRYING DUAL DYNAMIC PROGRAMMING

edamprofitlimit = P_max*(sum(damepricedata[find(damepricedata.>0)])-sum(damepricedata[find(damepricedata.<0)]));
ertmprofitlimit = P_max*(sum(rtmepricedata[find(rtmepricedata.>0)])-sum(rtmepricedata[find(rtmepricedata.<0)]));
regupprofitlimit = regup_max*(sum(damreguppricedata[find(damreguppricedata.>0)])-sum(damreguppricedata[find(damreguppricedata.<0)]));
regdownprofitlimit = regdown_max*(sum(damregdownpricedata[find(damregdownpricedata.>0)])-sum(damregdownpricedata[find(damregdownpricedata.<0)]));
thetalimit = -edamprofitlimit - ertmprofitlimit - regupprofitlimit - regdownprofitlimit;
function forwardmodel(ebat0,L)
    m = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m, -P_max <= Prtm[rtm] <= P_max)	                #Net Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam <= P_max)    	                #Net Power sold to the day-ahead market, kW
    @expression(m, Pnet[i in rtm], Prtm[i] + Pdam)              #Net power discharged from battery in all 5-min interval, kW
    @variable(m, 0 <= ebat[rtm] <= ebat_max)                	#Energy stored in the battery at the end of each real time interval
    @variable(m, 0 <= regupdam <= regup_max)                       #Amount of regulation up, kW
    @variable(m, 0 <= regdowndam <= regdown_max)                   #Amount of regulation down, kW
    @variable(m, suppliedload[rtm] >= 0)
    @variable(m, unmetload[rtm] >= 0)
    @variable(m, profitErtm[rtm])# >= 0)				#Profit from the real time market, USD
    @variable(m, profitEdam)# >= 0)	        		#Profit from the day ahead market, USD
    @variable(m, profitregupdam)# >= 0)			        #Profit from the day ahead market, USD
    @variable(m, profitregdowndam)# >= 0)	        		#Profit from the day ahead market, USD
    @variable(m, profittotal)# >= 0)		                	#Total profit in the day, USD
    @variable(m, unmetcost)
    @variable(m, theta >= thetalimit)
    @constraint(m, InitialEnergy, ebat[1] == ebat0 - 1/eff*Pnet[1]*dtrtm - suppliedload[1]*dtrtm)	#Inital energy in the battery

    #    @constraint(m, EndSOC[i=rtm[end]], soc[i] >= socend)		#Constraint on SOC at the end of the day

    @constraint(m, rtmEBalance[i in rtm[2:end]], ebat[i] == ebat[i-1] - 1/eff*Pnet[i]*dtrtm - suppliedload[i]*dtrtm)	#Dynamics constraint

    # @constraint(m, RTMRamp[i in rtm[2:end],k in dam, s in S], rampmin*dtrtm <= Pnet[i,k,s]  - Pnet[i-1,k,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

    # @constraint(m, DAMRamp[i=rtm[1],k in dam[2:end],iend=rtm[end],s in S], rampmin*dtrtm <= Pnet[i,k,s] - Pnet[iend,k-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

    @constraint(m, RegUp[i in rtm], Pnet[i] + regupdam <= P_max)	#Constraint on total power

    @constraint(m, RegDown[i in rtm], Pnet[i] - regdowndam >= -P_max)	#Constraint on total power

    @constraint(m, UnmetLoad[i in rtm], suppliedload[i] + unmetload[i] >=  L[i][1])#load[i,s_realized][1])

    @constraint(m, BoundSupplied[i in rtm], suppliedload[i] <= L[i][1])

    @constraint(m, BoundUnmet[i in rtm], unmetload[i] <= L[i][1])

    @constraint(m, RTMEProfits[i in rtm], profitErtm[i] == rtmepr[i]*Prtm[i]*dtrtm)	#Economic calculation
    @constraint(m, DAMEProfits, profitEdam == damepr*Pdam*dtdam)        	#Economic calculation

    @constraint(m, DAMregupProfits, profitregupdam == damreguppr*regupdam)
    @constraint(m, DAMregdownProfits, profitregdowndam == damregdownpr*regdowndam)

    @constraint(m, TotalProfit, profittotal ==
                        sum{profitErtm[i], i in rtm} + profitEdam
                        + profitregupdam + profitregdowndam)

    @constraint(m, UnmetCost, unmetcost == sum{rtmepr[i]*unmetload[i], i in rtm})

    @objective(m, Min, -profittotal + unmetcost + theta)

    return m;

end



function scenariosubproblem(ebat0,s)
    m = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m, -P_max <= Prtm[rtm] <= P_max)	                #Net Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam <= P_max)    	                #Net Power sold to the day-ahead market, kW
    @expression(m, Pnet[i in rtm], Prtm[i] + Pdam)              #Net power discharged from battery in all 5-min interval, kW
    @variable(m, 0 <= ebat[rtm] <= ebat_max)                	#Energy stored in the battery at the end of each real time interval
    @variable(m, 0 <= regupdam <= regup_max)                       #Amount of regulation up, kW
    @variable(m, 0 <= regdowndam <= regdown_max)                   #Amount of regulation down, kW
    @variable(m, suppliedload[rtm] >= 0)
    @variable(m, unmetload[rtm] >= 0)
    @variable(m, profitErtm[rtm])# >= 0)				#Profit from the real time market, USD
    @variable(m, profitEdam)# >= 0)	        		#Profit from the day ahead market, USD
    @variable(m, profitregupdam)# >= 0)			        #Profit from the day ahead market, USD
    @variable(m, profitregdowndam)# >= 0)	        		#Profit from the day ahead market, USD
    @variable(m, profittotal)# >= 0)		                	#Total profit in the day, USD
    @variable(m, unmetcost)
    @variable(m, theta >= thetalimit)
    @constraint(m, InitialEnergy, ebat[1] == ebat0 - 1/eff*Pnet[1]*dtrtm - suppliedload[1]*dtrtm)	#Inital energy in the battery

    #    @constraint(m, EndSOC[i=rtm[end]], soc[i] >= socend)		#Constraint on SOC at the end of the day

    @constraint(m, rtmEBalance[i in rtm[2:end]], ebat[i] == ebat[i-1] - 1/eff*Pnet[i]*dtrtm - suppliedload[i]*dtrtm)	#Dynamics constraint

    # @constraint(m, RTMRamp[i in rtm[2:end],k in dam, s in S], rampmin*dtrtm <= Pnet[i,k,s]  - Pnet[i-1,k,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

   # @constraint(m, DAMRamp[i=rtm[1],k in dam[2:end],iend=rtm[end],s in S], rampmin*dtrtm <= Pnet[i,k,s] - Pnet[iend,k-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

    @constraint(m, RegUp[i in rtm], Pnet[i] + regupdam <= P_max)	#Constraint on total power

    @constraint(m, RegDown[i in rtm], Pnet[i] - regdowndam >= -P_max)	#Constraint on total power

    @constraint(m, UnmetLoad[i in rtm], suppliedload[i] + unmetload[i] >=  load[i,s][1])

    @constraint(m, BoundSupplied[i in rtm], suppliedload[i] <= load[i,s][1])

    @constraint(m, BoundUnmet[i in rtm], unmetload[i] <= load[i,s][1])

    @constraint(m, RTMEProfits[i in rtm], profitErtm[i] == rtmepr[i]*Prtm[i]*dtrtm)	#Economic calculation
    @constraint(m, DAMEProfits, profitEdam == damepr*Pdam*dtdam)        	#Economic calculation

    @constraint(m, DAMregupProfits, profitregupdam == damreguppr*regupdam)
    @constraint(m, DAMregdownProfits, profitregdowndam == damregdownpr*regdowndam)

    @constraint(m, TotalProfit, profittotal ==
                        sum{profitErtm[i], i in rtm} + profitEdam
                        + profitregupdam + profitregdowndam)

    @constraint(m, UnmetCost, unmetcost == sum{rtmepr[i]*unmetload[i], i in rtm})


    @objective(m, Min, -profittotal + unmetcost + theta)

    return m;

end








#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
dam = 1:nhours_horizon;

rtmepr = Vector(nrtm);	    	#Real Time Market price, $/MWh
damepr = nothing;
damreguppr = nothing;
damregdownpr = nothing;


ebat0 = ebat_max;		  #Initial State of charge, 100 means fully charged
m_f = Array{JuMP.Model}(nhours_planning)
m_b = Array{JuMP.Model}(NS,nhours_planning)


dual = Vector(NS);
obj = Vector(NS)
v = Vector(NS)
lowerbound = 0; upperbound = 1e10;
policy_cost = Vector();

# Solving root node to bigin forward pass for first iteration
node0problem = Model(solver = GurobiSolver(OutputFlag=0,Threads=2));
@variable(node0problem, theta >= thetalimit)
@objective(node0problem, Min, theta)
solve(node0problem)
lowerbound = getobjectivevalue(node0problem);

ebat0 = ebat_max
tic()
# First iteration of Forward and backward passes
j=1;
println("Iteration Number : $(j)")
# Forward pass for first iteration starting
for p in 1:nhours_planning # Forward pass
    #Load and price data
    load = loaddata1[(p-1)*nrtm+(1:nrtm),:];	#Load, MW
    rtmepr = rtmpricedata[(p-1)*nrtm+(1:nrtm),4];	    	#Real Time Market price, $/MWh
    damepr = dampricedata[p,4];	    	#Day Ahead Market Selling price, $/MWh
    damreguppr = dampricedata[p,5];	    	#Day Ahead Market Regulation up price, $/MWh
    damregdownpr = dampricedata[p,6]; 	#Day Ahead Market Regulation down price, $/MWh
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
# Forward pass for first iteration ended and got upperbound from first iteration
# Backward pass starting for first iteration
v_avg = 0; pi_avg = 0;
ebat0_b = nothing;
for p in nhours_planning:-1:1 # Backward pass
  load = loaddata1[(p-1)*nrtm+(1:nrtm),:];	#Load, MW
  rtmepr = rtmpricedata[(p-1)*nrtm+(1:nrtm),4];	    	#Real Time Market price, $/MWh
  damepr = dampricedata[p,4];	    	#Day Ahead Market Selling price, $/MWh
  damreguppr = dampricedata[p,5];	    	#Day Ahead Market Regulation up price, $/MWh
  damregdownpr = dampricedata[p,6]; 	#Day Ahead Market Regulation down price, $/MWh
  if p ==1
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
end # Backward pass

theta = getvariable(node0problem,:theta);
@constraint(node0problem, ThetaConst, theta >= v_avg + pi_avg*ebat_max)
status0 = solve(node0problem);
lowerbound = getobjectivevalue(node0problem);
# Backward pass for first iteration ended at root node
# Iterations 2,3,... start
while abs(upperbound - lowerbound) >= 0.1 # Iterations
j = j+1;
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
# Backward pass for iteration j ended at root node
println("Current upperbound = $(upperbound), Current lowerbound = $(lowerbound)")
end # Iteration j ends

time_taken_dual_dynamic = toc();



































#= Start comment here


println("Total Profits: ", sum(dailyprofit),"\n")



totalpowerplot = reshape(totalpower,nrtmpoints*ndays);
upperbandplot = reshape(upperband,nrtmpoints*ndays);
lowerbandplot = reshape(lowerband,nrtmpoints*ndays);





#################### Plotting #####################

timertm = 0:dtrtm:ndam*ndays;
timedam = 0:dtdam:ndam*ndays;

loadplot = push!(loadplot,loadplot[end]);
rtmeprplot = push!(rtmeprplot,rtmeprplot[end]);
dameprplot = push!(dameprplot,dameprplot[end]);
damregupprplot = push!(damregupprplot,damregupprplot[end]);
damregdownprplot = push!(damregdownprplot,damregdownprplot[end]);
#ebatplot = push!(ebatplot,ebatplot[end]);
ebatplot = [ebat_max,ebatplot];
#socplot = push!(socplot,socplot[end]);
socplot = [soc0,socplot];
Prtmplot = push!(Prtmplot,Prtmplot[end]);
Pdamplot = push!(Pdamplot,Pdamplot[end]);
regupdamplot = push!(regupdamplot,regupdamplot[end]);
regdowndamplot = push!(regdowndamplot,regdowndamplot[end]);
rtmeprofitplot = push!(rtmeprofitplot,rtmeprofitplot[end]);
dameprofitplot = push!(dameprofitplot,dameprofitplot[end]);
damregprofitplot = push!(damregprofitplot,damregprofitplot[end]);



totalpowerplot = push!(totalpowerplot,totalpowerplot[end]);
upperbandplot = push!(upperbandplot,upperbandplot[end]);
lowerbandplot = push!(lowerbandplot,lowerbandplot[end]);




function cumul(A)
    C = zeros(length(A)-1);
    for i in 1:length(A)-1
        C[i] = sum(A[1:i]);
    end
    C = push!(C,C[end]);
    return C;
end

function damtortm(A)
    A = A[1:ndampoints*ndays]
    B = zeros(nrtmpoints*ndays);
    for i in 1:ndampoints*ndays
        B[(i-1)*12+1:i*12] = repmat([A[i]],12);
    end
    B = push!(B,B[end])
    return B
end

rtmeprofitplot = cumul(rtmeprofitplot);
dameprofitplot = cumul(dameprofitplot);
totaleprofitplot = rtmeprofitplot + damtortm(dameprofitplot);
damregprofitplot = cumul(damregprofitplot);
totalregprofitplot = damtortm(damregprofitplot);
totalprofitplot = totaleprofitplot + totalregprofitplot;


figure()
plot(timertm,loadplot,drawstyle="steps-post",label="Load (MW)",LineWidth=1.5)
#hold(true)
#plot(timertm,[loadmeanplot,loadmeanplot[end]])
xlabel("Time (hr)",size=20)
ylabel("Load (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timertm[(1:nrtmpoints+1)],rtmeprplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="RTM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("RTM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("RTM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timertm[(1:nrtmpoints+1)],Prtmplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam[(1:ndampoints+1)],dameprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="DAM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam[(1:ndampoints+1)],Pdamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam[(1:ndampoints+1)],damregupprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timedam[(1:ndampoints+1)],damregdownprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam[(1:ndampoints+1)],regupdamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in DAM",LineWidth=1.5)
hold(true)
plot(timedam[(1:ndampoints+1)],regdowndamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in DAM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm[(1:nrtmpoints+1)],socplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",label="SOC",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("SOC (%)",size=20)
grid()
tick_params(labelsize=18)
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")


figure()
plot(timertm[(1:nrtmpoints+1)],totalpowerplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="Total Power (MW)",LineWidth=1.5)
hold(true)
plot(timertm[(1:nrtmpoints+1)],lowerbandplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="red",label="Lower band (MW)",LineWidth=1.5)
plot(timertm[(1:nrtmpoints+1)],upperbandplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="green",label="Upper band (MW)",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("Power (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")




figure()
plot(timertm/24,rtmeprofitplot,drawstyle="steps-post",color="blue",label="RTM",LineWidth=1.5)
xlabel("Time (day)",size=20), ylabel("Revenues from Energy (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,dameprofitplot,drawstyle="steps-post",color="cyan",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totaleprofitplot,drawstyle="steps-post",color="red",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
xlabel("Time (day)",size=20), ylabel("Revenues from Reg (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,damregprofitplot,drawstyle="steps-post",color="cyan",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totalregprofitplot,drawstyle="steps-post",color="red",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm/24,totalprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
xlabel("Time (day)",size=14), ylabel("Total Revenue (\$)",size=20)
grid();tick_params(labelsize=18);#xlim(0,370)



#End comment here
=#







#=







figure()
plot(timertm,loadplot,drawstyle="steps-post",label="Load (MW)",LineWidth=1.5)
#hold(true)
#plot(timertm,[loadmeanplot,loadmeanplot[end]])
xlabel("Time (hr)",size=20)
ylabel("Load (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timertm,rtmeprplot,drawstyle="steps-post",color="blue",label="RTM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("RTM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("RTM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timertm,Prtmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam,dameprplot,drawstyle="steps-post",color="blue",label="DAM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam,edamplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam,damregupprplot,drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timedam,damregdownprplot,drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam,regupdamplot,drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in DAM",LineWidth=1.5)
hold(true)
plot(timedam,regdowndamplot,drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in DAM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm,socplot,drawstyle="steps-post",label="SOC",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("SOC (%)",size=20)
grid()
tick_params(labelsize=18)
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")



figure()
plot(timertm/24,rtmeprofitplot,drawstyle="steps-post",label="RTM",LineWidth=1.5)
xlabel("Time (day)",size=20), ylabel("Revenues from Energy (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,dameprofitplot,drawstyle="steps-post",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totaleprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
xlabel("Time (day)",size=20), ylabel("Revenues from Reg (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,damregprofitplot,drawstyle="steps-post",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totalregprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm/24,totalprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
xlabel("Time (day)",size=14), ylabel("Total Revenue (\$)",size=20)
grid();tick_params(labelsize=18);#xlim(0,370)



function addBasic


m = Model()

addBasicConstraints(m)

















=#
