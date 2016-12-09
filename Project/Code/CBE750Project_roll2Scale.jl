
using JuMP
using PyPlot
using Gurobi

close("all")

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
eff = 0.9;                  #Discharging Efficiency of battery
soc0 = 100;		  #Initial State of charge, 100 means fully charged
socend = 100;		  #State of charge at the end of the day

rollh = 1; #Roll the horizon every rollh day
ntotaldays = 30;

m = nothing
day = nothing
dam = nothing
rtm = nothing

for p in 1:ntotaldays

    if p >= 335
        ntotaldays = 365 - p;
        end
    
ndays = ntotaldays - (p-1)*rollh;

dailyprofit = zeros(ndays);
loadplot = zeros(ndays*nrtmpoints);
rtmeprplot = zeros(ndays*nrtmpoints);
dameprplot = zeros(ndays*ndampoints);
damregupprplot = zeros(ndays*ndampoints);
damregdownprplot = zeros(ndays*ndampoints);
ebatplot = zeros(ndays*nrtmpoints);
socplot = zeros(ndays*nrtmpoints);
Prtmplot = zeros(ndays*nrtmpoints);
Pdamplot = zeros(ndays*ndampoints);
regupdamplot = zeros(ndays*ndampoints);
regdowndamplot = zeros(ndays*ndampoints);
dameprofitplot = zeros(ndays*ndampoints);
damregprofitplot = zeros(ndays*ndampoints);

#Load and price data
load = loaddata[2:nrtmpoints+1,2+p-1+(1:ndays)]/1000;	#Load, MW

eprrtm = rtmpricedata[(nrtmpoints*(p-1)+1):nrtmpoints*((p-1)+ndays),4];	    	#Real Time Market price, $/MWh    
eprdam = dampricedata[(ndampoints*(p-1)+1):ndampoints*((p-1)+ndays),4];	    	#Day Ahead Market Selling price, $/MWh    
regupprdam = dampricedata[(ndampoints*(p-1)+1):ndampoints*((p-1)+ndays),5];	    	#Day Ahead Market Regulation up price, $/MWh
regdownprdam = dampricedata[(ndampoints*(p-1)+1):ndampoints*((p-1)+ndays),6]; 		#Day Ahead Market Regulation down price, $/MWh


#Reshape the data to matrices
rtmepr = reshape(eprrtm,nrtm,ndam,ndays);
damepr = reshape(eprdam,ndam,ndays);
damreguppr = reshape(regupprdam,ndam,ndays);
damregdownpr = reshape(regdownprdam,ndam,ndays);
load = reshape(load,nrtm,ndam,ndays);


#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
dam = 1:ndam;
day = 1:ndays;

################ Model ##################

m = Model(solver = GurobiSolver())

    @variable(m, -P_max <= Prtm[rtm,dam,day] <= P_max)	                #Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam[dam,day] <= P_max)    	                #Power sold to the day ahead market, kW
    @expression(m, Pnet[i=rtm,k=dam,l=day], Prtm[i,k,l] + Pdam[k,l])    #Net power discharged from battery in all 5-min interval, kW 
    @variable(m, 0 <= ebat[rtm,dam,day] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval
    @variable(m, 0 <= soc[rtm,dam,day] <= 100)		        #SOC of the battery at the end of each real time interval      
    @variable(m, 0 <= regupdam[dam,day] <= regup_max)                 #Amount of regulation up, kW
    @variable(m, 0 <= regdowndam[dam,day] <= regdown_max)             #Amount of regulation down, kW
    @variable(m, profitErtm[rtm,dam,day])# >= 0)				        #Profit from the real time market, USD
    @variable(m, profitEdam[dam,day])# >= 0)	        			#Profit from the day ahead market, USD		
    @variable(m, profitregupdam[dam,day])# >= 0)			        #Profit from the day ahead market, USD
    @variable(m, profitregdowndam[dam,day])# >= 0)	        		#Profit from the day ahead market, USD
    @variable(m, profittotal)# >= 0)		        	#Total profit in the day, USD		




    @constraint(m, InitialEnergy, ebat[1,1,1] == soc0/100*ebat_max - (eff*Pnet[1,1,1] - load[1,1,1])*dtrtm)	#Inital energy in the battery
        
    @constraint(m, DefSOC[i=rtm,k=dam,l=day], soc[i,k,l] == ebat[i,k,l]/ebat_max*100)			#Define SOC

#    @constraint(m, EndSOC[i=rtm[end],k=dam[end],l=day], soc[i,k,l] >= socend)		#Constraint on SOC at the end of the day
    
    @constraint(m, rtmEBalance[i=rtm[2:end],k=dam,l=day], ebat[i,k,l] == ebat[i-1,k,l] - (eff*Pnet[i,k,l] - load[i,k,l])*dtrtm)	#Dynamics constraint
    
    @constraint(m, damEBalance[i=rtm[1],k=dam[2:end],iend=rtm[end],l=day], ebat[i,k,l] == ebat[iend,k-1,l] - (eff*Pnet[i,k,l] - load[i,k,l])*dtrtm)	#Dynamics constraint

    @constraint(m, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l=day[2:end]], ebat[i,k,l] == ebat[iend,kend,l-1] - (eff*Pnet[i,k,l] - load[i,k,l])*dtrtm)	#Dynamics constraint

   # @constraint(m, RTMRamp[i=rtm[2:end],k=dam,l=day], rampmin*dtrtm <= Pnet[i,k,l]  - Pnet[i-1,k,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time    

    @constraint(m, RegUp[i=rtm,k=dam,l=day], Pnet[i,k,l] + regupdam[k,l] <= P_max)	#Constraint on total power

    @constraint(m, RegDown[i=rtm,k=dam,l=day], Pnet[i,k,l] - regdowndam[k,l] >= -P_max)	#Constraint on total power

  println(Prtm)    
    @constraint(m, RTMEProfits[i=rtm,k=dam,l=day], profitErtm[i,k,l] == rtmepr[i,k,l]*dtrtm*Prtm[i,k,l])	#Economic calculation

    @constraint(m, DAMEProfits[k=dam,l=day], profitEdam[k,l] == damepr[k,l]*Pdam[k,l])	#Economic calculation

    @constraint(m, DAMregupProfits[k=dam,l=day], profitregupdam[k,l] == damreguppr[k,l]*regupdam[k,l])

    @constraint(m, DAMregdownProfits[k=dam,l=day], profitregdowndam[k,l] == damregdownpr[k,l]*regdowndam[k,l])

    @constraint(m, TotalProfit, profittotal == sum(profitErtm) + sum(profitEdam) + sum(profitregupdam) + sum(profitregdowndam))
    
    @objective(m, Max, profittotal);


#    print(m)

    status = solve(m)

    soc0 = getValue(soc[nrtm,ndam,p]);


    println("\nTotal Profits on day ", getValue(profittotal),"\n" )

n3 = [nrtm,ndam,ndays]
n2 = [ndam,ndays]

loadplot = reshape(load,nrtmpoints*ndays);
rtmeprplot = reshape(rtmepr,nrtmpoints*ndays);
dameprplot = reshape(damepr,ndampoints*ndays);
damregupprplot = reshape(damreguppr,ndampoints*ndays);
damregdownprplot = reshape(damregdownpr,ndampoints*ndays);
ebatplot = reshape(convertToArray3(ebat,n3),nrtmpoints*ndays);
socplot = reshape(convertToArray3(soc,n3),nrtmpoints*ndays);
Prtmplot = reshape(convertToArray3(Prtm,n3),nrtmpoints*ndays);
Pdamplot = reshape(convertToArray2(Pdam,n2),ndampoints*ndays);
regupdamplot = reshape(convertToArray2(regupdam,n2),ndampoints*ndays);
regdowndamplot = reshape(convertToArray2(regdowndam,n2),ndampoints*ndays);
rtmeprofitplot = reshape(convertToArray3(profitErtm,n3),nrtmpoints*ndays);
dameprofitplot = reshape(convertToArray2(profitEdam,n2),ndampoints*ndays);
damregprofitplot = reshape(convertToArray2(profitregupdam,n2),ndampoints*ndays) + reshape(convertToArray2(profitregdowndam,n2),ndampoints*ndays);


totalpower = zeros(nrtm,ndam,ndays);
totalregup = zeros(nrtm,ndam,ndays);
totalregdown = zeros(nrtm,ndam,ndays);
upperband = zeros(nrtm,ndam,ndays);
lowerband = zeros(nrtm,ndam,ndays);
totalpowerplot = zeros(nrtm,ndam,ndays);
upperbandplot = zeros(nrtm,ndam,ndays);
lowerbandplot = zeros(nrtm,ndam,ndays);

Prtmarray = convertToArray3(Prtm,n3);
Pdamarray = convertToArray2(Pdam,n2);
regupdamarray = convertToArray2(regupdam,n2);
regdowndamarray = convertToArray2(regdowndam,n2);


for i in rtm
        for k in dam
            for l in day
                println("i= ",i)
                println("k= ",k)
                totalpower[i,k,l] = Prtmarray[i,k,l] + Pdamarray[k,l];
                totalregup[i,k,l] = regupdamarray[k,l];
                totalregdown[i,k,l] = regdowndamarray[k,l];
                upperband[i,k,l] = totalpower[i,k,l] + totalregup[i,k,l];
                lowerband[i,k,l] = totalpower[i,k,l] - totalregdown[i,k,l];
            end
        end
end


totalpowerplot = reshape(totalpower,nrtmpoints*ndays);    
upperbandplot = reshape(upperband,nrtmpoints*ndays);
lowerbandplot = reshape(lowerband,nrtmpoints*ndays);


    
#################### Plotting #####################

timertm = 0:dtrtm:ndam*ndays;
timedam = 0:dtdam:ndam*ndays;

end

#=

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
totalpowerplot = push!(totalpowerplot,totalpowerplot[end]);
upperbandplot = push!(upperbandplot,upperbandplot[end]);
lowerbandplot = push!(lowerbandplot,lowerbandplot[end]);




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
plot(timedam,Pdamplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
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
subplot(211)
plot(timertm,ebatplot[1:length(ebatplot)],drawstyle="steps-post",label="Energy of battery (MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("Energy (MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("Energy of Battery (MWh)")

subplot(212)
plot(timertm,socplot,drawstyle="steps-post",label="SOC",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("SOC (%)",size=20)
grid()
tick_params(labelsize=18)
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")


figure()
plot(timertm,totalpowerplot,drawstyle="steps-post",color="blue",label="Total Power (MW)",LineWidth=1.5)
hold(true)
plot(timertm,lowerbandplot,drawstyle="steps-post",color="red",label="Lower band (MW)",LineWidth=1.5)
plot(timertm,upperbandplot,drawstyle="steps-post",color="green",label="Upper band (MW)",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("Power (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")

=#
