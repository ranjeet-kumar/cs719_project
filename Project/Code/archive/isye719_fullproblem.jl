
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
nrtm = Int64(dtdam/dtrtm);		#No. of real time intervals in each hour

nrtmpoints = ndam*nrtm;		        #Total number of points in RTM data in one day
ndampoints = ndam;			#Total number of points in hourly data in one day

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
regup_max = 0.5*P_max;    #Regulation Up Capacity, MW
regdown_max = 0.5*P_max;  #Regulation Up Capacity, MW
rampmin = -100;	          #Lower bound for ramp discharge, MW/5min
rampmax = 100;  	  #Upper bound for ramp discharge, MW/5min
eff = 1;                  #Discharging Efficiency of battery
soc0 = 100;		  #Initial State of charge, 100 means fully charged
socend = 100;		  #State of charge at the end of the day
ndays = 365;

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


#Load and price data
load = loaddata[2:nrtmpoints+1,2+(1:ndays)]/1000;	#Load, MW

eprrtm = rtmpricedata[1:nrtmpoints*ndays,4];	    	#Real Time Market price, $/MWh    
eprdam = dampricedata[1:ndampoints*ndays,4];	    	#Day Ahead Market Selling price, $/MWh    
regupprdam = dampricedata[1:ndampoints*ndays,5];	    	#Day Ahead Market Regulation up price, $/MWh
regdownprdam = dampricedata[1:ndampoints*ndays,6]; 		#Day Ahead Market Regulation down price, $/MWh


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




    @constraint(m, InitialEnergy, ebat[1,1,1] == soc0/100*ebat_max - 1/eff*Pnet[1,1,1]*dtrtm - load[1,1,1]*dtrtm)	#Inital energy in the battery
        
    @constraint(m, DefSOC[i=rtm,k=dam,l=day], soc[i,k,l] == ebat[i,k,l]/ebat_max*100)			#Define SOC

#    @constraint(m, EndSOC[i=rtm[end],k=dam[end],l=day], soc[i,k,l] >= socend)		#Constraint on SOC at the end of the day
    
    @constraint(m, rtmEBalance[i=rtm[2:end],k=dam,l=day], ebat[i,k,l] == ebat[i-1,k,l] - 1/eff*Pnet[i,k,l]*dtrtm - load[i,k,l]*dtrtm)	#Dynamics constraint
    
    @constraint(m, damEBalance[i=rtm[1],k=dam[2:end],iend=rtm[end],l=day], ebat[i,k,l] == ebat[iend,k-1,l] - 1/eff*Pnet[i,k,l]*dtrtm - load[i,k,l]*dtrtm)	#Dynamics constraint

    @constraint(m, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l=day[2:end]], ebat[i,k,l] == ebat[iend,kend,l-1] - 1/eff*Pnet[i,k,l]*dtrtm - load[i,k,l]*dtrtm)	#Dynamics constraint

   # @constraint(m, RTMRamp[i=rtm[2:end],k=dam,l=day], rampmin*dtrtm <= Pnet[i,k,l]  - Pnet[i-1,k,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time    

   # @constraint(m, DAMRamp[i=rtm[1],k=dam[2:end],iend=rtm[end],l=day], rampmin*dtrtm <= Pnet[i,k,l] - Pnet[iend,k-1,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time    

   # @constraint(m, DAYRamp[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l=day[2:end]], rampmin*dtrtm <= Pnet[i,k,l] - Pnet[iend,kend,l-1] <= rampmax*dtrtm)   #Ramp discharge constraint at each time    

    @constraint(m, RegUp[i=rtm,k=dam,l=day], Pnet[i,k,l] + regupdam[k,l] <= P_max)	#Constraint on total power

    @constraint(m, RegDown[i=rtm,k=dam,l=day], Pnet[i,k,l] - regdowndam[k,l] >= -P_max)	#Constraint on total power

    @constraint(m, RTMEProfits[i=rtm,k=dam,l=day], profitErtm[i,k,l] == rtmepr[i,k,l]*Prtm[i,k,l]*dtrtm)	#Economic calculation        
    @constraint(m, DAMEProfits[k=dam,l=day], profitEdam[k,l] == damepr[k,l]*Pdam[k,l]*dtdam)	#Economic calculation
  
    @constraint(m, DAMregupProfits[k=dam,l=day], profitregupdam[k,l] == damreguppr[k,l]*regupdam[k,l])
    @constraint(m, DAMregdownProfits[k=dam,l=day], profitregdowndam[k,l] == damregdownpr[k,l]*regdowndam[k,l])

    @constraint(m, TotalProfit, profittotal == sum(profitErtm) + sum(profitEdam) + sum(profitregupdam) + sum(profitregdowndam))
    
    @objective(m, Min, -profittotal)


#    print(m)

    status = solve(m)

    println("\nTotal Profits on day ", getvalue(profittotal),"\n" )

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
