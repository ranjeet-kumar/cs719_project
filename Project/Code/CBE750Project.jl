using JuMP
using Gurobi
using PyPlot
close("all")

ndays = 365;

############# Funtions to convert JuMP returned variables to arrays ################

function convertToArray(x)
	y = getValue(x)
	n = length(y)
	a = zeros(n)
	for i = 1:n
		a[i] = y[i]
	end
	return a
end

function convertToArray2(A)	
	AA = getValue(A)	
	n = [4,24]
	m = (n[1],n[2])
	B = zeros(m)	
	for i in 1:n[1]
		for j in 1:n[2]
			B[i,j] = AA[i,j]
		end
	end	
	return B
end

function convertToArray3(A)
	AA = getValue(A)	
	n = [3,4,24]
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
qhmpricedata = readcsv("AggregatedData_HASP_ALTA31GT_7_B1.csv")
dampricedata = readcsv("AggregatedData_DAM_ALTA31GT_7_B1.csv")

dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day ahead market (hourly) intervals in a day
dtqhm = 15/60;				#time interval of each qhour [hour]
nqhm = Int64(dtdam/dtqhm);		#No. of quarter hours in an hour
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtqhm/dtrtm);		#No. of real time intervals in each quarter-hour [hour]

nrtmpoints = ndam*nqhm*nrtm;		#Total number of points in RTM data
nqhmpoints = ndam*nqhm;			#Total number of points in quarter hourly data
ndampoints = ndam;			#Total number of points in hourly data	

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
e_max = 1;	          #Maximum power, MW
regup_max = 0.5*e_max;    #Regulation Up Capacity, MW
regdown_max = 0.5*e_max;  #Regulation Up Capacity, MW
rampmin = -100;	          #Lower bound for ramp discharge, MW/s
rampmax = 100;  	  #Upper bound for ramp discharge, MW/s
eff = 0.9;                  #Discharging Efficiency of battery
soc0 = 100;		  #Initial State of charge, 100 means fully charged
socend = 100;		  #State of charge at the end of the day

dailyprofit = zeros(ndays);
loadplot = zeros(ndays*nrtmpoints);
rtmeprplot = zeros(ndays*nrtmpoints);
qhmeprplot = zeros(ndays*nqhmpoints);
dameprplot = zeros(ndays*ndampoints);
rtmregupprplot = zeros(ndays*nrtmpoints);
qhmregupprplot = zeros(ndays*nqhmpoints);
damregupprplot = zeros(ndays*ndampoints);
rtmregdownprplot = zeros(ndays*nrtmpoints);
qhmregdownprplot = zeros(ndays*nqhmpoints);
damregdownprplot = zeros(ndays*ndampoints);
ebatplot = zeros(ndays*nrtmpoints);
socplot = zeros(ndays*nrtmpoints);
ertmplot = zeros(ndays*nrtmpoints);
eqhmplot = zeros(ndays*nqhmpoints);
edamplot = zeros(ndays*ndampoints);
reguprtmplot = zeros(ndays*nrtmpoints);
regupqhmplot = zeros(ndays*nqhmpoints);
regupdamplot = zeros(ndays*ndampoints);
regdownrtmplot = zeros(ndays*nrtmpoints);
regdownqhmplot = zeros(ndays*nqhmpoints);
regdowndamplot = zeros(ndays*ndampoints);
rtmeprofitplot = zeros(ndays*nrtmpoints);
qhmeprofitplot = zeros(ndays*nqhmpoints);
dameprofitplot = zeros(ndays*ndampoints);
rtmregprofitplot = zeros(ndays*nrtmpoints);
qhmregprofitplot = zeros(ndays*nqhmpoints);
damregprofitplot = zeros(ndays*ndampoints);

totalpower = zeros(nrtm,nqhm,ndam,ndays);
totalregup = zeros(nrtm,nqhm,ndam,ndays);
totalregdown = zeros(nrtm,nqhm,ndam,ndays);
upperband = zeros(nrtm,nqhm,ndam,ndays);
lowerband = zeros(nrtm,nqhm,ndam,ndays);
totalpowerplot = zeros(nrtm,nqhm,ndam,ndays);
upperbandplot = zeros(nrtm,nqhm,ndam,ndays);
lowerbandplot = zeros(nrtm,nqhm,ndam,ndays);






j=1;
for p in 1:ndays

#Load and price data
load = loaddata[2:nrtmpoints+1,2+p]/1000;	#Load, MW

eprrtm = rtmpricedata[(p-1)*nrtmpoints+[1:nrtmpoints],4];	    	#Real Time Market price, $/MWh    
eprqhm = qhmpricedata[(p-1)*nqhmpoints+[1:nqhmpoints],4];	    	#Quarter Hourly Market Selling price, $/MWh    
eprdam = dampricedata[(p-1)*ndampoints+[1:ndampoints],4];	    	#Day Ahead Market Selling price, $/MWh    
regupprrtm = rtmpricedata[(p-1)*nrtmpoints+[1:nrtmpoints],5];	    	#Real Time Market Regulation up price, $/MWh
regdownprrtm = rtmpricedata[(p-1)*nrtmpoints+[1:nrtmpoints],6]; 		#Real Time Market Regulation down price, $/MWh
regupprqhm = qhmpricedata[(p-1)*nqhmpoints+[1:nqhmpoints],5];	    	#Quarter Hourly Market Regulation up price, $/MWh
regdownprqhm = qhmpricedata[(p-1)*nqhmpoints+[1:nqhmpoints],6];    	#Quarter Hourly Market Regulation down price, $/MWh
regupprdam = dampricedata[(p-1)*ndampoints+[1:ndampoints],5];	    	#Day Ahead Market Regulation up price, $/MWh
regdownprdam = dampricedata[(p-1)*ndampoints+[1:ndampoints],6]; 		#Day Ahead Market Regulation down price, $/MWh

#Reshape the data to matrices
rtmepr = reshape(eprrtm,nrtm,nqhm,ndam);
qhmepr = reshape(eprqhm,nqhm,ndam);
damepr = reshape(eprdam,ndam);
damreguppr = reshape(regupprdam,ndam);
damregdownpr = reshape(regdownprdam,ndam);
rtmreguppr = reshape(regupprrtm,nrtm,nqhm,ndam);
rtmregdownpr = reshape(regdownprrtm,nrtm,nqhm,ndam);
qhmreguppr = reshape(regupprqhm,nqhm,ndam);
qhmregdownpr = reshape(regdownprqhm,nqhm,ndam);
load = reshape(load,nrtm,nqhm,ndam);

#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
qhm = 1:nqhm;
dam = 1:ndam;

################ Model ##################

m = Model(solver = GurobiSolver())

    @defVar(m, -e_max <= ertm[rtm,qhm,dam] <= e_max)	                #Net Power sold to the real time market, kW
    @defVar(m, -e_max <= eqhm[qhm,dam] <= e_max)	                #Net Power sold to the real time market, kW
    @defVar(m, -e_max <= edam[dam] <= e_max)    	                #Net Power sold to the real time market, kW
    @defVar(m, 0 <= ebat[rtm,qhm,dam] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval
    @defVar(m, 0 <= soc[rtm,qhm,dam] <= 100)		        #SOC of the battery at the end of each real time interval      
    @defVar(m, 0 <= reguprtm[rtm,qhm,dam] <= regup_max)         #Amount of regulation up, kW
    @defVar(m, 0 <= regdownrtm[rtm,qhm,dam] <= regdown_max)     #Amount of regulation down, kW
    @defVar(m, 0 <= regupqhm[qhm,dam] <= regup_max)             #Amount of regulation up, kW
    @defVar(m, 0 <= regdownqhm[qhm,dam] <= regdown_max)         #Amount of regulation down, kW
    @defVar(m, 0 <= regupdam[dam] <= regup_max)                 #Amount of regulation up, kW
    @defVar(m, 0 <= regdowndam[dam] <= regdown_max)             #Amount of regulation down, kW
    @defVar(m, profitertm[rtm,qhm,dam])# >= 0)				        #Profit from the real time market, USD
    @defVar(m, profiteqhm[qhm,dam])# >= 0)  				#Profit from the quarter hourly market, USD
    @defVar(m, profitedam[dam])# >= 0)	        			#Profit from the day ahead market, USD		
    @defVar(m, profitreguprtm[rtm,qhm,dam])# >= 0)	        		#Profit from the real time market, USD
    @defVar(m, profitregupqhm[qhm,dam])# >= 0)		        	#Profit from the quarter hourly market, USD
    @defVar(m, profitregupdam[dam])# >= 0)			        #Profit from the day ahead market, USD
    @defVar(m, profitregdownrtm[rtm,qhm,dam])# >= 0)			        #Profit from the real time market, USD
    @defVar(m, profitregdownqhm[qhm,dam])# >= 0)   			#Profit from the quarter hourly market, USD
    @defVar(m, profitregdowndam[dam])# >= 0)	        		#Profit from the day ahead market, USD
    @defVar(m, profittotal)# >= 0)		        	#Total profit in the day, USD		

    @addConstraint(m, InitialEnergy, ebat[1,1,1] == soc0/100*ebat_max - (eff*(ertm[1,1,1] + eqhm[1,1] + edam[1]) - load[1,1,1])*dtrtm)	#Inital energy in the battery
        
    @addConstraint(m, DefSOC[i=rtm,j=qhm,k=dam], soc[i,j,k] == ebat[i,j,k]/ebat_max*100)			#Define SOC

#    @addConstraint(m, EndSOC[i=rtm[end],j=qhm[end],k=dam[end]], soc[i,j,k] >= socend)		#Constraint on SOC at the end of the day
    
    @addConstraint(m, rtmEBalance[i=rtm[2:end],j=qhm,k=dam], ebat[i,j,k] == ebat[i-1,j,k] - (eff*(ertm[i,j,k] + eqhm[j,k] + edam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint
    
    @addConstraint(m, qhmEBalance[i=rtm[1],j=qhm[2:end],k=dam,iend=rtm[end]], ebat[i,j,k] == ebat[iend,j-1,k] - (eff*(ertm[i,j,k] + eqhm[j,k] + edam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint
    
    @addConstraint(m, damEBalance[i=rtm[1],j=qhm[1],k=dam[2:end],iend=rtm[end],jend=qhm[end]], ebat[i,j,k] == ebat[iend,jend,k-1] - (eff*(ertm[i,j,k] + eqhm[j,k] + edam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint

   # @addConstraint(m, RTMRamp[i=rtm[2:end],j=qhm,k=dam], rampmin*dtrtm <= ertm[i,j,k]  - ertm[i-1,j,k] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

   # @addConstraint(m, QHMRamp[i=rtm[1],j=qhm[2:end],k=dam,iend=rtm[end]], rampmin*dtrtm <= (eqhm[i,j,k] + eqhm[j,k]) - (ebat[iend,j-1,k] + eqhm[j-1,k]) <= rampmax*dtrtm)   #Ramp charging rate constraint at each time
    
   # @addConstraint(m, DAMRamp[i=rtm[1],j=qhm[1],k=dam[2:end],iend=rtm[end],jend=qhm[end]], rampmin*dtrtm <= (ertm[i,j,k] + eqhm[j,k] + edam[k]) - (ertm[iend,jend,k-1] + eqhm[jend,k-1] + edam[k-1]) <= rampmax*dtrtm)   #Ramp discharge rate constraint at each time
    

    @addConstraint(m, rtmRegUp[i=rtm,j=qhm,k=dam], ertm[i,j,k] + eqhm[j,k] + edam[k] + reguprtm[i,j,k] + regupqhm[j,k] + regupdam[k] <= e_max)	#Constraint on total power

    @addConstraint(m, rtmRegDown[i=rtm,j=qhm,k=dam], ertm[i,j,k] + eqhm[j,k] + edam[k] - regdownrtm[i,j,k] - regdownqhm[j,k] - regdowndam[k] >= -e_max)	#Constraint on total power

    @addConstraint(m, qhmRegDown[j=qhm,k=dam], eqhm[j,k] + edam[k] - regdownqhm[j,k] - regdowndam[k] >= -e_max)	#Constraint on total power

    @addConstraint(m, dammRegDown[k=dam], edam[k] - regdowndam[k] >= -e_max)	#Constraint on total power

    @addConstraint(m, RTMEProfits[i=rtm,j=qhm,k=dam], profitertm[i,j,k] == rtmepr[i,j,k]*ertm[i,j,k]*dtrtm)	#Economic calculation    
    @addConstraint(m, QHMEProfits[j=qhm,k=dam], profiteqhm[j,k] == qhmepr[j,k]*eqhm[j,k]*dtqhm)	#Economic calculation
    @addConstraint(m, DAMEProfits[k=dam], profitedam[k] == damepr[k]*edam[k]*dtdam)	#Economic calculation
  
    @addConstraint(m, RTMregupProfits[i=rtm,j=qhm,k=dam], profitreguprtm[i,j,k] == rtmreguppr[i,j,k]*reguprtm[i,j,k])
    @addConstraint(m, QHMregupProfits[j=qhm,k=dam], profitregupqhm[j,k] == qhmreguppr[j,k]*regupqhm[j,k])
    @addConstraint(m, DAMregupProfits[k=dam], profitregupdam[k] == damreguppr[k]*regupdam[k])

    @addConstraint(m, RTMregdownProfits[i=rtm,j=qhm,k=dam], profitregdownrtm[i,j,k] == rtmregdownpr[i,j,k]*regdownrtm[i,j,k])
    @addConstraint(m, QHMregdownProfits[j=qhm,k=dam], profitregdownqhm[j,k] == qhmregdownpr[j,k]*regdownqhm[j,k])
    @addConstraint(m, DAMregdownProfits[k=dam], profitregdowndam[k] == damregdownpr[k]*regdowndam[k])

    @addConstraint(m, TotalProfit, profittotal == sum(profitertm) + sum(profiteqhm) + sum(profitedam) + sum(profitreguprtm) + sum(profitregupqhm) + sum(profitregupdam) + sum(profitregdownrtm) + sum(profitregdownqhm) + sum(profitregdowndam))
    
    @setObjective(m, Max, profittotal);
    
#    print(m)

    status = solve(m)

    soc0 = getValue(soc[nrtm,nqhm,ndam]);

    println("\nTotal Profits on day ", p, ": ", getValue(profittotal),"\n")
    dailyprofit[j] = getValue(profittotal);


loadplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(load,nrtmpoints);
rtmeprplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(rtmepr,nrtmpoints);
qhmeprplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(qhmepr,nqhmpoints);
dameprplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(damepr,ndampoints);
rtmregupprplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(rtmreguppr,nrtmpoints);
qhmregupprplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(qhmreguppr,nqhmpoints);
damregupprplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(damreguppr,ndampoints);
rtmregdownprplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(rtmregdownpr,nrtmpoints);
qhmregdownprplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(qhmregdownpr,nqhmpoints);
damregdownprplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(damregdownpr,ndampoints);
ebatplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(ebat),nrtmpoints);
socplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(soc),nrtmpoints);
ertmplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(ertm),nrtmpoints);
eqhmplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(convertToArray2(eqhm),nqhmpoints);
edamplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(convertToArray(edam),ndampoints);
reguprtmplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(reguprtm),nrtmpoints);
regupqhmplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(convertToArray2(regupqhm),nqhmpoints);
regupdamplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(convertToArray(regupdam),ndampoints);
regdownrtmplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(regdownrtm),nrtmpoints);
regdownqhmplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(convertToArray2(regdownqhm),nqhmpoints);
regdowndamplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(convertToArray(regdowndam),ndampoints);
rtmeprofitplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(profitertm),nrtmpoints);
qhmeprofitplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(convertToArray2(profiteqhm),nqhmpoints);
dameprofitplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(convertToArray(profitedam),ndampoints);
rtmregprofitplot[(j-1)*nrtmpoints+1:j*nrtmpoints] = reshape(convertToArray3(profitreguprtm),nrtmpoints) + reshape(convertToArray3(profitregdownrtm),nrtmpoints);
qhmregprofitplot[(j-1)*nqhmpoints+1:j*nqhmpoints] = reshape(convertToArray2(profitregupqhm),nqhmpoints) + reshape(convertToArray2(profitregdownqhm),nqhmpoints);
damregprofitplot[(j-1)*ndampoints+1:j*ndampoints] = reshape(convertToArray(profitregupdam),ndampoints) + reshape(convertToArray(profitregdowndam),ndampoints);

n4 = [nrtm,nqhm,ndam,ndays]
n3 = [nqhm,ndam,ndays]
n2 = [ndam,ndays]

ertmarray = convertToArray3(ertm);
eqhmarray = convertToArray2(eqhm);
edamarray = convertToArray(edam);
reguprtmarray = convertToArray3(reguprtm);
regupqhmarray = convertToArray2(regupqhm);
regupdamarray = convertToArray(regupdam);
regdownrtmarray = convertToArray3(regdownrtm);
regdownqhmarray = convertToArray2(regdownqhm);
regdowndamarray = convertToArray(regdowndam);


for i in rtm
    for l in qhm
        for k in dam
                totalpower[i,l,k,j] = ertmarray[i,l,k] + eqhmarray[l,k] + edamarray[k];
                totalregup[i,l,k,j] = reguprtmarray[i,l,k] + regupqhmarray[l,k] + regupdamarray[k];
                totalregdown[i,l,k,j] = regdownrtmarray[i,l,k] + regdownqhmarray[l,k] + regdowndamarray[k];
                upperband[i,l,k,j] = totalpower[i,l,k] + totalregup[i,l,k];
                lowerband[i,l,k,j] = totalpower[i,l,k] - totalregdown[i,l,k];
        end
    end
end


    j = j+1;

end

println("Total Profits: ", sum(dailyprofit),"\n")



totalpowerplot = reshape(totalpower,nrtmpoints*ndays);    
upperbandplot = reshape(upperband,nrtmpoints*ndays);
lowerbandplot = reshape(lowerband,nrtmpoints*ndays);





#################### Plotting #####################

timertm = 0:dtrtm:ndam*ndays;
timeqhm = 0:dtqhm:ndam*ndays;
timedam = 0:dtdam:ndam*ndays;

loadplot = push!(loadplot,loadplot[end]);
rtmeprplot = push!(rtmeprplot,rtmeprplot[end]);
qhmeprplot = push!(qhmeprplot,qhmeprplot[end]);
dameprplot = push!(dameprplot,dameprplot[end]);
rtmregupprplot = push!(rtmregupprplot,rtmregupprplot[end]);
qhmregupprplot = push!(qhmregupprplot,qhmregupprplot[end]);
damregupprplot = push!(damregupprplot,damregupprplot[end]);
rtmregdownprplot = push!(rtmregdownprplot,rtmregdownprplot[end]);
qhmregdownprplot = push!(qhmregdownprplot,qhmregdownprplot[end]);
damregdownprplot = push!(damregdownprplot,damregdownprplot[end]);
#ebatplot = push!(ebatplot,ebatplot[end]);
ebatplot = [ebat_max,ebatplot];
#socplot = push!(socplot,socplot[end]);
socplot = [soc0,socplot];
ertmplot = push!(ertmplot,ertmplot[end]);
eqhmplot = push!(eqhmplot,eqhmplot[end]);
edamplot = push!(edamplot,edamplot[end]);
reguprtmplot = push!(reguprtmplot,reguprtmplot[end]);
regupqhmplot = push!(regupqhmplot,regupqhmplot[end]);
regupdamplot = push!(regupdamplot,regupdamplot[end]);
regdownrtmplot = push!(regdownrtmplot,regdownrtmplot[end]);
regdownqhmplot = push!(regdownqhmplot,regdownqhmplot[end]);
regdowndamplot = push!(regdowndamplot,regdowndamplot[end]);
rtmeprofitplot = push!(rtmeprofitplot,rtmeprofitplot[end]);
qhmeprofitplot = push!(qhmeprofitplot,qhmeprofitplot[end]);
dameprofitplot = push!(dameprofitplot,dameprofitplot[end]);
rtmregprofitplot = push!(rtmregprofitplot,rtmregprofitplot[end]);
qhmregprofitplot = push!(qhmregprofitplot,qhmregprofitplot[end]);
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

function qhmtortm(A) 
    A = A[1:nqhmpoints*ndays]
    B = zeros(nrtmpoints*ndays);
    for i in 1:nqhmpoints*ndays
        B[(i-1)*3+1:i*3] = repmat([A[i]],3);
    end
    B = push!(B,B[end])
    return B
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
qhmeprofitplot = cumul(qhmeprofitplot);
dameprofitplot = cumul(dameprofitplot);
totaleprofitplot = rtmeprofitplot + qhmtortm(qhmeprofitplot) + damtortm(dameprofitplot);
rtmregprofitplot = cumul(rtmregprofitplot);
qhmregprofitplot = cumul(qhmregprofitplot);
damregprofitplot = cumul(damregprofitplot);
totalregprofitplot = rtmregprofitplot + qhmtortm(qhmregprofitplot) + damtortm(damregprofitplot);
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
plot(timertm[(1:nrtmpoints+1)],ertmplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timeqhm[(1:nqhmpoints+1)],qhmeprplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="blue",label="QHM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("QHM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("QHM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timeqhm[(1:nqhmpoints+1)],eqhmplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in QHM",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("QHM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in QHM")
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
plot(timedam[(1:ndampoints+1)],edamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timeqhm[(1:nqhmpoints+1)],qhmregupprplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timeqhm[(1:nqhmpoints+1)],qhmregdownprplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("QHM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("QHM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timeqhm[(1:nqhmpoints+1)],regupqhmplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in QHM",LineWidth=1.5)
hold(true)
plot(timeqhm[(1:nqhmpoints+1)],regdownqhmplot[91*nqhmpoints+(1:nqhmpoints+1)],drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in QHM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("QHM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in QHM")
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
ylim(-1.1*e_max,1.1*e_max)
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
plot(timeqhm/24,qhmeprofitplot,drawstyle="steps-post",color="green",label="QHM",LineWidth=1.5)
#xlabel("Time (hr)",size=20), ylabel("QHM Revenue (\$)",size=20)
#xlim(0,24)
#grid();tick_params(labelsize=18)
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
plot(timeqhm/24,qhmregprofitplot,drawstyle="steps-post",color="green",label="QHM",LineWidth=1.5)
#xlabel("Time (hr)",size=20), ylabel("QHM Revenue (\$)",size=20)
#xlim(0,24)
#grid();tick_params(labelsize=18)
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
plot(timertm,ertmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timeqhm,qhmeprplot,drawstyle="steps-post",color="blue",label="QHM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("QHM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("QHM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timeqhm,eqhmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in QHM",LineWidth=1.5)
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("QHM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in QHM")
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
ylim(-1.1*e_max,1.1*e_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timeqhm,qhmregupprplot,drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timeqhm,qhmregdownprplot,drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("QHM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("QHM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timeqhm,regupqhmplot,drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in QHM",LineWidth=1.5)
hold(true)
plot(timeqhm,regdownqhmplot,drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in QHM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("QHM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in QHM")
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
plot(timeqhm/24,qhmeprofitplot,drawstyle="steps-post",label="QHM",LineWidth=1.5)
#xlabel("Time (hr)",size=20), ylabel("QHM Revenue (\$)",size=20)
#xlim(0,24)
#grid();tick_params(labelsize=18)
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
plot(timeqhm/24,qhmregprofitplot,drawstyle="steps-post",label="QHM",LineWidth=1.5)
#xlabel("Time (hr)",size=20), ylabel("QHM Revenue (\$)",size=20)
#xlim(0,24)
#grid();tick_params(labelsize=18)
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

=#
