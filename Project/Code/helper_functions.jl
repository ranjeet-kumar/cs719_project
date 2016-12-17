######## Funtions to convert JuMP returned dictionaries to arrays ########
function convertToArray(x)
	y = getvalue(x)
	n = length(y)
	a = zeros(n)
	for i = 1:n
		a[i] = y[i]
	end
	return a
end

function convertToArray2(AA,n)
	m = (n[1],n[2])
	B = zeros(m)
	for i in 1:n[1]
		for j in 1:n[2]
			B[i,j] = AA[i,j]
		end
	end
	return B
end

function convertToArray3(AA,n)
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

function convertToArray4(AA,n)
	m = (n[1],n[2],n[3],n[4])
	B = zeros(m)
	for i in 1:n[1]
		for j in 1:n[2]
			for k in 1:n[3]
				for l in 1:n[4]
					B[i,j,k,l] = AA[i,j,k,l]
				end
			end
		end
	end
	return B
end
######## Funtion to get a cumulative sum vector ########
function cumul(A)
    C = zeros(length(A));
    for i in 1:length(A)
        C[i] = sum(A[1:i]);
    end
    return C;
end
######## Funtion to convert a DAM vector to RTM vector ########
function damtortm(A)
		nrtm = 12;
    B = zeros(length(A)*nrtm);
    for i in 1:length(A)
        B[(i-1)*nrtm+1:i*nrtm] = repmat([A[i]],nrtm);
    end
    return B;
end
####### Generate scenarios using Ledoit-Wolf Covariance Estimator ###########
function generate_weekly_loads_scenarios(load,estimation_weeks,ndays_planning,NS,filename)
# Input load is a 288x365 matrix of daily load data at 5 minutes resolution
# Generates a csv file with NS columns corresponding to each scenario of ...
# ... load profile (in kW) for one week more than the planning ndays
# Also returns a matrix of the same

# Reshaping the data as a matrix of 52 weekly profiles
(nrtmpoints,ndays) = size(load);
weekly_ndays = 7;
reshape_nrows = Int64(weekly_ndays*nrtmpoints);
reshape_ncolumns = Int64(floor(nrtmpoints*ndays/reshape_nrows));
loadvec = vec(load);
load_data = loadvec[1:reshape_nrows*reshape_ncolumns];
load_data_weekly = reshape(load_data,reshape_nrows,reshape_ncolumns);
load_estimation_data = load_data_weekly[:,estimation_weeks];

# Using Ledoit-Wolfe Sample Covariance Estimator
(p,n) = size(load_estimation_data);
load_estimation_mean = mean(load_estimation_data,2);
X = load_estimation_data-repmat(load_estimation_mean,1,n); # Columns are 168*1 random vectors with mean 0 and covariance Sigma
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
nweeks_planning = Int64(ceil((ndays_planning/weekly_ndays)));
loadNSdata = zeros(reshape_nrows*(nweeks_planning+1),NS);
R = chol(Snstar);
for j in 1:nweeks_planning+1
    loadNSdata[(j-1)*reshape_nrows+(1:reshape_nrows),:]  = repmat(load_estimation_mean,1,NS) + R'*randn(p,NS);
end
loadNSdata[loadNSdata.<=0] = minimum(loadvec);
writecsv(filename,loadNSdata)
return loadNSdata;
end
# Generating a random path
function generate_sample_paths(load,NS,filename_paths)
(nrtm,ndam,ndays_planning,N) = size(load);
nhours_planning = ndam*ndays_planning;
paths	= Matrix(nhours_planning,NS);
for k in 1:nhours_planning
  paths[k,:] = (Array{Int64,1}(rand(S,NS)))';
end
	writecsv(filename_paths,paths);
	return paths,loadperm;
end
