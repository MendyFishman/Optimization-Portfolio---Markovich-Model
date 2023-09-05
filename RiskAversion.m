function RiskAversion(K)

data_file = readmatrix("data.xlsx");
Security_mat =  data_file(1:end,2:end);
num = size(Security_mat,2);
Period_division = size(Security_mat,1);

yield_matrix = rand(Period_division-1,num); %create random yeild matrix 
i = 2; %start from the second price to calculate the yield
k = 1; 
while k<=num %create the yeild matrix:
   while i<=Period_division
        yield = (Security_mat(i,k)-Security_mat(i-1,k))/Security_mat(i-1,k); %calculate the yield
        yield_matrix(i-1,k) = yield; %insert the calculated yield to the yield matrix
        i = i + 1;
   end
    i = 2;
    k = k + 1;   
end

%calculate the yields' expectancy:
p = 1/(Period_division-1); %calculate the probability
j = 1;
while j<=num
    expectancy_vec(j,1) = p*ones(1,Period_division-1)*yield_matrix(1:end,j:j);
    j = j + 1;
end

%calculate the yields' variance:
i = 1;
while i<=num    %var(x) = cov(x,x) = sigma(Pi*(Xi-E(Xi))^2)
    E_Xi = expectancy_vec(i,1); %get the expectancy of the i investment
    vec_RA_E = yield_matrix(1:end,i:i)-E_Xi*ones(Period_division-1,1); %(RAi-E(RA)) **bulid vector
    sigma = p*(transpose(vec_RA_E)*vec_RA_E);  %sigma(Pi*(Xi-E(Xi))^2)
    var_vec(i,1) = sigma; %create the variance vector
    i = i+1;
end

%compute variance matrix:
%cov(w1,w2) =
%E[(w1-E(w1))(w2-E(w2))]=[i=1,...,n:sum(pi*(w1i-E(w1)(w2i-E(w2)))

i = 1;
while i<=num 
    j = 1;
    while j<=num  
        k = 1;
        sum = 0;
        while k<=(Period_division-1)
            sum = sum + p*(yield_matrix(k,i)-expectancy_vec(i))*(yield_matrix(k,j)-expectancy_vec(j));
            k = k+1;
        end
        Cov_mat(i,j) = sum;
        j = j+1;
    end
    i = i+1;
end

%create e vector:
i = 1;
while i<=num
    e(i,1)=1;
    i=i+1;
end


t=1;
while t<50
    random_K = 1*rand();
    if t == 49
        random_K = 0;
    end
    K_vec(t,1) = random_K;
    cvx_begin
        variables w_Long_only(num) ;
        maximize transpose(expectancy_vec)*w_Long_only-0.5*random_K*transpose(w_Long_only)*Cov_mat*w_Long_only;
        subject to
            transpose(w_Long_only)*e==1;
            i=1;
            while i<=num
                w_Long_only(i)>=0;
                i = i+1;
            end
    cvx_end
    sigma_Long_only = transpose(w_Long_only)*Cov_mat*w_Long_only;
    mu_Long_only = transpose(expectancy_vec)*w_Long_only;
    sigma_vec(t,1) = sigma_Long_only;
    mu_vec(t,1) = mu_Long_only;
    Z_vec(t,1) = transpose(expectancy_vec)*w_Long_only-0.5*random_K*transpose(w_Long_only)*Cov_mat*w_Long_only;
    t=t+1;
end

figure(1)
plot(K_vec,sigma_vec,'r.');
xlabel('K');
%ylabel('Variance of the Portfolio');
hold on
plot(K_vec,mu_vec,'b*');
legend('Variance of the Portfolio','expected return of the Portfolio')

figure(2)
plot(K_vec,Z_vec,'g.');
xlabel('K');
ylabel('Object Function');

%check det(C)
if det(Cov_mat)==0
     quit(0,"force");
end

%result by the CVX:
cvx_begin
    variables w_Long_only(num) ;
    maximize transpose(expectancy_vec)*w_Long_only-0.5*K*transpose(w_Long_only)*Cov_mat*w_Long_only;
    subject to
        transpose(w_Long_only)*e==1;
        i=1;
        while i<=num
            w_Long_only(i)>=0;
            i = i+1;
        end
cvx_end

sigma_Long_only = transpose(w_Long_only)*Cov_mat*w_Long_only;
mu_Long_only = transpose(expectancy_vec)*w_Long_only;


%print the components of problem:
yield_matrix
expectancy_vec
Cov_mat

%print results by CVX:
w_Long_only
sigma_Long_only
mu_Long_only

end