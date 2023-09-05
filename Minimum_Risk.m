function Minimum_Risk(M)

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

%check det(C)
if det(Cov_mat)==0
     quit(0,"force");
end

%result by the CVX:
cvx_begin
    variables w_Long_only(num) ;
    minimize 0.5*transpose(w_Long_only)*Cov_mat*w_Long_only;
    subject to
        transpose(expectancy_vec)*w_Long_only>=M;
        transpose(w_Long_only)*e==1;
        i=1;
        while i<=num
            w_Long_only(i)>=0;
            i = i+1;
        end
cvx_end

sigma_Long_only = transpose(w_Long_only)*Cov_mat*w_Long_only;
mu_Long_only = transpose(expectancy_vec)*w_Long_only;

%result by the algorithm:
inverse_Cov_mat = inv(Cov_mat);
%after computing in the word:
a = transpose(expectancy_vec)*inverse_Cov_mat*expectancy_vec;
b = transpose(expectancy_vec)*inverse_Cov_mat*e;
c = transpose(e)*inverse_Cov_mat*e;
special_matrix = [a b;b c];
inverse_special_mat = (1/(a*c-b^2))*[c -b;-b a];
inverse_special_mat = inv(special_matrix);
lambda_vec = inverse_special_mat*[M;1];
w = lambda_vec(1)*inverse_Cov_mat*expectancy_vec + lambda_vec(2)*inverse_Cov_mat*e;

%print the components of problem:
yield_matrix
expectancy_vec
Cov_mat


k=1;
sum = 0;
while k<=num
    sum = sum + w(k);
    k=k+1;
end
if sum < 0 || sum > 1
    disp('not existed');
else
    w_Short_Long = w
    sigma_Short_Long = transpose(w)*Cov_mat*w
    mu_Short_Long = transpose(expectancy_vec)*w
end
if M==0
    disp('You entered return expectancy equal to zero');
end

%generate weight vector to plot graph:
i=1;
while i<=200
    k=1;
    while k<num
        genterate_w(k,1) = rand();
        k=k+1;
    end
    genterate_w(k,1) = 0;
    genterate_w(k,1) = 1 - ones(1,num)*genterate_w; %the last coordinate complete for one
    sigma_to_plot(i,1) = transpose(genterate_w)*Cov_mat*genterate_w;
    mu_to_plot(i,1) = transpose(expectancy_vec)*genterate_w;
    i=i+1;
end


figure(1)
plot(sigma_to_plot,mu_to_plot,'r.');
xlabel('Variance of the Portfolio');
ylabel('expected return of the Portfolio');
hold on
y = -2*M:0.001:2*M;
if M==0
 y = -0.1:0.001:0.1;   
end
y_square = y.^2;
x = (1/(a*c-b^2))*(y_square.*c-y.*2*b+a);
plot(x,y,'b-');

figure(2)
plot(mu_to_plot,sigma_to_plot,'r.');
y = 0:0.01:M;
x = M*ones(1,size(y,2));
hold on
xlabel('expected return of the Portfolio');
ylabel('Variance of the Portfolio');
plot(x,y,'g-');
hold on
x = -2*M:0.001:2*M;
if M==0
 x = -0.1:0.001:0.1;   
end
x_square = x.^2;
y = (1/(a*c-b^2))*(x_square.*c-x.*2*b+a);
plot(x,y,'b-');


%print Long results:
w_Long_only
mu_Long_only
sigma_Long_only

end
