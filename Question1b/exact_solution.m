function u_ex=exact_solution(x,k)
if (k==1)
    u_ex=-1/12*x^4+1/12*x;
end

if(k==2)
    u_ex=-x^4+x^3;
end