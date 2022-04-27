function u_ex=exact_solution(x,k)
if (k==1)
    u_ex=-1/12*x^4+37/12;
end

if(k==2)
    u_ex=x^3-x^4+3;
end