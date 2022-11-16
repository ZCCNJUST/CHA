%% We thank Prof. Zhihua Xie of Jiangxi University of Science and Technology for providing MATLAB support.


function ch=cchdm(m)
N=log2(m);
x=hadamard(m);
a=sqrt(m);b=a;
cchdmm=zeros(m);
for i=1:m
    row=reshape(x(i,:),[a,b]);
    num1=0;num2=0;
    for j=1:a-1
        if row(1,j)~= row(1,j+1);
            num1=num1+1;
        end
        if row(j,1)~=row(j+1,1)
            num2=num2+1;
        end
        num(i)=(num1+1)*(num2+1);  
    end
end
 
[B,index]=sort(num);
 
for k=1:m
    cchdmm(k,:)=x(index(k),:);
end
ch=cchdmm;
end
