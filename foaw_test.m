% Test program illustrates how to use the function "foaw".
% The input is "u" and the output is "v". 
% "v" is supposed to be the time derivative of "u",
% but "u" is contaminated by noise.
clear;
close all;
T   = 0.002;  %%%% Sampling interval
UM  = 0.5  ;  %%%% Error tolerance
NN  = 10000;  %%%% Length of the simulation.
tseq =zeros(1,NN);
useq_ideal =zeros(1,NN);
vseq_ideal =zeros(1,NN);
useq =zeros(1,NN);
for k=1:NN
	tseq(k)=T*k;
	useq_ideal(k)= sin(T*k)+ sin(3.*T*k);
	useq(k)= useq_ideal(k)+0.5*(rand-0.5);
	vseq_ideal(k)=(useq_ideal(k)-useq_ideal(max(1,k-1)))/T;
end
vseq =zeros(1,length(useq));
buf = zeros(1,50); %%%% You have to set the size of buffer here. The larger the smoother.
for k=1:NN
	[vseq(k),buf]=foaw(useq(k),buf,UM,T); %%%% This is how to use it.
end
plot(tseq,vseq,'b');
hold on;
plot(tseq,vseq_ideal,'m');
hold on;
plot(tseq,useq,'r');
hold on;
plot(tseq,useq_ideal,'c');
hold off;
