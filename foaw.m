function [b_out, u_buf] = foaw(u, u_buf_old, UM, T)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This is an implementation of the best-fit-FOAW (First-Order-Adaptive-Windowing), proposed by the following paper:
	% - Janabi-Sharifi, Hayward, Chen: "Discrete-time Adaptive Windowing for Velocity Estimation", IEEE-TCST, 2000.
	%    DOI: 10.1109/87.880606
	% It is using slightly different expressions but is analytically identical to Janabi-Sharifi's original work.
	% This implementation is done by Ryo Kikuuwe solely based on the idea presented in the above original work.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% u         : current input
	%%% u_buf_old : the input buffer with a certain length from the previous timestep.
	%%% UM        : Error tolerance
	%%% T         : Sampling interval
	%%% b_out     : the best-fit slope, which is supposed to be the derivative of u.
	%%% u_buf     : the input buffer that should be used in the next timestep.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	NN     = length(u_buf_old);
	u_buf  = [u,u_buf_old(1,1:NN-1)] ;
	b_out  = 0;
	for i=2:NN
		sum_u  = 0;
		sum_tu = 0;
		for k=1:i
			t = -(k-1)*T;
			u =  u_buf(k);  
			sum_u  = sum_u  + u  ;
			sum_tu = sum_tu + t*u;
		end
		sum_t  = -   T*i*(i-1)/2         ; %%% Sum_{k=0}^i (-(k-1)*T)   
		sum_tt =   T*T*i*(i-1)*(2*i-1)/6 ; %%% Sum_{k=0}^i (-(k-1)*T)^2 
		tb     = sum_t/i;
		ub     = sum_u/i;
		Vtt    = sum_tt-sum_t*tb;
		Vtu    = sum_tu-sum_t*ub;
		b      = Vtu/Vtt   ;
		a      = ub - b*tb ;
		out_of_band = 0;
		for k=1:i
			if abs(a+b*(-(k-1)*T)-u_buf(k)) > UM 
				out_of_band = 1; 
				break
			end
		end
		if out_of_band == 1
			break
		end
		b_out = b; 
	end
end
