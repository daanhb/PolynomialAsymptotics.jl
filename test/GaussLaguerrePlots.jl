# Initialising
using PyPlot;
using FastGaussQuadrature; # For besselroots.jl

include("gausslaguerreTest.jl"); 


# Test for accuracy of the weights
nt = 200; alphat = 0.7; (xd,wd) = gausslaguerreTest(nt,alphat,"rec");
for hT = 1:4; (xex,wex) = gausslaguerreTest(nt, alphat, "exp", 1.0, 1, 2*hT); semilogy(xex, abs(wd-wex)./wd, label="$hT terms"); end; 
legend(loc="upper center"); xlabel(L"$x_k$"); ylabel(L"Relative error on $w_k$");


# Timings 
ns = Array{Int64,1}(round(10.^((6:24)/3))); 
lns = length(ns);
rests = zeros(5,lns); alphat = 0.0; nt = 100;

figure; 
maxTim = 10; # Stop if the previous gausslaguerreTest took over 10 seconds.

(x,w) = gausslaguerreTest(nt,alphat,"GW");
for ni = 1:lns
	rests[2,ni] = @elapsed gausslaguerreTest(ns[ni],alphat,"GW");
	print("\n", ni, " GW ", rests[2,ni]);
	if(rests[2,ni] > maxTim)	break;	end
end
loglog(ns, rests[2,:], label="Golub-Welsch");

(x,w) = gausslaguerreTest(nt,alphat,"rec");
for ni = 1:lns
	rests[1,ni] = @elapsed gausslaguerreTest(ns[ni],alphat,"rec");
	print("\n", ni, " rec ", rests[1,ni]);
	if(rests[1,ni] > maxTim)	break;	end
end
loglog(ns, rests[1,:], label="Newton on recurrence");

(x,w) = gausslaguerreTest(nt,alphat,"GLR");
for ni = 1:lns
	rests[3,ni] = @elapsed gausslaguerreTest(ns[ni],alphat,"GLR");
	print("\n", ni, " GLR ", rests[3,ni]);
	if(rests[3,ni] > maxTim)	break;	end
end
loglog(ns, rests[3,:], label="Glaser-Liu-Rokhlin");


(x,w) = gausslaguerreTest(nt,alphat,"gen");
for ni = 1:lns
	rests[4,ni] = @elapsed gausslaguerreTest(ns[ni],alphat,"gen");
	print("\n", ni, " gen ", rests[4,ni]);
	if(rests[4,ni] > maxTim)	break;	end
end
loglog(ns, rests[4,:], label="Newton on a.e.");


(x,w) = gausslaguerreTest(nt,alphat,"exp");
for ni = 1:lns
	rests[5,ni] = @elapsed gausslaguerreTest(ns[ni],alphat,"exp");
	print("\n", ni, " exp ", rests[5,ni]);
	if(rests[5,ni] > maxTim)	break;	end
end
loglog(ns, rests[5,:], label="Asymptotic expansions");

loglog(10.^(3:8), [0.216; 0.311; 0.391; 0.510; 2.71; 23.9], label="Bremer (Fortran)")

legend(loc="best");
xlabel(L"$n$");
ylabel("Time (s)")


