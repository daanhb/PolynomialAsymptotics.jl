
binDig = 88;
set_bigfloat_precision(binDig);
# big or BigFloat can be used on 2.0 and integers because they are exactly representable, and big(pi) should also (about) have all digits.
relacc = BigFloat(2.0)^(-binDig*big(4)/5);


for wei = 1:4

if wei == 1
	function Q(z) big(7)/10*z^3+big(3)/2; end
	alpha = big(28)/10;
	maxN = 128;
	zt = Complex{BigFloat}(1,-1)/big(100);
        function betan(n) big(n)^(1/3)*(big(3)/2*big(+7)/10*gamma(big(3.5))/gamma(big(4))/sqrt(big(pi)))^(-1/big(3)); end;
	label = "Mon0p7ThirdDegImag"
	begIn = BigFloat(0.0);

elseif wei == 2
	function Q(z) BigFloat(1.0)*z^4 -BigFloat(3)*z^2; end
	alpha = BigFloat(0.0);
	zt = big(31)/100;
	maxN = 128*2+1
        function betan(n) big(n)^(1/2)*(big(2)/2*big(+1.0)*gamma(big(2.5))/gamma(big(3))/sqrt(big(pi)))^(-1/big(2)); end
	label = "Hermite_x^4"
	begIn = BigFloat(-Inf);

elseif wei == 3
	function Q(z) BigFloat(exp(z)); end
	alpha = BigFloat(-1/2);
	zt = BigFloat(6)/10;
	maxN = 128;
        function betan(n) big(log(4*big(pi)*n) - log(log(4*big(pi)*n))); end
	label = "Exp"
	begIn = BigFloat(0);

elseif wei == 4
	function Q(z) big(1)*z^6 -big(21)/10*z^5 + big(3)*z^3 -big(6)*z^2 +big(9); end
	alpha = big(11)/10;
	maxN = 128;    
	zt = Complex{BigFloat}(-1,-2);
        function betan(n) big(n)^(1/6)*(big(6)/2*big(1)*gamma(big(6.5))/gamma(big(7))/sqrt(big(pi)))^(-1/big(6)); end;
	label = "genPoly"
	begIn = BigFloat(0);
end

printtoc = 30; 
anPs = zeros(BigFloat,maxN+2,1);
bnm1Ps = zeros(BigFloat,maxN+2,1);
bnm1Ps[1] = sqrt(quadgk(x -> x^alpha*exp(-Q(x)), begIn, BigFloat(+Inf), reltol=relacc)[1]);


function pn(n,z) 
        pnm2 = 0;
        if n < 0
            return 0*z;
	end
        pnm1 = 1/bnm1Ps[1];
        for j = 1:n
            pnew = ((z-anPs[j])*pnm1 -bnm1Ps[j]*pnm2 )/bnm1Ps[j+1];
            pnm2 = pnm1;
            pnm1 = pnew;
	end
        return pnm1;
end

function saveAll()
# open(string("/home/path/bnm1Ps", label, ".txt"), "w") do g # to specify another path
open(string("bnm1Ps", label, ".txt"), "w") do g
        write(g, string(bnm1Ps));
end
open(string("anPs", label, ".txt"), "w") do g
        write(g, string(anPs));
end
# writedlm(string("/bnm1Ps", label, ".txt"), bnm1Ps); # Also possible

myVal = zeros(Complex{Float64}, convert(Int64,floor(log(maxN)/log(2))),1);
for ni = 1:length(myVal) 
        myVal[ni] = pn(2^ni, zt*betan(2^ni));
end

open(string("result", label, ".txt"), "w") do g
        write(g, string(label, ": evaluated polys are: ", myVal) );
end
if wei == 3
        for ni = 1:length(myVal)
                myVal[ni] = bnm1Ps[2^ni+1];
        end
        open(string("coefs", label, ".txt"), "w") do g
                write(g, string(label, ": Recurrence coefficients are: ", myVal) );
        end
end

end


start = time();
ticc = start;

for nt = 1:(maxN+1)
	ea = big(0);
        if wei == 2 # We know the integrand is odd so the integral is zero, to speed up the computations.
		anPs[nt] = big(0);
	else
		(anPs[nt], ea) = quadgk(x-> x^alpha*exp(-Q(x))*x*(pn(nt-1,x))^2, begIn, BigFloat(+Inf), reltol=relacc, abstol = 2^(-binDig/6*5)); # Add abstol for when values are zero (ex. Hermite with odd integrand)
	end
        (tb, eb) = quadgk(x -> x^alpha*exp(-Q(x))*x*pn(nt-1,x)*((x-anPs[nt])*pn(nt-1,x) -bnm1Ps[nt]*pn(nt-2,x)), begIn, BigFloat(+Inf), reltol=relacc, abstol = 2^(-binDig/6*5));
	bnm1Ps[nt+1] = sqrt(tb);
	if (abs(ea) > max(abs(anPs[nt]*2*relacc), 2^(-binDig/6*5) ) ) || (abs(eb) > max(abs(tb*2*relacc), 2^(-binDig/6*5) ) )
		warn(ea, ", ", anPs[nt], " Integrals not computed accurately: ", eb, " ", tb);
	end
        if (time() > ticc + printtoc) || (nt == 10)
            ticc = time();
#            print(wei, " =wei,nt= ", nt, " ", now(), " ", Float64(anPs[nt]), " ", Float64(bnm1Ps[nt+1]), " ", Libc.strftime((maxN^2-nt^2)/nt^2*(ticc-start)+ticc), "\n"); # Print time info that is easier to read if Libc and now() are available.
            saveAll();
            print(wei, " =wei,nt= ", nt, " ", ticc, " ", convert(Float64,anPs[nt]), " ", convert(Float64,bnm1Ps[nt+1]), " ", (maxN^2-nt^2)*(ticc-start)/nt^2, "\n");
	end
end


saveAll();
print("Ended weight ", wei, ".\n")

end # loop over nonstd weights

# std Hermite:
myVal= zeros(BigFloat, 9,1)
for ni = 1:9
        n = 2^ni;
	x = sqrt(388/big(100)*n);
        pnm2 = big(0);
        pnm1 = big(pi)^(-1/4);
        for j = 1:(2*n)
            pnew = 2*x*pnm1/sqrt(big(2)*j) -(j-1)*pnm2/sqrt(j*max(1,big(j-1)));
            pnm2 = pnm1;
            pnm1 = pnew;
	end
        myVal[ni] = pnm1;
end
open(string("resultStdNormHermite.txt"), "w") do g
        write(g, string("Std normalized Hermite: evaluated polys are: ", myVal) );

end


