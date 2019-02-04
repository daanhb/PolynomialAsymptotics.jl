# In Julia, do: include("onTheFlyGH.jl"); compare();
using FastGaussQuadrature
using Printf

function compare()

	m = 6
	exa = sqrt(pi)/1^((m+1)/2)*prod((2:2:m).-1)/2^(m/2) # equals \int_{-\infty}^\infty p^m exp(-[k=1]*p^2) dp for even m
	ns = 2 .^(8:2:18)
	print(exa, ", ns=", ns, "\n")
	rests = NaN*ones(4,length(ns))
	s1 = "Time (s) ASY  "
	s2 = "Time (s) cref{Aherm}??  " #Add and check \cref in .tex
	s3 = "Error ASY  "
	s4 = "Error cref{Aherm}??  "
	for ni = 1:length(ns)
		n = ns[ni]
		rests[1,ni] = @elapsed igh = gh(2*n)
		s1 = string(s1, " & ", @sprintf("%7.1e", rests[1,ni] ))
		rests[2,ni] = @elapsed ig = onTheFly(n, ceil(Int64, 17*sqrt(n*1.0)) )
		s2 = string(s2, " & ", @sprintf("%7.1e", rests[2,ni] ))
		@time ig = onTheFly(n, ceil(Int64,17*sqrt(n*1.0)) )
		@time igh = gh(2*n)
		@test (ig-exa)/exa < 1e-10
		@test (igh-exa)/exa < 1e-10
		# print(n, ", ", (ig-exa)/exa, "=ig, igh=", (igh-exa)/exa, "\n")
		rests[3,ni] = (igh-exa)/exa
		s3 = string(s3, " & ", @sprintf("%7.1e", rests[3,ni] ))
		rests[4,ni] = (ig-exa)/exa
		s4 = string(s4, " & ", @sprintf("%7.1e", rests[4,ni] ))
	end
	print(rests, "\n\n asdf \n", s1, "\n", s2, "\n", s3, "\n", s4, "\n") # Prints the table
end

function onTheFly(n::Int64, mi::Int64)

	alpha = -0.5 # Expect Julia to inline this into the following formulae
	d = 1/(4*n+2*alpha+2) # Could specify all variables to be d::Float64
	s = 0.0
	csn = ceil(Int64,sqrt(n*1.0))
	for k=min(n,mi):-1:(csn+1)
		p = (4*n-4*k+3.0)/(4*n+1)
		t = pi^2*(p-1)^2/16
		for i = 1:6 # Newton iterations for t
			t = t-(p*pi+2*sqrt(t-t^2) -acos(2*t-1))/2*sqrt(t/(1-t))
		end
		x = 0.0
		w = 0.0
		# We aim for extremely high n so only four terms suffice
		if n < 1265
			x = x -d^5/181440 *(1-t)^2/t^2 *(10797500*(1-t)^(-10) - 43122800*(1-t)^(-9) + 66424575*(1-t)^(-8) -48469876*(1-t)^(-7) + 193536*alpha^6 + 16131880*(1-t)^(-6) + 80*(315*alpha^4 - 630*alpha^2 -221)*(1-t)^(-4) - 1727136*(1-t)^(-5) - 967680*alpha^4 - 320*(63*alpha^4 - 126*alpha^2 +43)*(1-t)^(-3)  + 384*(945*alpha^6 - 4620*alpha^4 + 6405*alpha^2 - 1346)*(1-t)^(-2) +1354752*alpha^2   - 23040*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-1) -285696)
			w = w - (1-t)^3/90720/t^3 *d^6 *(43190000*(1-t)^(-12) -204917300*(1-t)^(-11) + 393326325*(1-t)^(-10)  - 386872990*(1-t)^(-9) + 201908326*(1-t)^(-8) +80*(315*alpha^4 - 630*alpha^2 + 53752)*(1-t)^(-6)  - 50986344*(1-t)^(-7) - 320*(189*alpha^4 -378*alpha^2 - 89)*(1-t)^(-5) + 480*(63*alpha^4 - 126*alpha^2 + 43)*(1-t)^(-4)  -384*(315*alpha^6 - 1470*alpha^4 + 1995*alpha^2 - 416)*(1-t)^(-3) + 2304*(21*alpha^6 -105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-2) )
		end
		if n < 22026
			x = x +d^3*(1-t)/t/720*(1600*(1-t)^(-6) - 3815*(1-t)^(-5) + 480*alpha^4 +2814*(1-t)^(-4) - 576*(1-t)^(-3)   - 960*alpha^2 - 48*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-1) -16*(1-t)^(-2) + 224)
			w = w +(1-t)^2/720/t^2*d^4*(8000*(1-t)^(-8) - 24860*(1-t)^(-7) + 27517*(1-t)^(-6) - 12408*(1-t)^(-5) + 1712*(1-t)^(-4) +16*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-2) + 32*(1-t)^(-3))
		end
		x = sqrt(t*(4*n+1) -1/(3*(4*n+1))*(5/4*(1-t)^(-2) - (1-t)^(-1) - 1/4) + x)
		s = s+((-x)^6 + (x)^6)/x*exp(-x^2)*pi*sqrt(t/(1-t))*(1 + (2*t + 3)/6/(4*n+1)^2/(t-1)^3 +w)
	end
	for k=csn:-1:1
		jak = pi*(k-1/2)
		x = 0.0
		w = 0.0
		if n < 259
			x= x + (10644*jak^8 + 60*(887*alpha^2 - 2879)*jak^6 + (125671*alpha^4 -729422*alpha^2 + 1456807)*jak^4 + 3*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 - 2201939)*jak^2 + 2*(107959*alpha^8 - 1146220*alpha^6 + 5095482*alpha^4 -10087180*alpha^2 + 6029959) )*d^8/42525
	   		w = w + (215918*alpha^8 + 53220*jak^8 + 240*(887*alpha^2 - 2879)*jak^6 -2292440*alpha^6 + 3*(125671*alpha^4 - 729422*alpha^2 + 1456807)*jak^4 + 10190964*alpha^4  + 6*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 -2201939)*jak^2 - 20174360*alpha^2 + 12059918)/42525*d^8
		end
		if n < 1265
			x = x + (657*jak^6 +36*jak^4*(73*alpha^2-181) +2*jak^2*(2459*alpha^4 -10750*alpha^2 +14051) + 4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6/2835
			w = w + (1493*alpha^6 + 657*jak^6 + 27*(73*alpha^2 - 181)*jak^4 - 9303*alpha^4  + (2459*alpha^4 -10750*alpha^2 + 14051)*jak^2 + 19887*alpha^2 - 12077)*4/2835*d^6
		end
		if n < 22026
			x = x + (11*jak^4 +3*jak^2*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4/45
			w = w + (46*alpha^4 + 33*jak^4 +6*jak^2*(11*alpha^2 -19) -140*alpha^2 +94)/45*d^4
		end
		x = pi*(k-1/2)*sqrt(1/(4*n+1) +(pi^2*(k-1/2)^2 -3/2)/(3*(4*n+1)^3) + x*d )
		s = s+((-x)^6 + (x)^6)*exp(-x^2)/(x*(n+1/4)*(besselj(-3/2, pi*(k-1/2)) )^2)*(1 + (2*pi^2*(k-1/2)^2 -3/2)/(3*(4*n+1)^2) +w)/2
	end
	return s
end


function gh(n)
	(x,w) = gausshermite(2n)
	sum(w .* x.^6)
end
